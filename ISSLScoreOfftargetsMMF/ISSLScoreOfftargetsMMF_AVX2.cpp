#include <fstream>
#include "ISSLScoreOfftargetsMMF.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
using namespace boost::iostreams;
using ScoreRecord = std::array<char,32>;


// Char to binary encoding
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
// Binary to char encoding
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

constexpr size_t SCORE_BUFFER_RECORDS = 1000;

// Scoring methods
otScoreMethod scoreMethod;
bool calcCfd = false;
bool calcMit = false;

struct SelectedBucket {
    uint64_t compactOffsetBytes;
    uint64_t numCandidates;
};

uint64_t sequenceToSignature(const std::string& seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[seq[j]]) << (j * 2);
    }
    return signature;
}

string signatureToSequence(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}

static uint64_t readLittleEndian(const uint8_t *data, unsigned bytes) {
    uint64_t value = 0;
    for (unsigned i = 0; i < bytes; ++i)
        value |= uint64_t(data[i]) << (8 * i);
    return value;
}

static uint64_t readPlanNumber(std::istream &input) {
    uint64_t value = 0;
    if (!(input >> value)) throw std::runtime_error("Truncated shard plan");
    return value;
}

void accumulate_scores(double& totScoreMit, double& totScoreCfd, uint64_t mismatches, uint64_t dist, uint64_t occurrences, uint64_t searchSignature, uint64_t offTargetSignature) {
    if (dist >= 0 && dist <= 4) {
        // Begin calculating MIT score
        if (calcMit) {
            if (dist > 0) {
                totScoreMit += precalculatedMITScores.at(mismatches) * (double)occurrences;
            }
        }

        // Begin calculating CFD score
        if (calcCfd) {
            /** "In other words, for the CFD score, a value of 0
                *      indicates no predicted off-target activity whereas
                *      a value of 1 indicates a perfect match"
                *      John Doench, 2016.
                *      https://www.nature.com/articles/nbt.3437
            */
            double cfdScore = 0;
            if (dist == 0) {
                cfdScore = 1;
            }
            else {
                cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM

                for (size_t pos = 0; pos < 20; pos++) {
                    size_t mask = pos << 4;

                    /** Create the mask to look up the position-identity score
                        *      In Python... c2b is char to bit
                        *       mask = pos << 4
                        *       mask |= c2b[sgRNA[pos]] << 2
                        *       mask |= c2b[revcom(offTaret[pos])]
                        *
                        *      Find identity at `pos` for search signature
                        *      example: find identity in pos=2
                        *       Recall ISSL is inverted, hence:
                        *                   3'-  T  G  C  C  G  A -5'
                        *       start           11 10 01 01 10 00
                        *       3UL << pos*2    00 00 00 11 00 00
                        *       and             00 00 00 01 00 00
                        *       shift           00 00 00 00 01 00
                        */
                    uint64_t searchSigIdentityPos = searchSignature;
                    searchSigIdentityPos &= (3ULL << (pos * 2));
                    searchSigIdentityPos = searchSigIdentityPos >> (pos * 2);
                    searchSigIdentityPos = searchSigIdentityPos << 2;

                    /** Find identity at `pos` for offtarget
                        *      Example: find identity in pos=2
                        *      Recall ISSL is inverted, hence:
                        *                  3'-  T  G  C  C  G  A -5'
                        *      start           11 10 01 01 10 00
                        *      3UL<<pos*2      00 00 00 11 00 00
                        *      and             00 00 00 01 00 00
                        *      shift           00 00 00 00 00 01
                        *      rev comp 3UL    00 00 00 00 00 10 (done below)
                        */
                    uint64_t offtargetIdentityPos = offTargetSignature;
                    offtargetIdentityPos &= (3ULL << (pos * 2));
                    offtargetIdentityPos = offtargetIdentityPos >> (pos * 2);

                    /** Complete the mask
                        *      reverse complement (^3UL) `offtargetIdentityPos` here
                        */
                    mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

                    if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
                        cfdScore *= cfdPosPenalties[mask];
                    }
                }
            }
            totScoreCfd += cfdScore * (double)occurrences;
        }

    }
}

static void storeResultWord(
    char *destination, uint64_t value, unsigned bytes) {
    for (unsigned i = 0; i < bytes; ++i) {
        const unsigned char byte =
            static_cast<unsigned char>(value >> (8 * i));
        std::memcpy(destination + i, &byte, 1);
    }
}

static void storeResultDouble(char *destination, double value) {
    static_assert(sizeof(double) == 8 &&
                  std::numeric_limits<double>::is_iec559,
                  "Score output requires IEEE-754 binary64");
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    storeResultWord(destination, bits, 8);
}

static ScoreRecord makeScoreRecord(
    uint64_t querySignature, uint64_t globalId,
    double mitContribution, double cfdContribution, unsigned idBits) {
    if (idBits != 32 && idBits != 64)
        throw std::runtime_error("Unsupported global ID width");
    if (idBits == 32 && globalId > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("Global ID exceeds selected output width");
    ScoreRecord record{};
    storeResultWord(record.data(), querySignature, 8);
    storeResultWord(record.data() + 8, globalId, idBits / 8);
    storeResultDouble(record.data() + 16, mitContribution);
    storeResultDouble(record.data() + 24, cfdContribution);
    return record;
}

static void flushScoreRecords(
    std::ofstream &output, std::vector<ScoreRecord> &buffer) {
    if (buffer.empty()) return;
    static_assert(sizeof(ScoreRecord) == 32, "Unexpected result padding");
    output.write(reinterpret_cast<const char *>(buffer.data()),
                 static_cast<std::streamsize>(buffer.size() * 32));
    if (!output)
        throw std::runtime_error("Failed writing Mapper results");
    buffer.clear();
}

static void emitContribution(
    std::ofstream &output, std::vector<ScoreRecord> &buffer,
    unsigned idBits, uint64_t searchSignature, uint64_t globalId,
    uint64_t offTargetSignature, uint64_t occurrences,
    uint64_t mismatches, uint64_t dist) {
    if (dist > 4) return;
    double mitContribution = 0.0, cfdContribution = 0.0;
    accumulate_scores(mitContribution, cfdContribution,
        mismatches, dist, occurrences, searchSignature, offTargetSignature);
    buffer.push_back(makeScoreRecord(
        searchSignature, globalId, mitContribution, cfdContribution, idBits));
    if (buffer.size() >= SCORE_BUFFER_RECORDS)
        flushScoreRecords(output, buffer);
}

bool seenOfftargetAlready(uint64_t* offtargetTogglesTail, uint64_t signatureId) {
    /** Prevent assessing the same off-target for multiple slices */
    uint64_t seen = 0;
    uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
    seen = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;

    if (!seen) {
        /** Mark the current off-target as seen for this slice **/
        *ptrOfftargetFlag |= (1ULL << (signatureId % 64));
    }

    return seen;
}

uint64_t LEB128Decode(const uint8_t *ptr, uint32_t &bytesUsed)
{
    uint64_t result = 0;
    int shift = 0;
    uint8_t byte = 0;
    bytesUsed = 0;

    do
    {
        byte = ptr[bytesUsed];
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;
        shift += 7;
        bytesUsed++;
    } while (byte & 0x80);

    return result;
}

uint64_t read40BitValue(const uint8_t *p)
{
    uint64_t v = 0;
    for (int b = 0; b < 5; b++)
        v |= static_cast<uint64_t>(p[b]) << (b * 8);
    return v;
}

int main(int argc, char** argv)
{
    auto startLoading = std::chrono::high_resolution_clock::now();

    if (argc < 9) {
        fprintf(stderr, "Usage: %s [candidates] [query file] [plan] [idBits] [max distance] [score-threshold] [score-method] [output prefix]\n", argv[0]);
        exit(1);
    }

    const string candidatePath = argv[1];
    const string queryPath = argv[2];
    const string planPath = argv[3];
    const uint64_t parsedBits = std::stoull(argv[4]);
    const unsigned idBits = (unsigned)(parsedBits);
    size_t thresholdEnd = 0;

    int maxDist = atoi(argv[5]);
    double threshold = atof(argv[6]);
    const string method = argv[7];

    // No early exiting so 'and', 'or' and 'avg' have the same scoring behaviour: to calculate both mit and cfd
    if (method != "and" && method != "or" && method != "avg" &&
    method != "mit" && method != "cfd")
    throw std::runtime_error("Invalid scoring method");
    calcMit = method != "cfd";
    calcCfd = method != "mit";

    const string prefixID = argv[8];
    const uint64_t recordBytes = 12 + idBits / 8;

    const auto fileBytes = std::filesystem::file_size(candidatePath);
    if (fileBytes > std::numeric_limits<size_t>::max() || fileBytes % recordBytes)
        throw std::runtime_error("Invalid candidate file size");

    const uint64_t candidateFileSize = static_cast<uint64_t>(fileBytes);
    mapped_file_source candidateMapping;
    const uint8_t *candidates = nullptr;

    if (candidateFileSize != 0)
    {
        candidateMapping.open(candidatePath);
        if (!candidateMapping.is_open() || candidateMapping.size() != fileBytes)
            throw std::runtime_error("Cannot map complete candidate file");
        candidates = reinterpret_cast<const uint8_t *>(candidateMapping.data());           
    }

    auto endLoading = std::chrono::high_resolution_clock::now();
    auto startProcessing = std::chrono::high_resolution_clock::now();

    //TODO: rewrite
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    constexpr size_t seqLength = 20;
    constexpr size_t seqLineLength = seqLength + 1;
    std::filesystem::path queryFile(queryPath);
    size_t fileSize = std::filesystem::file_size(queryFile);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    FILE* fp = fopen(queryPath.c_str(), "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    vector<double> querySignatureMitScores(queryCount);
    vector<double> querySignatureCfdScores(queryCount);

    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

    /** Binary encode query sequences */
    #pragma omp parallel
    {
    #pragma omp for
        for (int i = 0; i < queryCount; i++) {
            char* ptr = &queryDataSet[i * seqLineLength];
            uint64_t signature = sequenceToSignature(ptr, 20);
            querySignatures[i] = signature;
        }
    }

    /*
    Shard Plan is structured as whitespace separated unsigned integers.
        shardId sliceId bucketCount
        bucketId compactOffsetBytes candidateCount  // repeated for each bucket
        query count
        queryBucketId                               // repeated for each query-bucket tuple

        The bucket rows describe the compacted ranges in candidates.bin:
            compactOffsetBytes is the byte offset of the bucket's first record
            candidateCount is the number of records in the range for the compactOffsetBytes.
        
        The final queryBucketId list maps each query sequence to one of the listed buckets, with 
        its order matching the query.txt file given to ensure queryBucketIds[searchIdx] selects
        the correct candidate range while scoring querySignatures[searchIdx]
    */

    std::ifstream shardRead(planPath);
    if (!shardRead) 
        throw std::runtime_error("Cannot open shard plan");

    const uint64_t shardId = readPlanNumber(shardRead);
    const uint64_t sliceId = readPlanNumber(shardRead);

    if (shardId != sliceId)
        throw std::runtime_error("Shard/slice identity mismatch");

    const uint64_t bucketCount = readPlanNumber(shardRead);
    std::unordered_map<uint64_t, SelectedBucket> selectedBuckets;
    uint64_t expectedOffset = 0;

    // Read selected bucket rows into a lookup so each query can find its 
    // candidate range in candidates.bin by bucketId.
    // bucket Id -> offset, numCandidates
    for (uint64_t i = 0; i < bucketCount; ++i) {
        const uint64_t bucketId = readPlanNumber(shardRead);
        const uint64_t offset = readPlanNumber(shardRead);
        const uint64_t count = readPlanNumber(shardRead);

        if (offset != expectedOffset || offset > candidateFileSize ||
            count > (candidateFileSize - offset) / recordBytes)
            throw std::runtime_error("Invalid candidate bucket range");

        if (!selectedBuckets.emplace(bucketId, SelectedBucket{offset, count}).second)
            throw std::runtime_error("Duplicate bucket in plan");

        expectedOffset += count * recordBytes;
    }
    if (expectedOffset != candidateFileSize)
    throw std::runtime_error("Plan does not cover candidate file");

    if (readPlanNumber(shardRead) != queryCount)
        throw std::runtime_error("Plan/query count mismatch");

    vector<uint64_t> queryBucketIds(queryCount);
    for (size_t queryIdx = 0; queryIdx < queryCount; ++queryIdx) {
        queryBucketIds[queryIdx] = readPlanNumber(shardRead);
        if (selectedBuckets.find(queryBucketIds[queryIdx]) == selectedBuckets.end())
            throw std::runtime_error("Query references a bucket absent from plan");
    }
    shardRead >> std::ws;
    if (shardRead.bad() || shardRead.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("Unexpected trailing plan data");
    shardRead.close();

    // One temporary file and one bounded buffer per actual OpenMP worker.
    int workerCount = 0;
    vector<string> workerFiles(omp_get_max_threads());
    vector<std::exception_ptr> workerErrors(workerFiles.size());

    /** Begin scoring */
    /** Begin scoring */
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        #pragma omp single
        workerCount = omp_get_num_threads();

        std::ofstream threadOut;
        vector<ScoreRecord> writeBuffer;
        try {
            workerFiles[tid] = "mapper_" + prefixID + "_thread_" +
                            std::to_string(tid) + ".bin";
            threadOut.open(workerFiles[tid], std::ios::binary | std::ios::trunc);
            if (!threadOut)
                throw std::runtime_error("Cannot open Mapper worker output");
            writeBuffer.reserve(SCORE_BUFFER_RECORDS);
        } catch (...) {
            workerErrors[tid] = std::current_exception();
        }

        /** For each candidate guide */
        // TODO: update to openMP > v2 (Use clang compiler)
        #pragma omp for
        for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            if (workerErrors[tid]) continue;
            try {
                auto searchSignature = querySignatures[searchIdx];
                __m256i searchSignatureVec = _mm256_set1_epi64x((int64_t)searchSignature);

                const auto &plan = selectedBuckets.at(queryBucketIds[searchIdx]);
                const uint64_t signaturesInSlice = plan.numCandidates;
                if (signaturesInSlice == 0) continue;
                const uint8_t *candidate = candidates + plan.compactOffsetBytes;
                alignas(32) uint64_t sigBuf[4], idBuf[4], occBuf[4];

                /** For each off-target signature in slice */
                size_t j = 0;
                for (; signaturesInSlice - j >= 4; j +=4) {
                    __m256i offTargetsVec;
                    __m256i signatureIdVec;
                    __m256i occurencesVec;

                    for (int lane = 0; lane < 4; lane++)
                    {
                        sigBuf[lane] = readLittleEndian(candidate, 8);
                        candidate += 8;
                        idBuf[lane] = readLittleEndian(candidate, idBits / 8);
                        candidate += idBits / 8;
                        occBuf[lane] = readLittleEndian(candidate, 4);
                        candidate += 4;
                    }

                    offTargetsVec = _mm256_load_si256((__m256i*)sigBuf);
                    signatureIdVec = _mm256_load_si256((__m256i*)idBuf);
                    occurencesVec  = _mm256_load_si256((__m256i*)occBuf);

                    __m256i xoredSignaturesVec = _mm256_xor_si256(searchSignatureVec, offTargetsVec);
                    __m256i evenBitsVec = _mm256_and_si256(xoredSignaturesVec, _mm256_set1_epi64x(0xAAAAAAAAAAAAAAAAULL));
                    __m256i oddBitsVec = _mm256_and_si256(xoredSignaturesVec, _mm256_set1_epi64x(0x5555555555555555ULL));
                    __m256i mismatchesVec = _mm256_or_si256(_mm256_srli_epi64(evenBitsVec, 1), oddBitsVec);

                    alignas(32) uint64_t mismatchesArr[4];
                    alignas(32) uint64_t signatureIdArr[4];
                    alignas(32) uint64_t occurencesArr[4];

                    _mm256_store_si256((__m256i *)mismatchesArr, mismatchesVec);
                    _mm256_store_si256((__m256i *)signatureIdArr, signatureIdVec);
                    _mm256_store_si256((__m256i *)occurencesArr, occurencesVec);

                    uint64_t distArr[4];
                    for (int lane = 0; lane < 4; lane++)
                    {
                        distArr[lane] = popcount64(mismatchesArr[lane]);
                    }

                    for (int lane = 0; lane < 4; lane++) {
                        uint64_t dist = distArr[lane];
                        uint64_t mismatches = mismatchesArr[lane];

                        uint64_t offTargetSignature;
                        uint64_t signatureId;
                        uint32_t occurrences;

                        offTargetSignature = sigBuf[lane];
                        signatureId = idBuf[lane];
                        occurrences = (uint32_t)occBuf[lane];

                        emitContribution(threadOut, writeBuffer, idBits, searchSignature, signatureId, offTargetSignature, occurrences, mismatches, dist);
                    }
                }
                // Clean-up loop
                for (; j < signaturesInSlice; j++) {
                    uint64_t offTargetSignature = readLittleEndian(candidate, 8);
                    candidate += 8;
                    uint64_t signatureId = readLittleEndian(candidate, idBits / 8);
                    candidate += idBits / 8;
                    uint32_t occurrences = (uint32_t)readLittleEndian(candidate, 4);
                    candidate += 4;

                    uint64_t xoredSignatures = searchSignature ^ offTargetSignature;
                    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                    uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                    uint64_t mismatches = (evenBits >> 1) | oddBits;
                    uint64_t dist = popcount64(mismatches);

                    emitContribution(threadOut, writeBuffer, idBits, searchSignature, signatureId, offTargetSignature, occurrences, mismatches, dist);
                }
            }catch(...)
            {
                workerErrors[tid] = std::current_exception();
            }
        }
        try{
            if (!workerErrors[tid])
                flushScoreRecords(threadOut, writeBuffer);
            if (threadOut.is_open()){
                threadOut.close();
                if (!threadOut)
                    throw std::runtime_error("Failed closing Mapper worker output");
            }

        } catch (...) {
            if (!workerErrors[tid]) workerErrors[tid] = std::current_exception();
        }
    } // All worker streams are closed

    for (int tid = 0; tid < workerCount; ++tid) {
        if (workerErrors[tid]) std::rethrow_exception(workerErrors[tid]);
    }

    // Merge only the files belonging to this invocation's actual worker team.
    const string shardFILEOUT = prefixID + "_shard_" + std::to_string(shardId) + ".bin";
    std::ofstream shardOut(shardFILEOUT, std::ios::binary | std::ios::trunc);

    if (!shardOut) throw std::runtime_error("Cannot open final shard output");

    std::array<char, 64 * 1024> copyBuffer;

    for (int tid = 0; tid < workerCount; ++tid) {
        std::ifstream in(workerFiles[tid], std::ios::binary);

        if (!in) throw std::runtime_error("Cannot open Mapper worker file for merge");

        for (;;) {
            in.read(copyBuffer.data(), static_cast<std::streamsize>(copyBuffer.size()));
            const std::streamsize copied = in.gcount();
            if (copied > 0) {
                shardOut.write(copyBuffer.data(), copied);
                if (!shardOut) throw std::runtime_error("Failed writing merged shard");
            }
            if (in.bad()) throw std::runtime_error("Failed reading Mapper worker file");
            if (in.eof()) break; // Empty worker files are valid.
            if (in.fail()) throw std::runtime_error("Failed reading Mapper worker file");
        }

        in.clear(); // Clear the normal EOF state before checking close.
        in.close();

        if (!in) throw std::runtime_error("Failed closing merged worker input");
        if (std::remove(workerFiles[tid].c_str()) != 0)
            throw std::runtime_error("Failed removing merged worker file");
    }
    shardOut.close();
    if (!shardOut) throw std::runtime_error("Failed closing final shard output");
    return 0;


    /** Print global scores to stdout */
    for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
        auto querySequence = signatureToSequence(querySignatures[searchIdx], 20);
        printf("%s\t", querySequence.c_str());
        if (calcMit)
            printf("%f\t", querySignatureMitScores[searchIdx]);
        else
            printf("-1\t");

        if (calcCfd)
            printf("%f\n", querySignatureCfdScores[searchIdx]);
        else
            printf("-1\n");

    }

}
