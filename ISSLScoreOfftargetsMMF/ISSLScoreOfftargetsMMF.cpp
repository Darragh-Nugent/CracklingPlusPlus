#include "ISSLScoreOfftargetsMMF.hpp"
#include <immintrin.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;
using namespace boost::iostreams;

// Char to binary encoding
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
// Binary to char encoding
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

// Scoring methods
otScoreMethod scoreMethod;
bool calcCfd = false;
bool calcMit = false;
double maximum_sum;

uint64_t sequenceToSignature(const char* seq, uint64_t seqLen)
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

void accumulate_scores(double& totScoreMit, double& totScoreCfd, bool& checkNextOfftargets, uint64_t mismatches, uint64_t dist, uint64_t occurrences, uint64_t searchSignature, uint64_t offTargetSignature) {
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

        /** Stop calculating global score early if possible */
        if (scoreMethod == otScoreMethod::mitAndCfd) {
            if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                checkNextOfftargets = false;
                return;
            }
        }
        if (scoreMethod == otScoreMethod::mitOrCfd) {
            if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                checkNextOfftargets = false;
                return;
            }
        }
        if (scoreMethod == otScoreMethod::avgMitCfd) {
            if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                checkNextOfftargets = false;
                return;
            }
        }
        if (scoreMethod == otScoreMethod::mit) {
            if (totScoreMit > maximum_sum) {
                checkNextOfftargets = false;
                return;
            }
        }
        if (scoreMethod == otScoreMethod::cfd) {
            if (totScoreCfd > maximum_sum) {
                checkNextOfftargets = false;
                return;
            }
        }
    }
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

int main(int argc, char** argv)
{
    auto startLoading = std::chrono::high_resolution_clock::now();

    if (argc < 4) {
        fprintf(stderr, "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method]\n", argv[0]);
        exit(1);
    }

    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);

    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);
    maximum_sum = (10000.0 - threshold * 100) / threshold;

    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    string argScoreMethod = argv[5];
    if (!argScoreMethod.compare("and")) {
        scoreMethod = otScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("or")) {
        scoreMethod = otScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("avg")) {
        scoreMethod = otScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("mit")) {
        scoreMethod = otScoreMethod::mit;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("cfd")) {
        scoreMethod = otScoreMethod::cfd;
        calcCfd = true;
    }
    else
    {
        fprintf(stderr, "Invalid scoring method. Acceptable options are: 'and', 'or', 'avg', 'mit', 'cfd'");
        exit(1);
    }

    /** Begin reading the binary encoded ISSL, structured as:
     *  - The header (3 items)
     *  - Slice masks
     *  - Size of slice 1 lists
     *  - Off-target signatures of slice 1 lists
     *  - Id/occurrences of slice 1 lists
     *  ...
     *  - Size of slice N lists (N being the number of slices)
     *  - Off-target signatures of slice N lists
     *  - Id/occurrences of slice N lists
     */
    mapped_file_source isslFp;
    isslFp.open(argv[1]);

    if (!isslFp.is_open())
    {
        throw std::runtime_error("Error reading index: could not open file\n");
    }
    const void* inFileFp = static_cast<const void*>(isslFp.data());

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the number of slices
     */
    const size_t* headerPtr = static_cast<const size_t*>(inFileFp);
    size_t offtargetsCount = *headerPtr++;
    size_t seqLength = *headerPtr++;
    size_t sliceCount = *headerPtr++;

    /** Read the slice masks and generate 2 bit masks */
    const uint64_t* sliceMasksPtr = static_cast<const uint64_t*>(headerPtr);
    vector<vector<uint64_t>> sliceMasks;
    for (size_t i = 0; i < sliceCount; i++)
    {
        vector<uint64_t> mask;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (*sliceMasksPtr & (1ULL << j))
            {
                mask.push_back(j);
            }
        }
        sliceMasks.push_back(mask);
        sliceMasksPtr++;
    }

    /** The contents of the slices. Stored by slice
    * Contains:
    *   - Size of each list within the slice stored contiguously
    *   - The off-target signatures of all the lists stored contiguously
    *   - The id/occurrences of all the lists stored contiguously
    */
    vector<const size_t*> allSlicelistSizes(sliceCount);
    vector<const uint64_t*> allSliceSignatures(sliceCount);
    vector<const uint64_t*> allSliceIdOccurrences(sliceCount);
    const size_t* listSizePtr = static_cast<const size_t*>(sliceMasksPtr);
    for (size_t i = 0; i < sliceCount; i++)
    {
        allSlicelistSizes[i] = listSizePtr;
        const uint64_t* signaturePtr = static_cast<const uint64_t*>(static_cast<const void*>(listSizePtr + (1ULL << (sliceMasks[i].size() * 2))));
        allSliceSignatures[i] = signaturePtr;
        const uint64_t* idOccPtr = signaturePtr + offtargetsCount;
        allSliceIdOccurrences[i] = idOccPtr;
        listSizePtr = static_cast<const size_t*>(static_cast<const void*>(idOccPtr + offtargetsCount));
    }

    /** Prevent assessing an off-target site for multiple slices
     *
     *      Create enough 1-bit "seen" flags for the off-targets
     *      We only want to score a candidate guide against an off-target once.
     *      The least-significant bit represents the first off-target
     *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
     *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
     */
    uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

    /** Start constructing index in memory
     *
     *      To begin, reverse the contiguous storage of the slices,
     *         into the following:
     *
     *         + Slice 0 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |----  ...
     *         |
     *         + Slice 1 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |---- ...
     *         | ...
     */

    vector<vector<const uint64_t*>> sliceListsSig(sliceCount);
    vector<vector<const uint64_t*>> sliceListsIdOcc(sliceCount);
    // Assign sliceLists size based on each slice length
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceListsSig[i] = vector<const uint64_t*>(1ULL << (sliceMasks[i].size() * 2));
        sliceListsIdOcc[i] = vector<const uint64_t*>(1ULL << (sliceMasks[i].size() * 2));
    }

    for (size_t i = 0; i < sliceCount; i++) {
        const uint64_t* sigList = allSliceSignatures[i];
        const uint64_t* idOccList = allSliceIdOccurrences[i];
        size_t sliceLimit = 1ULL << (sliceMasks[i].size() * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            sliceListsSig[i][j] = sigList;
            sliceListsIdOcc[i][j] = idOccList;
            sigList += allSlicelistSizes[i][j];
            idOccList += allSlicelistSizes[i][j];
        }
    }

    auto endLoading = std::chrono::high_resolution_clock::now();
    auto startProcessing = std::chrono::high_resolution_clock::now();

    //TODO: rewrite
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    size_t seqLineLength = seqLength + 1;
    std::filesystem::path queryFile(argv[2]);
    size_t fileSize = std::filesystem::file_size(queryFile);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    FILE* fp = fopen(argv[2], "rb");
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

    /** Begin scoring */
    #pragma omp parallel
    {
        vector<uint64_t> offtargetToggles(numOfftargetToggles);
        uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;
        /** For each candidate guide */
        // TODO: update to openMP > v2 (Use clang compiler)
        #pragma omp for
        for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

            auto searchSignature = querySignatures[searchIdx];
            __m512i searchSignatureVec = _mm512_set1_epi64((int64_t)searchSignature);

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            bool checkNextOfftargets = true;

            size_t sliceLimitOffset = 0;
            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                vector<uint64_t>& sliceMask = sliceMasks[i];

                uint64_t searchSlice = 0ULL;
                for (int j = 0; j < sliceMask.size(); j++)
                {
                    searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                }

                size_t idx = sliceLimitOffset + searchSlice;

                size_t signaturesInSlice = allSlicelistSizes[i][searchSlice];
                const uint64_t* sigOffset = sliceListsSig[i][searchSlice];
                const uint64_t* idOccOffset = sliceListsIdOcc[i][searchSlice];

                /** For each off-target signature in slice */
                size_t j = 0;
                for (; j + 8 <= signaturesInSlice; j += 8) {
                    _mm_prefetch((const char*)&sigOffset[j + 8], _MM_HINT_T0);
                    _mm_prefetch((const char*)&idOccOffset[j + 8], _MM_HINT_T0);

                    __m512i offTargetsVec = _mm512_loadu_si512((__m512i*)&sigOffset[j]);
                    __m512i idOccVec = _mm512_loadu_si512((__m512i*)&idOccOffset[j]);
                    __m512i signatureIdVec = _mm512_and_si512(idOccVec, _mm512_set1_epi64(0xFFFFFFFFULL));
                    __m512i occurencesVec = _mm512_srli_epi64(idOccVec, 32);

                    __m512i xoredSignaturesVec = _mm512_xor_si512(searchSignatureVec, offTargetsVec);
                    __m512i evenBitsVec = _mm512_and_si512(xoredSignaturesVec, _mm512_set1_epi64(0xAAAAAAAAAAAAAAAAULL));
                    __m512i oddBitsVec = _mm512_and_si512(xoredSignaturesVec, _mm512_set1_epi64(0x5555555555555555ULL));
                    __m512i mismatchesVec = _mm512_or_si512(_mm512_srli_epi64(evenBitsVec, 1), oddBitsVec);
                    __m512i distVec = _mm512_popcnt_epi64(mismatchesVec);

                    uint64_t mismatchesArr[8];
                    uint64_t distArr[8];

                    _mm512_storeu_si512((__m512i *)mismatchesArr, mismatchesVec);
                    _mm512_storeu_si512((__m512i *)distArr, distVec);

                    for (int lane = 0; lane < 8; lane++) {
                        uint64_t dist = distArr[lane];
                        uint64_t mismatches = mismatchesArr[lane];

                        uint64_t offTargetSignature = sigOffset[j + lane];
                        uint64_t idOccVal = idOccOffset[j + lane];
                        uint64_t signatureId = idOccVal & 0xFFFFFFFFULL;
                        uint32_t occurrences = (uint32_t)(idOccVal >> 32);

                        if (seenOfftargetAlready(offtargetTogglesTail, signatureId)) continue;

                        accumulate_scores(totScoreMit, totScoreCfd, checkNextOfftargets, mismatches, dist, occurrences, searchSignature, offTargetSignature);
                        if (!checkNextOfftargets) break;
                    }
                    if (!checkNextOfftargets) break;
                }
                // Clean-up loop
                for (; j < signaturesInSlice && checkNextOfftargets; j++) {
                    uint64_t offTargetSignature = sigOffset[j];
                    uint64_t idOccVal = idOccOffset[j];
                    uint64_t signatureId = idOccVal & 0xFFFFFFFFULL;
                    uint32_t occurrences = (uint32_t)(idOccVal >> 32);

                    if (seenOfftargetAlready(offtargetTogglesTail, signatureId)) continue;

                    uint64_t xoredSignatures = searchSignature ^ offTargetSignature;
                    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                    uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                    uint64_t mismatches = (evenBits >> 1) | oddBits;
                    uint64_t dist = popcount64(mismatches);

                    accumulate_scores(totScoreMit, totScoreCfd, checkNextOfftargets, mismatches, dist, occurrences, searchSignature, offTargetSignature);
                    if (!checkNextOfftargets) break;
                }
                if (!checkNextOfftargets) break;
                sliceLimitOffset += 1ULL << (sliceMasks[i].size() * 2);
            }
            querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
            querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

            memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
        }
    }

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
