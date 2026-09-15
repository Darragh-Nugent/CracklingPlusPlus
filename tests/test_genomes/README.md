# Automated Genome Benchmarking

This script automatically extracts off-targets, creates an ISSL index, and benchmarks off-target scoring across a set of genomes.

## Directory Structure

```text
CracklingPlusPlus/
    ├── build/
    │   ├── ExtractOfftargets/
    │   │   └── ExtractOfftargets
    │   ├── ISSLCreateIndex/
    │   │   └── ISSLCreateIndex
    │   └── ISSLScoreOfftargets/
    │       └── ISSLScoreOfftargets
    │
    └── tests/
        ├── test.py
        ├── test_genomes/
        │   ├── genome_1/
        │   │   ├── genome.fna
        │   │   ├── guides.txt
        │   │   └── guides.fa
        │   │
        │   └── genome_2/
        │       ├── genome.fna
        │       └── guides.txt
        │
        └── results.txt
```
