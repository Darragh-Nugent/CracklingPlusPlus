#!/usr/bin/env python3

import sys
import subprocess
import time
from pathlib import Path

# Configurable parameters
slice_config = "slice-8-20"
sequence_length = "20"
max_distance = "4"
threshold = "0"
score_method = "mit"

SCRIPT_DIR = Path(__file__).parent

test_genomes_dir = SCRIPT_DIR / "test_genomes"
RESULTS_FILE = SCRIPT_DIR / "results.txt"

off_target_scoring_binary = (SCRIPT_DIR / "../build/ISSLScoreOfftargets/ISSLScoreOfftargets").resolve()
extract_off_targets_binary = (SCRIPT_DIR / "../build/ExtractOfftargets/ExtractOfftargets").resolve()
issl_create_index_binary = (SCRIPT_DIR / "../build/ISSLCreateIndex/ISSLCreateIndex").resolve()
slice_config_path = (SCRIPT_DIR / f"../sample/{slice_config}.txt")

def score_off_targets(genome_folder: Path, guides: Path, number_of_guides: int | None) -> float:
    issl_index = "issl.index"
    print(f"Scoring {genome_folder.name} against {number_of_guides if number_of_guides is not None else "UNKNOWN"} guides...")
    start = time.perf_counter()
    subprocess.run(
        [off_target_scoring_binary, issl_index, guides, max_distance, threshold, score_method],
        cwd=genome_folder,
        stdout=subprocess.DEVNULL,
        check=True)
    runtime = time.perf_counter() - start
    print("Scoring complete.")
    return runtime

def extract_off_targets(genome_folder: Path, genome_file: Path) -> None:
    output_file = "off_targets.txt"
    print(f"Extracting off-targets for {genome_folder.name}...")
    subprocess.run(
        [extract_off_targets_binary, output_file, genome_file],
        cwd=genome_folder,
        stdout=subprocess.DEVNULL,
        check=True)
    print("Off-target extraction complete.")

def create_issl_index(genome_folder: Path) -> None:
    off_targets = "off_targets.txt"
    output_file = "index.issl"
    print(f"Creating index for {genome_folder.name}")
    subprocess.run(
        [issl_create_index_binary, off_targets, slice_config_path, sequence_length, output_file],
        cwd=genome_folder,
        stdout=subprocess.DEVNULL,
        check=True)
    print("Index creation complete.")

def count_lines(file: Path) -> int:
    result = subprocess.run(
        ["wc", "-l", str(file)],
        capture_output=True,
        text=True)

    return int(result.stdout.split()[0])

def write_results(genome_name: str, number_of_guides: int, runtime: float):
    with RESULTS_FILE.open("a") as f:
        f.write(
            f"{genome_name:<20} | {number_of_guides:>10,} | {runtime:>7.3f}s\n"
        )

if not test_genomes_dir.exists():
    sys.exit("Failed to locate test_genomes directory.") 

if RESULTS_FILE.exists():
    sys.exit("Results file already exists... exiting.") 

with RESULTS_FILE.open("w") as f:
    f.write(
        f"{'Genome':<20} | {'Guides':>10} | {'Runtime':>8}\n"
        f"{'-' * 20}-+-{'-' * 10}-+-{'-' * 8}\n"
    )

for genome_folder in test_genomes_dir.iterdir():
    if not genome_folder.is_dir():
        continue

    genome_name = genome_folder.name
    genome_file = list(genome_folder.glob("*.fna"))[0] # assume only one .fna file per genome dir

    if all(file.name != "off_targets.txt" for file in genome_folder.glob("*.txt")):
        extract_off_targets(genome_folder, genome_file)
    else:
        print("off_targets.txt file detected... skipping off-target extraction")

    if all(file.name != "index.issl" for file in genome_folder.glob("*.issl")):
        create_issl_index(genome_folder)
    else:
        print("index.issl file detected... skipping ISSL index creation")

    for guides in list(genome_folder.glob("*.txt")) + list(genome_folder.glob("*.fa")):
        if guides.name == "off_targets.txt":
            continue
        number_of_guides = count_lines(guides)
        runtime = score_off_targets(genome_folder, guides, number_of_guides)  
        write_results(genome_name, number_of_guides, runtime)      

sys.exit(0) 
