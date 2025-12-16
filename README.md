# SAMprojet

Small Python utility to parse and analyze SAM alignment files (SAM format).

This repository contains `src/sam_reader.py`, a command-line tool that:

- Parses SAM files (headers + alignment records).
- Classifies reads into: fully mapped, partially mapped, unmapped.
- Performs paired-end template analysis (one full + one unmapped, one full + one partial, fully mapped pairs).
- Computes reads-per-chromosome counts and MAPQ distribution.
- Supports filtering by mapping status and MAPQ, and writes a human-readable summary.

Recent changes (important)
- The CIGAR global summary is now computed in-memory from the parsed records (no intermediate files required).
- `globalPercentCigar(records)` returns a formatted string and is appended to the main output (file or stdout).
- Filtering helpers (`mapq_ok` and `status_ok`) are module-level functions to improve testability.

## Features

- Read and parse SAM files with headers and alignment records.
- Classify reads into: fully mapped, partially mapped, unmapped.
- Paired-end template analysis.
- Reads per chromosome counts.
- MAPQ distribution analysis.
- CLI options to filter records before analysis and control output.

## Requirements

- Python 3.8+ (Python 3.x recommended)
- No external packages required (uses stdlib: `argparse`, `re`, etc.)

## Files of interest

- `src/sam_reader.py` — main analysis script.
- `data/` — example SAM files (if present).
- `results/` — suggested place for summary outputs.

## CLI options (summary)

- `-i, --input <file.sam>`: input SAM file (required)
- `-o, --output <file>`: write the summary into `<file>` (if provided). The CIGAR summary is appended to this file.
- `--filter {all,fully,unmapped,partial}`: restrict analysis to a subset of reads (default: `all`).
- `--min-mapq <int>`: minimum MAPQ (inclusive) to include.
- `--max-mapq <int>`: maximum MAPQ (inclusive) to include.
- `--counts-only`: only print the high-level counts (total / unmapped / fully / partially) after filters.

Behavior notes:
- When `-o/--output` is provided, the script writes the main summary and then appends the computed CIGAR operation percentages to the same file.
- When no `-o` is provided, the full summary (including the CIGAR block) is printed to stdout.

## Usage examples

```powershell
# Basic summary printed to stdout (includes CIGAR block)
python .\src\sam_reader.py -i .\data\test.sam

# Save the full summary (summary + CIGAR block) to a file
python .\src\sam_reader.py -i .\data\test.sam -o .\results\summary.txt

# Print only counts (after filters)
python .\src\sam_reader.py -i .\data\test.sam --counts-only

# Count only fully mapped reads and save
python .\src\sam_reader.py -i .\data\test.sam --filter fully --counts-only -o .\results\fully_counts.txt

# Restrict by MAPQ (example)
python .\src\sam_reader.py -i .\data\test.sam --min-mapq 30 --max-mapq 254 -o .\results\mapq30_254_summary.txt
```

## Example output (full summary)

The appended CIGAR block looks like:

```
Global CIGAR operation percentages (across considered records):
  M: 85.12%
  I: 3.45%
  D: 1.02%
  S: 5.21%
  H: 0.00%
  N: 0.00%
  P: 0.00%
  X: 5.20%
  =: 0.00%
```

## Implementation notes

- The script parses the 11 required SAM fields and stores records as dictionaries.
- Numeric conversions for `POS`, `MAPQ`, `PNEXT`, and `TLEN` are robust (fall back to the original string if conversion fails).
- `readCigar` parses CIGAR strings into per-operation counts; `globalPercentCigar(records)` aggregates those counts and returns percentages.
- Filtering helpers `mapq_ok(rec, min_mapq, max_mapq)` and `status_ok(rec, filter_mode)` are available at module-level to enable unit testing.

## Extensions you may want

- `--dump-filtered <file>`: export the filtered SAM lines to a file for downstream processing.
- CSV/TSV export of summary tables.
- Add unit tests for `readCigar`, `globalPercentCigar`, and filtering helpers.

## License

See the license block at the top of `src/sam_reader.py` (GPL v3 or later in the file header).
