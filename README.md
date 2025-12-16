# SAMprojet

A small Python utility to parse and analyze SAM alignment files (SAMtools-format).

This repository contains `src/sam_reader.py`, a script that reads SAM files, classifies reads (fully mapped, partially mapped, unmapped), computes per-chromosome counts and MAPQ distributions, and supports filtering by mapping status and MAPQ.

## Features

- Read and parse SAM files with headers and alignment records.
- Classify reads into: fully mapped, partially mapped, unmapped.
- Paired-end template analysis (one full + one unmapped, one full + one partial, fully mapped pairs).
- Reads per chromosome counts.
- MAPQ distribution analysis.
- CLI options to filter records before analysis:
  - `--filter` (all | fully | unmapped | partial)
  - `--min-mapq` and `--max-mapq` to restrict MAPQ range
  - `--counts-only` to print only counts after filtering

## Requirements

- Python 3.8+ (or 3.x)
- No external packages required (uses only stdlib: `argparse`, `re`, etc.)

## Files of interest

- `src/sam_reader.py` — main analysis script.
- `data/` — example SAM files (if present).
- `results/` — output or summary files produced by the script.

## Usage

Open a terminal (PowerShell on Windows) and run:

```powershell
# Basic summary printed to stdout
python src/sam_reader.py -i data/test.sam

# Save the full summary to a file
python src/sam_reader.py -i data/test.sam -o results/summary.txt

# Print only counts (total, unmapped, fully mapped, partially mapped)
python src/sam_reader.py -i data/test.sam --counts-only

# Count only fully mapped reads
python src/sam_reader.py -i data/test.sam --filter fully --counts-only

# Analyze reads with MAPQ in [30,254]
python src/sam_reader.py -i data/test.sam --min-mapq 30 --max-mapq 254

# Combine filters and save
python src/sam_reader.py -i data/test.sam --filter partial --min-mapq 20 -o results/partial_mapq20_summary.txt
```

Notes:
- `--filter` defaults to `all` (no status restriction).
- MAPQ filtering excludes records with non-numeric MAPQ when numeric bounds are provided.

## Examples

After running with `--counts-only` you will get an output like:

```
Total considered reads: 10000
Unmapped reads: 1200
Fully mapped reads: 7000
Partially mapped reads: 1800
```

## Development / Implementation notes

- The script parses the 11 required SAM fields and stores records as dictionaries.
- Numeric conversions for `POS`, `MAPQ`, `PNEXT`, and `TLEN` are done robustly (fall back to original string on failure).
- A small `readCigar` helper parses CIGAR strings into operation counts.

If you want additional output modes (for example dumping filtered SAM lines to a file, or producing CSV output), I can add a `--dump-filtered <file>` option.

## License

See the license block at the top of `src/sam_reader.py` (GPL v3 or later in the file header).
