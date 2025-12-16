#!/usr/bin/python3
#-*- coding : utf-8 -*-


__authors__ = ("Thai-Son Nguyen" )
__contact__ = ("thai-son.nguyen01@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/16/2025"
__licence__ = """This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>."""


     
    ### OPTION LIST:
        ## -h, --help           : show help and exit
        ## -i, --input          : input SAM file (required)
        ## -o, --output         : output file to write the summary (default: stdout)
        ## --filter            : restrict analysis to a subset of reads; choices: all, fully, unmapped, partial (default: all)
        ## --min-mapq          : minimum MAPQ (inclusive) to include in analysis
        ## --max-mapq          : maximum MAPQ (inclusive) to include in analysis
        ## --counts-only       : only print counts (total/unmapped/fully/partially) after filters

    # Synopsis examples:
        ## SamReader.py -i <file.sam>
        ##     Basic analysis and summary printed to stdout.

        ## SamReader.py -i <file.sam> -o summary.txt
        ##     Save full summary to `summary.txt`.

        ## SamReader.py -i <file.sam> --counts-only
        ##     Print only counts (total, unmapped, fully mapped, partially mapped).

        ## SamReader.py -i <file.sam> --filter fully --counts-only
        ##     Count only reads that are fully mapped.

        ## SamReader.py -i <file.sam> --min-mapq 30 --max-mapq 254
        ##     Analyse only reads with MAPQ in [30,254].
  


############### IMPORT MODULES ###############

import os, sys, re, argparse

def checkFile(filePath):
    """
    Validate that the directory (if provided) exists, the file exists, and has a .sam extension;
    exit with an error message if any check fails.
    """
    import os, sys
    dirn = os.path.dirname(filePath) or '.'
    # If a directory component is present, ensure it exists
    if dirn and not os.path.exists(dirn):
        print("Error: The directory " + dirn + " does not exist.")
        sys.exit(1)
    # Ensure the file exists
    if not os.path.isfile(filePath):
        print("Error: The file " + filePath + " does not exist.")
        sys.exit(1)
    # Ensure correct extension
    if not filePath.endswith(".sam"):
        print("Error: The file " + filePath + " is not a .sam file.")
        sys.exit(1)
    return True
## 2/ Read, 

def read_sam(file_path):
    headers = []   # contains header lines
    records = []   # contains alignment records (list of lists)

    # 1. Open file using with open(...)
    with open(file_path, "r") as fd:

        # 2. Iterate line by line
        for line in fd:
            line = line.rstrip("\n")  # remove trailing newline

            # 3. Distinguish header
            if line.startswith("@"):
                headers.append(line)
                continue  # skip to next line

            # skip empty lines
            if line.strip() == "":
                continue

            # 4. Parse alignment fields (SAM has 11 mandatory fields)
            fields = line.split("\t")

            # OPTIONAL: validate number of columns (11)
            if len(fields) < 11:
                print(f"[WARNING] Invalid SAM alignment line ignored:\n{line}")
                continue

            records.append(fields)

    return headers, records



## 3/ Store,
def store_records(headers, records_raw):
    stored_records = []  # list storing alignment dictionaries

    for cols in records_raw:

        # Convert the 11 mandatory fields into a dictionary following the SAM specification
        # Use safe conversions for numeric fields (POS, MAPQ, PNEXT, TLEN)
        try:
            pos = int(cols[3])
        except Exception:
            pos = cols[3]

        try:
            mapq = int(cols[4])
        except Exception:
            mapq = cols[4]

        try:
            pnext = int(cols[7])
        except Exception:
            pnext = cols[7]

        try:
            tlen = int(cols[8])
        except Exception:
            tlen = cols[8]

        record = {
            "QNAME": cols[0],
            "FLAG": int(cols[1]),
            "RNAME": cols[2],
            "POS": pos,
            "MAPQ": mapq,
            "CIGAR": cols[5],
            "RNEXT": cols[6],
            "PNEXT": pnext,
            "TLEN": tlen,
            "SEQ": cols[9],
            "QUAL": cols[10]
        }

        stored_records.append(record)

    # Return a single data structure containing all parsed data
    stored_data = {
        "headers": headers,
        "records": stored_records
    }

    return stored_data


#### Convert the flag into binary ####
def flagBinary(flag) :
    """Return the FLAG value as a list of binary characters (MSB first).

    Pads to 12 bits to make testing specific FLAG bits easier.
    """
    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustment to 12 (maximal flag size)
        add = 12 - len(flagB) # difference between the maximal flag size (12) and the length of the binary flag
        for t in range(add):
            flagB.insert(0,'0') # pad with leading zeros until maximal flag size
    return flagB

def readCigar(cigar):
    """Parse a CIGAR string and return a dict of operation -> total count.

    Example: '10M1I5M' -> {'M': 15, 'I': 1}
    """
    ext = re.findall('\w', cigar)  # tokenize CIGAR string into letters and digits
    key = []
    value = []
    val = ""

    for i in range(0, len(ext)):  # for each character (digits or letters)
        if ext[i] in ('M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', '='):
            key.append(ext[i])
            value.append(val)
            val = ""
        else:
            val = val + ext[i]  # accumulate digit characters

    dico = {}
    n = 0
    for k in key:  # build dictionary summing counts for each operation
        if k not in dico:
            dico[k] = int(value[n])
        else:
            dico[k] += int(value[n])
        n += 1
    return dico

### Analyse the CIGAR = regular expression that summarise each read alignment ###

def globalPercentCigar(records):
    """
    Compute a global CIGAR operation percentage summary from `records`.

    - Skips unmapped records and entries with missing CIGAR ('*').
    - Returns a formatted string listing percentage for each CIGAR op.
    """
    ops = ['M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', '=']
    totals = {op: 0.0 for op in ops}
    total_counts = 0.0

    for rec in records:
        cigar = rec.get('CIGAR')
        # skip unmapped or missing CIGAR
        if not cigar or cigar == '*' or (rec.get('FLAG', 0) & 4):
            continue
        try:
            d = readCigar(cigar)
        except Exception:
            continue
        for op, val in d.items():
            try:
                cnt = float(val)
            except Exception:
                continue
            totals[op] = totals.get(op, 0.0) + cnt
            total_counts += cnt

    if total_counts == 0.0:
        return "No CIGAR information available for the selected records."

    lines = ["Global CIGAR operation percentages (across considered records):"]
    for op in ops:
        pct = (totals.get(op, 0.0) / total_counts) * 100.0
        lines.append(f"  {op}: {round(pct, 2)}%")

    return "\n".join(lines)


def analyse_fully_mapped(records):
    """
    Identify reads that are fully mapped.
    A fully mapped read = FLAG not unmapped (bit 0x4 == 0)
                       AND CIGAR contains only M (e.g., 100M).
    """
    fully_mapped = []

    for rec in records:
        flag = rec["FLAG"]
        cigar = rec["CIGAR"]

        # 1. If FLAG 0x4 is set -> read is unmapped -> skip
    
        
        #The flag &4 is a bitwise AND operation used to check if the 0x4 bit of the SAM FLAG is enabled. 
        #This bit indicates that the read is not mapped to the reference. 
        #If the operation returns a non-zero value, the read is considered unmapped and will be excluded from parsing. 
        if flag & 4: 
            continue

        # 2. Check that CIGAR contains only 'M' operations
        #    A valid fully-mapped CIGAR must match the pattern "[0-9]+M"
        import re
        if re.fullmatch(r"\d+M", cigar):
            fully_mapped.append(rec)

    return fully_mapped

def group_by_qname(records):
    """Group a list of records by their `QNAME`.

    Returns a dict mapping QNAME -> list_of_records.
    """
    groups = {}
    for rec in records:
        q = rec["QNAME"]
        if q not in groups:
            groups[q] = []
        groups[q].append(rec)
    return groups
def is_fully_mapped(rec):
    """Return True if record is mapped and its CIGAR contains only a single M operation.

    A fully-mapped read is defined here as having FLAG bit 0x4 unset and a CIGAR
    matching e.g. '100M'.
    """
    import re
    return (rec["FLAG"] & 4 == 0) and re.fullmatch(r"\d+M", rec["CIGAR"])


def is_partially_mapped(rec):
    """Return True for reads that are mapped but not fully mapped.

    Excludes unmapped reads (FLAG 0x4). For mapped reads, returns True when
    the CIGAR string is not a single contiguous match (not like '100M').
    """
    # Mapped but not fully mapped
    if rec["FLAG"] & 4:
        return False
    import re
    return not re.fullmatch(r"\d+M", rec["CIGAR"])

def analyse_paired_end(records):
    """
    Analyse paired-end templates:
    - one fully mapped + one unmapped
    - one fully mapped + one partially mapped
    - properly paired (both fully mapped)
    """

    pairs = group_by_qname(records)

    one_full_one_unmapped = 0
    one_full_one_partial = 0
    fully_mapped_pairs = 0

    for qname, segs in pairs.items():

        # skip templates that are not paired (FLAG 0x1)
        if len(segs) != 2:
            continue

        r1, r2 = segs[0], segs[1]

        flag1, flag2 = r1["FLAG"], r2["FLAG"]

        # Skip if not paired-end (FLAG 0x1 == 0)
        if not (flag1 & 1):
            continue

        # Evaluate mapping status
        r1_fully = is_fully_mapped(r1)
        r2_fully = is_fully_mapped(r2)

        r1_partial = is_partially_mapped(r1)
        r2_partial = is_partially_mapped(r2)

        r1_unmapped = flag1 & 4 != 0 #not logical
        r2_unmapped = flag2 & 4 != 0

        # Case 1: One fully mapped + one unmapped
        if (r1_fully and r2_unmapped) or (r2_fully and r1_unmapped):
            one_full_one_unmapped += 1
            continue

        # Case 2: One fully mapped + one partially mapped
        if (r1_fully and r2_partial) or (r2_fully and r1_partial):
            one_full_one_partial += 1
            continue

        # Case 3: Both fully mapped (proper pair)
        if r1_fully and r2_fully:
            fully_mapped_pairs += 1
            continue

    return {
        "one_full_one_unmapped": one_full_one_unmapped,
        "one_full_one_partial": one_full_one_partial,
        "fully_mapped_pairs": fully_mapped_pairs
    }

def analyse_reads_per_chromosome(records):
    """
    Count number of mapped reads per chromosome (RNAME = Column 3).
    Unmapped reads have RNAME='*' and should be excluded.
    """
    chrom_counts = {}

    for rec in records:
        rname = rec["RNAME"]

        # Skip unmapped (RNAME = "*")
        if rname == "*":
            continue

        if rname not in chrom_counts:
            chrom_counts[rname] = 0
        chrom_counts[rname] += 1

    return chrom_counts

def analyse_mapq_distribution(records):
    """
    Analyse MAPQ values across all records.
    Group into ranges required:
    - MAPQ < 30
    - 30 ≤ MAPQ < 255
    - MAPQ == 255
    """
    low = 0
    mid = 0
    high = 0

    for rec in records:
        mapq = rec["MAPQ"]

        # MAPQ is int (converted during store)
        if mapq < 30:
            low += 1
        elif mapq < 255:
            mid += 1
        else:
            high += 1

    return {
        "MAPQ_<30": low,
        "MAPQ_30_254": mid,
        "MAPQ_255": high
    }


def mapq_ok(rec, min_mapq=None, max_mapq=None):
    """
    Check whether a record's MAPQ value passes the numeric filters.

    - If MAPQ is numeric, it must satisfy min_mapq <= MAPQ <= max_mapq when those
      bounds are provided.
    - If MAPQ is non-numeric and any numeric filter is set, treat the record as
      failing the numeric filter (returns False). If no numeric filters set,
      non-numeric MAPQ is treated as acceptable.

    Returns True/False.
    """
    try:
        mq = int(rec.get('MAPQ', 0))
    except Exception:
        return (min_mapq is None) and (max_mapq is None)
    if min_mapq is not None and mq < min_mapq:
        return False
    if max_mapq is not None and mq > max_mapq:
        return False
    return True


def status_ok(rec, filter_mode='all'):
    """
    Check whether a record matches the requested mapping-status filter.

    filter_mode values:
    - 'all'     : accept any record
    - 'fully'   : only fully mapped reads (CIGAR =~ "^[0-9]+M$" and not unmapped)
    - 'unmapped': only unmapped reads (FLAG bit 0x4 set)
    - 'partial' : mapped but not fully mapped

    Returns True/False.
    """
    if filter_mode == 'all':
        return True
    if filter_mode == 'fully':
        return is_fully_mapped(rec)
    if filter_mode == 'unmapped':
        return (rec.get('FLAG', 0) & 4) != 0
    if filter_mode == 'partial':
        return is_partially_mapped(rec)
    return True

#### Summarise the results ####

def Summary(data):
    """
    Build a full analysis summary from parsed SAM records.
    Includes:
    - unmapped reads
    - fully mapped reads
    - partially mapped reads
    - paired-end analysis
    - reads per chromosome (RNAME)
    - MAPQ distribution
    """

    records = data["records"]

    # -----------------------------
    # 1) Basic Classifications
    # -----------------------------

    # Unmapped reads (FLAG 0x4)
    unmapped_reads = [rec for rec in records if rec["FLAG"] & 4]
    unmapped_count = len(unmapped_reads)

    # Fully mapped reads
    fully_mapped_reads = analyse_fully_mapped(records)
    fully_mapped_count = len(fully_mapped_reads)

    # Partially mapped reads (mapped but not fully)
    partially_mapped_reads = [
        rec for rec in records
        if not (rec["FLAG"] & 4) and not is_fully_mapped(rec)
    ]
    partially_mapped_count = len(partially_mapped_reads)

    # -----------------------------
    # 2) Paired-end analysis
    # -----------------------------
    paired_stats = analyse_paired_end(records)

    # -----------------------------
    # 3) Reads per chromosome (Q3)
    # -----------------------------
    chrom_stats = analyse_reads_per_chromosome(records)

    # -----------------------------
    # 4) MAPQ distribution (Q4)
    # -----------------------------
    mapq_stats = analyse_mapq_distribution(records)

    # -----------------------------
    # Build final summary dictionary
    # -----------------------------

    results = {
        # MAIN COUNTS
        "total_reads": len(records),

        "unmapped_count": unmapped_count,
        "fully_mapped_count": fully_mapped_count,
        "partially_mapped_count": partially_mapped_count,

        # PAIR STATS
        "one_full_one_unmapped_pairs": paired_stats["one_full_one_unmapped"],
        "one_full_one_partial_pairs": paired_stats["one_full_one_partial"],
        "fully_mapped_pairs": paired_stats["fully_mapped_pairs"],

        # CHROMOSOME DISTRIBUTION
        "reads_per_chromosome": chrom_stats,

        # MAPQ DISTRIBUTION
        "mapq_distribution": mapq_stats,

        # STORE LIST OF READS (FOR OPTIONAL SAVING)
        "unmapped_reads": unmapped_reads,
        "fully_mapped_reads": fully_mapped_reads,
        "partially_mapped_reads": partially_mapped_reads
    }

    return results

   

#### Main function ####

def main(argv=None):
    """Command-line entry point for the SAM reader analysis.

    - Parses arguments: input SAM file (required) and optional output file.
    - Runs the existing pipeline: checkFile -> read_sam -> store_records -> Summary
    - Prints the summary to stdout or writes it to `-o/--output` file.
    """
    parser = argparse.ArgumentParser(description="Simple SAM file analysis")
    parser.add_argument('-i', '--input', required=True, help='Input SAM file')
    parser.add_argument('-o', '--output', help='Output file to write summary (default: stdout)')
    parser.add_argument('--filter', choices=['all', 'fully', 'unmapped', 'partial'], default='all',
                        help='Restrict analysis to a subset of reads')
    parser.add_argument('--min-mapq', type=int, default=None, help='Minimum MAPQ to include')
    parser.add_argument('--max-mapq', type=int, default=None, help='Maximum MAPQ to include')
    parser.add_argument('--counts-only', action='store_true', help='Only print counts for unmapped/fully/partially mapped reads (after filters)')

    args = parser.parse_args(argv)

    input_file = args.input
    output_file = args.output
    filter_mode = args.filter
    min_mapq = args.min_mapq
    max_mapq = args.max_mapq
    counts_only = args.counts_only

    # Basic checks
    checkFile(input_file)

    # Read and parse
    headers, raw_records = read_sam(input_file)
    data = store_records(headers, raw_records)

    # Apply MAPQ and mapping-status filters to records before summarising
    records = data['records']

    # Use module-level helpers `mapq_ok` and `status_ok` so the filtering
    # logic is reusable and testable outside of `main`.
    filtered = [
        r for r in records
        if mapq_ok(r, min_mapq, max_mapq) and status_ok(r, filter_mode)
    ]

    # If the user wants only counts (and not a full summary), compute and print them
    if counts_only:
        total = len(filtered)
        unmapped_count = len([r for r in filtered if r.get('FLAG', 0) & 4])
        fully_count = len([r for r in filtered if is_fully_mapped(r)])
        partial_count = len([r for r in filtered if (not (r.get('FLAG', 0) & 4)) and not is_fully_mapped(r)])

        out_lines = [f"Total considered reads: {total}",
                     f"Unmapped reads: {unmapped_count}",
                     f"Fully mapped reads: {fully_count}",
                     f"Partially mapped reads: {partial_count}"]

        output_text = "\n".join(out_lines)

        if output_file:
            # Write counts first, then append the CIGAR summary into the same file
            try:
                with open(output_file, 'w') as outf:
                    outf.write(output_text)
                try:
                    cigar_text = globalPercentCigar(filtered)
                    with open(output_file, 'a') as outf:
                        outf.write("\n\n" + cigar_text)
                except Exception as e:
                    with open(output_file, 'a') as outf:
                        outf.write("\n\n" + f"Could not include global CIGAR summary: {e}")
                print(f"Counts written to {output_file}")
            except Exception as e:
                print(f"Error writing counts to {output_file}: {e}")
                print(output_text)
        else:
            # No output file: generate legacy Final_Cigar_table.txt and include it
            try:
                cigar_text = globalPercentCigar(filtered)
                output_text = output_text + "\n\n" + cigar_text
            except Exception as e:
                output_text = output_text + "\n\n" + f"Could not include global CIGAR summary: {e}"
            print(output_text)
        return

    # Otherwise compute a full summary on the filtered records
    filtered_data = {
        'headers': headers,
        'records': filtered
    }

    summary = Summary(filtered_data)

    # Build human-readable output
    lines = []
    lines.append(f"Total reads (after filters): {summary.get('total_reads', 0)}")
    lines.append(f"Unmapped reads: {summary.get('unmapped_count', 0)}")
    lines.append(f"Fully mapped reads: {summary.get('fully_mapped_count', 0)}")
    lines.append(f"Partially mapped reads: {summary.get('partially_mapped_count', 0)}")
    lines.append("")
    lines.append("Paired-end stats:")
    lines.append(f"  One full + one unmapped: {summary.get('one_full_one_unmapped_pairs', 0)}")
    lines.append(f"  One full + one partial: {summary.get('one_full_one_partial_pairs', 0)}")
    lines.append(f"  Fully mapped pairs: {summary.get('fully_mapped_pairs', 0)}")
    lines.append("")
    lines.append("MAPQ distribution:")
    md = summary.get('mapq_distribution', {})
    lines.append(f"  MAPQ < 30: {md.get('MAPQ_<30', 0)}")
    lines.append(f"  30 ≤ MAPQ < 255: {md.get('MAPQ_30_254', 0)}")
    lines.append(f"  MAPQ == 255: {md.get('MAPQ_255', 0)}")
    lines.append("")
    lines.append("Reads per chromosome:")
    for chrom, count in summary.get('reads_per_chromosome', {}).items():
        lines.append(f"  {chrom}: {count}")

    output_text = "\n".join(lines)

    if output_file:
        # Write main summary first, then append cigar summary into same file
        try:
            with open(output_file, 'w') as outf:
                outf.write(output_text)
            try:
                cigar_text = globalPercentCigar(filtered)
                with open(output_file, 'a') as outf:
                    outf.write("\n\n" + cigar_text)
            except Exception as e:
                with open(output_file, 'a') as outf:
                    outf.write("\n\n" + f"Could not include global CIGAR summary: {e}")
            print(f"Summary written to {output_file}")
        except Exception as e:
            print(f"Error writing to output file {output_file}: {e}")
            print("Falling back to printing summary to stdout:\n")
            print(output_text)
    else:
        # No output file: generate legacy Final_Cigar_table.txt and include it
        try:
            cigar_text = globalPercentCigar(filtered)
            output_text = output_text + "\n\n" + cigar_text
        except Exception as e:
            output_text = output_text + "\n\n" + f"Could not include global CIGAR summary: {e}"
        print(output_text)


if __name__ == '__main__':
    main()



