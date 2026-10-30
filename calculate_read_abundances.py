"""
Author: David S. Dooley
Date: 10/29/2026
Input(s): Sequence alignment map (SAM) file(s)
Output: Comma-separated values (CSV) file containing absolute and fractional mapped bases per reference.
"""

import time
import pysam
import pandas as pd
import argparse
from pathlib import Path

valid_operations = [0, 2, 7, 8] # These are valid CIGAR operations for mapped bases: M, D, =, X

def calculate_read_abundances(sam_files):
    master_df = pd.DataFrame()  # Initialize an empty DataFrame to hold the merged results
    for sam_file in sam_files:
        sample_name = sam_file.stem
        abs_dict = {}
        f = pysam.AlignmentFile(sam_file, "r")
        for read in f.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue # Ignore unmapped reads
            ref_name = read.reference_name
            ref_len = sum([len for op, len in read.cigartuples if op in valid_operations]) # Sum up lengths of valid CIGAR operations for read
            if ref_name not in abs_dict:
                abs_dict[ref_name] = 0
            abs_dict[ref_name] += ref_len

        # Convert to DF 
        abundance_df = pd.DataFrame.from_dict(abs_dict, orient='index', columns=[sample_name+'_mapped_bases'])
        # Calculate relative abundances
        abundance_df[sample_name + "_fraction_mapped_bases"] = abundance_df[sample_name+"_mapped_bases"] / abundance_df[sample_name+"_mapped_bases"].sum()
        master_df = pd.concat([master_df, abundance_df], axis=1)  # Merge with the master DataFrame
    return master_df

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        epilog="Example: python3 calculate_read_abundances.py <SAM file 1> <SAM file 2> ... <SAM file n> -o output.csv",
    )

    parser.add_argument(
        "inputs",
        type=Path,
        nargs="+",
        help="One or more input SAM files to process.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default="merged_output.csv",
        help="Path to the output file (Default: %(default)s).",
    )

    return parser.parse_args()

def main():
    # Start timer
    start_time = time.perf_counter()
    # Parse the CLI inputs
    args = parse_arguments()
    # Calculate read abundances for each SAM file
    abundance_df = calculate_read_abundances(args.inputs)
    # Save the results to a CSV file
    abundance_df.to_csv(args.output)
    stop_time = time.perf_counter()
    print(f"Reads tallied in {round(stop_time - start_time, 2)} s for {len(args.inputs)} SAM files.\nData saved to {args.output}")

if __name__ == "__main__":
    main()
