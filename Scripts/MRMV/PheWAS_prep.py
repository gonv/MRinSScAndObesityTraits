#!/usr/bin/env python3

"""
Prepare NEALE-style PheWAS summary statistics for clumping in MVMR.
This script is intended to process input GWAS files with colon-separated SNP IDs
and outputs a standardized version for downstream use in the pipeline.

Usage:
    python PheWAS_prep.py input_file.txt

Output:
    {input_file}_MRMV : formatted version for MR
"""

import sys
import pandas as pd

def main():
    if len(sys.argv) < 2:
        print("Usage: python PheWAS_prep.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = f"{input_file}_MRMV"

    # Read the input data (expect tab-delimited)
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Example processing: remove duplicates, fix SNP format, or subset columns
    # You should adapt this to the actual structure of your input files

    # Save reformatted file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Formatted file saved to {output_file}")

if __name__ == "__main__":
    main()