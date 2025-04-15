import numpy as np
import pandas as pd
import os, sys, argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="This script will edit multiple gff files with fasta information and add sequence-region headers.")
    parser.add_argument("-csv", "--csv", type=str, help="CSV file with columns: gff filename, fasta filename, output filename")
    parser.add_argument("-gff_dir", "--gff_dir", type=str, help="Directory containing GFF files")
    parser.add_argument("-fasta_dir", "--fasta_dir", type=str, help="Directory containing FASTA files")
    parser.add_argument("-output_dir", "--output_dir", type=str, help="Directory to save modified GFF files")
    return parser.parse_args()

def read_fasta(fasta_path):
    fasta_dict = {}
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                fasta_dict[header] = ""
            else:
                fasta_dict[header] += line.strip()
    return fasta_dict

def extract_sequence_region(gff_lines):
    sequence_regions = []
    pattern = re.compile(r'# Sequence Data: seqnum=\d+;seqlen=(\d+);seqhdr="([^"]+)"')
    
    for line in gff_lines:
        match = pattern.search(line)
        if match:
            seq_length = match.group(1)
            seq_header = match.group(2)
            sequence_regions.append(f"##sequence-region {seq_header} 1 {seq_length}\n")
            sequence_regions.append("##description\n")  # Add blank description line
    
    return sequence_regions

def read_gff(gff_path):
    with open(gff_path, "r") as f:
        gff_lines = f.readlines()
    return gff_lines

def write_gff(gff_lines, fasta_dict, output_path):
    sequence_regions = extract_sequence_region(gff_lines)
    
    with open(output_path, "w") as f:
        for line in gff_lines:
            if line.startswith("##gff-version"):
                f.write(line)
                f.writelines(sequence_regions)  # Insert sequence-region headers after version line
            elif not line.startswith("# Sequence Data") and not line.startswith("# Model Data"):
                f.write(line)
        
        f.write("##FASTA\n")
        for header in fasta_dict:
            f.write(header + "\n")
            f.write(fasta_dict[header] + "\n")

def process_files(csv_path, gff_dir, fasta_dir, output_dir):
    df = pd.read_csv(csv_path)
    for _, row in df.iterrows():
        gff_path = os.path.join(gff_dir, row['gff'])
        fasta_path = os.path.join(fasta_dir, row['fasta'])
        output_path = os.path.join(output_dir, row['output'])
        
        if os.path.exists(gff_path) and os.path.exists(fasta_path):
            fasta_dict = read_fasta(fasta_path)
            gff_lines = read_gff(gff_path)
            write_gff(gff_lines, fasta_dict, output_path)
        else:
            print(f"Skipping {row['gff']} or {row['fasta']} - File not found")

def main():
    args = parse_args()
    process_files(args.csv, args.gff_dir, args.fasta_dir, args.output_dir)
    
if __name__ == "__main__":
    main()
