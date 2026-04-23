
#!/usr/bin/env python3
"""
Parse IQTree .iqtree and .partlh/.sitelh files to produce a TSV with:
PartitionName or SiteIndex, Tree1, Tree2, ..., Delta(Tree1-Tree2) (if two trees).
"""

import csv
import sys
import os

def parse_iqtree_partitions(iqtree_file):
    """Extract partition names from IQTree log file."""
    partition_names = []
    with open(iqtree_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            # Partition rows have numeric ID in first column and name in last column
            if len(parts) >= 9 and parts[0].isdigit():
                partition_names.append(parts[-1])
    return partition_names

def parse_partlh(partlh_file):
    """Parse .partlh file: header + tree rows."""
    with open(partlh_file, "r") as f:
        lines = [line.strip().split() for line in f if line.strip()]
    # First line: num_trees and num_partitions
    try:
        num_trees, num_parts = map(int, lines[0])
    except ValueError:
        raise ValueError(f"Expected two integers in first line of {partlh_file}, got: {lines[0]}")
    # Next num_trees lines: tree name + likelihoods
    tree_ids = []
    tree_values = []
    for i in range(1, num_trees + 1):
        tree_ids.append(lines[i][0])
        tree_values.append(list(map(float, lines[i][1:])))
    return tree_ids, tree_values, num_parts

def parse_sitelh(sitelh_file):
    """Parse .sitelh file: first row = tree names, rest = site likelihoods."""
    with open(sitelh_file, "r") as f:
        lines = [line.strip().split() for line in f if line.strip()]
    tree_ids = lines[0]
    data_rows = [list(map(float, row)) for row in lines[1:]]
    return tree_ids, data_rows

def create_tsv(iqtree_file, likelihood_file, output_file):
    # Detect file type by name
    if ".partlh" in likelihood_file:
        file_type = "partition"
        partition_names = parse_iqtree_partitions(iqtree_file)
        tree_ids, tree_values, num_parts = parse_partlh(likelihood_file)

        header = ["PartitionName"] + tree_ids
        if len(tree_ids) == 2:
            header.append("Delta(Tree1-Tree2)")

        with open(output_file, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(header)

            for idx in range(num_parts):
                name = partition_names[idx] if idx < len(partition_names) else f"Partition_{idx+1}"
                row = [name] + [tree_values[t][idx] for t in range(len(tree_ids))]
                if len(tree_ids) == 2:
                    row.append(row[1] - row[2])  # Tree1 - Tree2
                writer.writerow(row)

    elif ".sitelh" in likelihood_file:
        file_type = "site"
        tree_ids, data_rows = parse_sitelh(likelihood_file)

        header = ["SiteIndex"] + tree_ids
        if len(tree_ids) == 2:
            header.append("Delta(Tree1-Tree2)")

        with open(output_file, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(header)

            for idx, values in enumerate(data_rows):
                row = [idx + 1] + values
                if len(values) == 2:
                    row.append(values[0] - values[1])
                writer.writerow(row)

    else:
        raise ValueError("File name must include .partlh or .sitelh to determine type.")

    print(f"TSV file '{output_file}' created successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python lnleval.py <iqtree_file> <likelihood_file> <output_tsv>")
        sys.exit(1)

    iqtree_file = sys.argv[1]
    likelihood_file = sys.argv[2]
    output_file = sys.argv[3]

    create_tsv(iqtree_file, likelihood_file, output_file)
