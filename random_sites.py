#!/usr/bin/env python3
"""

"""
import argparse
import dnaio
import random


def random_sites(s):
    s = s.upper()
    nickase_sites = sum(s.count(pattern) for pattern in ["GAGAC", "GTCTC"])

    # NOTE This simulates overlapping sites which we would not get with Nickase
    for start in random.sample(range(len(s)), k=nickase_sites):
        yield start, start + 5


def main():
    random.seed(0)
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("fasta")
    args = parser.parse_args()

    with dnaio.open(args.fasta) as f:
        for record in f:
            for start, end in random_sites(record.sequence):
                print(record.id, start, end, sep="\t")


if __name__ == "__main__":
    main()

