#!/usr/bin/env python3
"""
Find where the Nt.BsmAI nickase (restriction enzyme) would cut a reference
and output the sites in BED format.

Run once without and once with --rc to get sites both on the forward and
reverse strand.
"""
import argparse
import dnaio


def restriction_sites(s, rc):
    s = s.upper()
    pattern = "GAGAC" if rc else "GTCTC"
    m = len(pattern)
    start = 0
    while (pos := s.find(pattern, start)) != -1:
        yield (pos, pos + m)
        start = pos + m


def main():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--rc", action="store_true", default=False)
    parser.add_argument("fasta")
    args = parser.parse_args()

    rc = args.rc
    with dnaio.open(args.fasta) as f:
        for record in f:
            for start, end in restriction_sites(record.sequence, rc=rc):
                print(record.id, start, end, sep="\t")


if __name__ == "__main__":
    main()

