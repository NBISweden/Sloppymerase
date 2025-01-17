#!/usr/bin/env python3
"""
Find occurrences of any of the given motifs or its reverse complements
and output the sites in BED format.
"""
# /// script
# dependencies = [
#   "dnaio",
# ]
# ///
import argparse
import re

import dnaio


_REVCOMP_TRANS = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn", "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn"
)


def motif_sites(s, regex):
    for match in regex.finditer(s.upper()):
        yield match.start(), match.end()


def reverse_complement(s):
    return s.translate(_REVCOMP_TRANS)[::-1]


def main():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        "-m",
        dest="motifs",
        default=[],
        action="append",
        help="Motif to search for (also the reverse complement will be searched)",
    )
    parser.add_argument("fasta")
    args = parser.parse_args()

    motifs = args.motifs + [reverse_complement(motif) for motif in args.motifs]
    motifs = [motif.upper() for motif in motifs]
    regex = re.compile("(" + "|".join(motifs) + ")")

    with dnaio.open(args.fasta) as f:
        for record in f:
            for start, end in motif_sites(record.sequence, regex):
                print(record.id, start, end, sep="\t")


if __name__ == "__main__":
    main()
