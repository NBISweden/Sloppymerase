#!/usr/bin/env python3
"""
Find ssDNA breaks and print a BED file with detected break sites to stdout.

The input is a BAM or CRAM file with aligned reads.

The output BED file has a special "#gffTags" header that is not part of the BED
specification, but that IGV understands and makes it display some nice
annotations when hovering with the mouse over an annotated event.
"""
# /// script
# dependencies = [
#   "pyfaidx",
#   "pysam",
# ]
# ///
import sys
import argparse
from dataclasses import dataclass
from typing import Iterator

import pysam
from pyfaidx import Fasta


class CommandlineError(Exception):
    pass


@dataclass
class BreakEvent:
    start: int  # position of first mutated base (reference coordinate)
    end: int  # position of last mutated base
    count: int  # number of mutated bases
    start_unmodified: int  # position of last unmodified base preceding first modified base
    end_unmodified: int  # position of first unmodified base following last modified base
    record: pysam.AlignedSegment
    query_bases: list[str]
    base_qualities: list[int]
    is_revcomp: bool

    def region_tuple(self):
        return (self.record.reference_name, self.start, self.end)

    def bed_record(self, error_rate: float) -> str:
        """Format as BED record"""
        record = self.record
        base_qual = ",".join(str(q) for q in self.base_qualities)
        bases = "".join(b if b is not None else "-" for b in self.query_bases)
        number_passes = f"number_passes={record.get_tag('np')};" if record.has_tag("np") else ""
        formatted_tags = (
            f"Name={bases}{'/rc' if self.is_revcomp else ''};"
            f"read_name={record.query_name};"
            f"mapping_quality={record.mapping_quality};"
            f"bases={bases};"
            f"base_qualities={base_qual};"
            + number_passes
            + f"error_rate={error_rate:.2%}25;"  # %25 is a URL-encoded '%'
            f"mutated_region={record.reference_name}:{self.start + 1}-{self.end};"
            f"reverse_complement={'yes' if self.is_revcomp else 'no'}"
        )
        if self.is_revcomp:
            start, end = self.end, self.end_unmodified
        else:
            start, end = self.start_unmodified, self.start

        assert start <= end
        return f"{self.record.reference_name}\t{start}\t{end}\t{formatted_tags}"

    def count_mismatches(self) -> int:
        return self.count - self.query_bases.count(None)


class Statistics:
    def __init__(self):
        # Record statistics
        self.unfiltered_records: int = 0
        self.filtered_min_passes: int = 0
        self.filtered_max_error_rate: int = 0
        self.records: int = 0
        self.event_counts = [0, 0, 0]  # zero, one, two or more

        # Event filtering statistics
        self.unfiltered_events: int = 0
        self.events: int = 0  # Final number of events
        self.filtered_min_affected: int = 0
        self.filtered_min_mismatches: int = 0
        self.filtered_min_base_quality: int = 0


def main():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        "--region",
        "-r",
        help="Only work on reads within this region",
    )
    parser.add_argument(
        "--min-passes",
        metavar="N",
        type=int,
        help="Skip HiFi reads with fewer than N passes (requires 'np' tag)",
    )
    parser.add_argument(
        "--max-error-rate",
        "-e",
        metavar="RATE",
        type=float,
        default=0.1,
        help="Skip reads with error rate higher than RATE",
    )
    parser.add_argument(
        "--min-affected",
        "-N",
        metavar="N",
        type=int,
        default=5,
        help="Require N affected bases (mismatch or deletion)",
    )
    parser.add_argument(
        "--min-mismatches",
        "-n",
        metavar="N",
        type=int,
        default=3,
        help="Require at least N mismatching (not deleted) bases",
    )
    parser.add_argument(
        "--min-base-quality",
        "-b",
        metavar="QUAL",
        type=float,
        default=10,
        help="Require at least QUAL average base quality for mismatching bases",
    )
    parser.add_argument("ref", metavar="fasta", help="Indexed reference FASTA")
    parser.add_argument("bam")
    args = parser.parse_args()
    run(**vars(args))


def run(
    ref: str,
    bam: str,
    region: str | None,
    min_passes: int | None,
    max_error_rate: float | None,
    min_affected: int,
    min_mismatches: int,
    min_base_quality: float,
):
    if min_mismatches > min_affected:
        raise CommandlineError("--min-mismatches must not be larger than --min-affected")
    with open(ref) as f:
        if f.read(1) != ">":
            raise CommandlineError(
                f"file '{ref}' does not appear to be a FASTA file "
                f"as it does not start with the character '>'."
            )
    with (
        Fasta(ref) as fasta,
        pysam.AlignmentFile(bam) as af,
    ):
        stats = Statistics()

        print("#gffTags")
        for record in af.fetch(region=region):
            if record.is_secondary or record.is_supplementary or record.is_unmapped:
                continue
            stats.unfiltered_records += 1
            if (
                min_passes is not None
                and record.has_tag("np")
                and record.get_tag("np") < min_passes
            ):
                stats.filtered_min_passes += 1
                continue

            ref_record = fasta[record.reference_name]
            reference_sequence = ref_record[
                record.reference_start : record.reference_end
            ].seq.upper()
            errors = len(alignment_error_tuples(record, reference_sequence))
            error_rate = errors / len(reference_sequence)
            if max_error_rate is not None and error_rate > max_error_rate:
                stats.filtered_max_error_rate += 1
                continue

            for revcomp in (False, True):
                events = list(detect_break(record, reference_sequence, revcomp))
                stats.event_counts[min(len(events), 2)] += 1

                for event in events:
                    stats.unfiltered_events += 1
                    if event.count < min_affected:
                        stats.filtered_min_affected += 1
                        continue
                    if event.count_mismatches() < min_mismatches:
                        stats.filtered_min_mismatches += 1
                        continue
                    if mean([bq for bq in event.base_qualities if bq is not None]) < min_base_quality:
                        stats.filtered_min_base_quality += 1
                        continue
                    stats.events += 1
                    print(event.bed_record(error_rate))
            stats.records += 1

    def log(n, *args, **kwargs):
        print(f"{n:7}", *args, **kwargs, file=sys.stderr)

    print("Discarding reads with fewer than", min_passes if min_passes is not None else 0, "passes (np tag)", file=sys.stderr)
    print("Discarding reads with error rate higher than", max_error_rate, file=sys.stderr)
    print("Discarding events with fewer than", min_affected, "consecutive affected bases", file=sys.stderr)
    print("Discarding events with fewer than", min_mismatches, "mismatching bases", file=sys.stderr)

    print(file=sys.stderr)
    log(stats.unfiltered_records, "reads (alignments) in input file")
    log(stats.filtered_min_passes, "reads with too few passes (np tag) filtered out")
    log(stats.filtered_max_error_rate, "reads with too high error rate filtered out")
    log(stats.records, "reads remained after filtering and were analyzed for events")

    print(file=sys.stderr)
    log(stats.event_counts[0], "reads had no event")
    log(stats.event_counts[1], "reads had one event")
    log(stats.event_counts[2], "reads had two or more events")

    print(file=sys.stderr)
    log(stats.unfiltered_events, "events found")
    log(stats.filtered_min_affected, "events filtered because they had too few consecutive affected bases")
    log(stats.filtered_min_mismatches, "events filtered because they had too many deletions")
    log(stats.filtered_min_base_quality, "events filtered because they had too low average base quality")
    log(stats.events, "events reported after filtering")


def detect_break(
    record: pysam.AlignedSegment, reference_sequence: str, revcomp: bool
) -> Iterator[BreakEvent]:
    """Detect ssDNA break events on a single AlignedSegment"""

    base = "T" if revcomp else "A"
    mutated = False  # are we currently in a mutated region?
    prev_position = record.reference_start
    event = None
    for query_pos, ref_pos in aligned_pairs_without_softclips(record):
        if ref_pos is None:
            # Insertion
            continue

        reference_base = reference_sequence[ref_pos - record.reference_start]
        if reference_base != base:
            continue

        query_base = None if query_pos is None else record.query_sequence[query_pos]
        base_quality = None if query_pos is None else record.query_qualities[query_pos]
        cur_mutated = query_base != base
        if not mutated:
            if cur_mutated:
                # New break starts
                event = BreakEvent(
                    start=ref_pos,
                    end=ref_pos + 1,
                    count=1,
                    start_unmodified=prev_position,
                    end_unmodified=record.reference_end,  # fixed later
                    record=record,
                    query_bases=[query_base],
                    base_qualities=[base_quality],
                    is_revcomp=revcomp,
                )
        else:
            if cur_mutated:
                # Extend detected region
                event.end = ref_pos + 1
                event.count += 1
                event.query_bases.append(query_base)
                event.base_qualities.append(base_quality)
            else:
                # Unmutated replacement base: End of detected region
                pos = reference_sequence.find(base, event.end - record.reference_start)
                if pos != -1:
                    event.end_unmodified = record.reference_start + pos
                yield event
        mutated = cur_mutated
        prev_position = ref_pos
    if mutated:
        yield event


def alignment_error_tuples(
    record: pysam.AlignedSegment, reference_sequence: str
) -> list[tuple[int, int, str, str]]:
    """
    Similar to get_aligned_pairs(), but excludes soft-clipped positions and
    positions where query and reference are identical, and returns a tuple
    (query_pos, ref_pos, query_base, ref_base)
    """
    result = []
    for query_pos, ref_pos in aligned_pairs_without_softclips(record):
        ref_base = (
            None
            if ref_pos is None
            else reference_sequence[ref_pos - record.reference_start]
        )
        query_base = None if query_pos is None else record.query_sequence[query_pos]
        if ref_base != query_base:
            result.append((query_pos, ref_pos, query_base, ref_base))
    return result


def aligned_pairs_without_softclips(record: pysam.AlignedSegment):
    clip_start, clip_end_length = soft_clip_lengths(record.cigartuples)
    aligned_pairs = record.get_aligned_pairs()[clip_start:]
    if clip_end_length > 0:
        aligned_pairs = aligned_pairs[:-clip_end_length]
    return aligned_pairs


def soft_clip_lengths(cigar: list[tuple[int, int]]) -> tuple[int, int]:
    """
    >>> soft_clip_lengths([(4, 99), (0, 5), (4, 22)])
    (99, 22)
    >>> soft_clip_lengths([(4, 99)])
    (99, 0)
    >>> soft_clip_lengths([(0, 5), (4, 22)])
    (0, 22)
    >>> soft_clip_lengths([(0, 5)])
    (0, 0)
    """
    clip_start, clip_end_length = 0, 0
    op, length = cigar[0]
    if op == 4:  # S
        clip_start = length
    if len(cigar) > 1:
        op, length = cigar[-1]
        if op == 4:
            clip_end_length = length
    return clip_start, clip_end_length


def median(values):
    n = len(values)
    assert n > 0
    if n % 2 != 0:
        return sorted(values)[n // 2]
    else:
        return sum(sorted(values)[n // 2 - 1 : n // 2 + 1]) / 2

    #assert median([1]) == 1
    #assert median([2,3]) == 2.5
    #assert median([3,2,3,2]) == 2.5
    #assert median([10,3,10,5,3]) == 5


def mean(values) -> float:
    if not values:
        return 0
    return sum(values) / len(values)


if __name__ == "__main__":
    main()
