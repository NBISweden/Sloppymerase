# Single-stranded DNA break detection using an engineered error-prone DNA polymerase

This software and pipeline accompany the paper
*Sequence-templated erroneous end-labelling to map single-stranded DNA breaks,
using an engineered error-prone DNA polymerase*
by Wenson et al. (2024)

## Overview

This repository contains

- A [Snakemake](https://snakemake.readthedocs.io/) workflow for replicating the results from the paper.
- The ssDNA break detection script itself (`detectbreak.py`),
  which can be used independently of the workflow.

## Running the full workflow

1. Clone this repository
2. Create a Conda environment:

       conda env create -f environment.yml

3. Activate the Conda environment:

       conda activate sloppymerase

4. Make the raw data available in directory `raw/`.
5. Run the workflow:

       snakemake --profile slurm

## Running the break detection script

The break detection script is a single file called `detectbreak.py`.
To run it and the required postprocessing step, the following software needs
to be made available: Python, pysam, pyfaidx, tabix and bedtools.

This can be done with Conda, for example:

    conda create -n ssdna python pysam pyfaidx tabix bedtools
    conda activate ssdna

The script needs to be provided with the path to the reference FASTA file
and a BAM or CRAM file with mapped reads.
It prints a BED file to standard output:

    python detectbreak.py reference.fasta mapped.bam > tmp.bed

This will create a `tmp.bed` file and output some statistics to the terminal.

The following commands need to be run in order to be able to open the
BED file in IGV
(this will add the required header, sort the events and remove duplicates):

    ( echo '#gffTags' ; bedtools sort -i tmp.bed | uniq ) | bgzip > breaks.bed.gz
    tabix -f breaks.bed.gz
    rm tmp.bed

Now the `breaks.bed.gz` file can be opened in IGV.


# Background

The program searches the aligned reads in the input BAM file for ssDNA
breaks (called "events" below).

Since each aligned read is searched independently,
the same event may be reported multiple times
if there exist multiple reads covering it.
These duplicates can be removed by a post-processing step
(`bedtools sort` followed by `uniq`).


## Filtering criteria

- Supplementary alignments are ignored.
- Alignments whose error rate is above a threshold are discarded.
  The error rate is computed as the sum of the number of substitutions, insertions and deletions divided by the length of the reference sequence to which the read was aligned.
- Alignments where the number of passes (SAM `np` tag) is below a threshold are discarded.
  This is only relevant for PacBio HiFi reads.
- Events with too few affected bases are discarded.
- Events where too few of the affected bases are mismatches are discarded.
  (Important for Nanopore reads where otherwise events with only indels are found.)

# BED annotations

For each aligned read in the input alignment file (BAM or CRAM),
the program inspects the read alignments to find regions of mutated `A` or `T` bases in each read.

If `A` bases are mutated, the break event is determined to start just after the last non-mutated `A`
before the mutated region and to end just before the first mutated `A`.

If `T` bases are mutated, the break event is determined to start just after the last mutated `A`
and ends just before the first nonmutated `A` following it.

The BED track contains all break events. Each event has a "name" that is displayed under the event
(in IGV). The name contains the mutated bases and also indicates whether the event was detected in
the forward direction (mutated `A`) or in the reverse direction.

For example, the name `GC/rc` means that the mutated regions encompasses two bases, of which the
first was replaced with `G` and the second with `C`. The `/rc` indicates that this is an event on
the reverse complement. If a base is deleted, it is shown as `-` (dash).

Additional annotations are displayed by IGV when hovering with the mouse above an event. For
example:

    GC/rc
    chr1:20280160-20280159
    Name: GC/rc
    read_name: m54259_220609_004305/66061181/ccs
    mapping_quality: 60
    base_qualities: 93,93
    bases: GC
    number_passes: 58
    error_rate: 0.53%
    mutated_region: chr1:20280157-20280159
    reverse_complement: yes

`Name`: The "name" as described above.

`read_name`: The name of the read.

`mapping_quality`: Confidence that the mapping location is correct. A mapping quality of 0 means
that the read is a multimapper that could have been mapped to other, equally good positions on
the reference. 60 is very high.

`bases: GC`: Mutated bases. Here, the mutated region spans two `T` bases that both have
been mutated to `G`. The `GG` is even displayed directly under the event (without hovering) and
therefore allows to quickly assess

`base_qualities`: Sequencing quality of the mutated bases. Low base qualities could indicate that a
base was mutated because of a sequencing error.

`number_passes`: For PacBio HiFi reads. The number of times the sequencer read the fragment.
Each pass increases sequencing quality.
A low number of passes means that each sequenced base has a higher error rate.

`error_rate`: The error rate for the entire read (this currently includes even detected mutated
regions.)

`mutated_region`: The detected mutated region.
