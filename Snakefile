from glob import glob

configfile: "breaks.yaml"

localrules: symlink_pacbio, index_bed, sort_bed, bgzip_bed, subtract_illumina_controls


REF = "ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

DATASETS = [f"illumina/TK-2746-{i}.filtered" for i in range(3, 9)] + \
    [f"nanopore/barcode{i:02d}" for i in range(1, 11)] + \
    [f"pacbio/ps_396_00{i}" for i in (1, 2)]


rule:
    input:
        [f"{ds}.bed.gz.tbi" for ds in DATASETS]

rule map_illumina:
    output: bam="illumina/{name}.bam"
    input:
        ref=REF,
        r1_fastq="raw/illumina/{name}_R1.fastq.gz",
        r2_fastq="raw/illumina/{name}_R2.fastq.gz",
    log:
        bam="illumina/{name}.bam.log",
        sort="illumina/{name}.bam.sort.log",
    threads: 20
    shell:
        "strobealign -t {threads} {input.ref} {input.r1_fastq} {input.r2_fastq} 2> {log.bam}"
        "| samtools view -h -e '[NM]>1'"
        "| samtools sort -@ 4 -o {output.bam}.tmp.bam 2> {log.sort}"
        " && mv {output.bam}.tmp.bam {output.bam}"


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """
    Use this function as sorting key to sort in 'natural' order.

    >>> names = ['file10', 'file1.5', 'file1']
    >>> sorted(names, key=natural_sort_key)
    ['file1', 'file1.5', 'file10']

    Source: http://stackoverflow.com/a/16090640/715090
    """
    return [int(text) if text.isdigit() else text.lower()
        for text in re.split(_nsre, s)]


rule map_nanopore:
    output: bam="nanopore/{name}.bam"
    input:
        ref=REF,
        fastq=lambda wildcards: sorted(glob(f"raw/nanopore/{wildcards.name}/*.fastq.gz"), key=natural_sort_key)
    log:
        bam="nanopore/{name}.bam.log",
        sort="nanopore/{name}.bam.sort.log",
    threads: 20
    shell:
        "fastq=nanopore/{wildcards.name}.fastq; "
        "rm -f ${{fastq}}; "
        "mkfifo ${{fastq}}; "
        "zcat {input.fastq} > ${{fastq}} &"
        "minimap2 -ax map-ont -t {threads} {input.ref} ${{fastq}} 2> {log.bam}"
        "| samtools sort -@ 4 -o {output.bam}.tmp.bam 2> {log.sort}; "
        "mv {output.bam}.tmp.bam {output.bam}; "
        "rm ${{fastq}}"

rule symlink_pacbio:
    # TODO we may want to re-map instead
    output: bam="pacbio/{name}.bam"
    input: "raw/pacbio/{name}.mapped.bam"
    shell: "ln -rs {input} {output}"

rule detect_breaks:
    output: bed=temp("{tech}/{name}.unsorted.bed")
    input:
        bam="{tech}/{name}.bam",
        bai="{tech}/{name}.bam.bai",
        ref=REF,
    log: "{tech}/{name}.breaks.log"
    params:
        conf=lambda wildcards: config[wildcards.tech]
    shell:
        "python detectbreak.py"
        " --min-affected={params.conf[min_affected]}"
        " --min-mismatches={params.conf[min_mismatches]}"
        " {input.ref}"
        " {input.bam}"
        " > {output}.tmp"
        " 2> {log}"
        " && mv {output.bed}.tmp {output.bed}"

rule sort_bed:
    output: bed="{name}.bed.gz"
    input: bed="{name}.unsorted.bed"
    shell:
        "( echo '#gffTags' ;"
        " sort -k1,1 -k2,2n -u {input.bed} | uniq"
        ")"
        " | bgzip > {output.bed}.tmp"
        " && mv {output.bed}.tmp {output.bed}"

rule subtract_illumina_controls:
    output: bed=temp("illumina/{name}.filtered.bed")
    input:
        bed="illumina/{name}.bed.gz",
        control1="illumina/TK-2746-1.bed.gz",
        control2="illumina/TK-2746-2.bed.gz",
    shell:
        "bedtools subtract -A -header -a {input.bed} -b {input.control1} > {output.bed}.tmp1.bed; "
        "bedtools subtract -A -header -a {output.bed}.tmp1.bed -b {input.control2} > {output.bed}.tmp2.bed; "
        "rm {output.bed}.tmp1.bed; "
        "mv {output.bed}.tmp2.bed {output.bed}"

rule index_bed:
    output: tbi="{name}.bed.gz.tbi"
    input: bed="{name}.bed.gz"
    shell:
        "tabix -f {input.bed}"

rule index_bam:
    output: "{name}.bam.bai"
    input: "{name}.bam"
    shell: "samtools index {input}"

rule bgzip_bed:
   output: bed="{name}.bed.gz"
   input: bed="{name}.bed"
   shell: "bgzip < {input.bed} > {output.bed}"
