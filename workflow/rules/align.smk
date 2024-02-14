rule bt1_index:
    input:
        ref =rules.get_genome.output
    output:
        directory('results/csem_mosaics/indices/bt1/')
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/bowtie:1.3.1--py310h4070885_4"
    resources:
        time=40,
        mem=24000,
        cpus=12
    shell:
        """
        mkdir -p {output} &&
        cp {input.ref} {output} &&
        bowtie-build --threads {threads} {input.ref} {output}/genome
        """


rule bt1_align:
    """
    https://deweylab.biostat.wisc.edu/csem/
    trying to mimic the bt1 approach used for the CSEM paper
    `-q -v 2 -a -m 99 -p 8 -S mm9 GATA1.fq GATA1.sam`
    v is mismatches
    m is maximum number of alignments for a read to be included at all
    a is to report all alignments
    """
    input:
        r = rules.fastp.output.fq,
        idx = rules.bt1_index.output,
    output:
        al = temp("results/csem_mosaics/alignments/{sample}.sam")
    threads:
        12
    log:
        "results/csem_mosaics/logs/bt1_align/{sample}.log"
    singularity:
        "docker://quay.io/biocontainers/bowtie:1.3.1--py310h4070885_4"
    resources:
        time=40,
        mem=48000,
        cpus=24
    shell:
        """
        bowtie --threads {threads} -q -v 2 -a -m 99 -p 8 --sam -x {input.idx}/genome {input.r} {output.al} 2> {log}
        """


rule picard_sort_by_name:
    input:
        rules.bt1_align.output.al
    output:
        bam = temp("results/csem_mosaics/alignments/{sample}.nsort.bam"),
    singularity:
        "docker://quay.io/biocontainers/picard:3.0.0--hdfd78af_1"
    resources:
        time=120,
        mem=24000,
        cpus=1
    shell:
        """
        picard SortSam -I {input} -O {output.bam} -SO queryname 
        """

rule picard_markdup:
    """
    namesorted input is important so that multiple alignments of the same read are adjacent and therefore removable
    see picard markduplicates docs for more details
    """
    input:
        rules.picard_sort_by_name.output.bam
    output:
        bam = temp("results/csem_mosaics/alignments/{sample}.markdup.bam"),
        stats = "results/csem_mosaics/alignments/{sample}.markdup.stats"
    singularity:
        "docker://quay.io/biocontainers/picard:3.0.0--hdfd78af_1"
    resources:
        time=40,
        mem=24000,
        cpus=1
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {input} -M {output.stats} -O {output.bam} 
        """