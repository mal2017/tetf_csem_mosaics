rule csem_bam2wig:
    input:
        bam = rules.csem.output.bam,
    output:
        wig = temp("results/csem_mosaics/csem/{sample}.wig"),
    params:
        path = config.get("CSEM_DIR")
    resources:
        time=20,
        mem=24000,
        cpus=1
    shell:
        """
        {params.path}/csem-bam2wig {input.bam} {output.wig} CSEM_CHIP --extend-reads 200
        """

rule fa_and_fai:
    input:
        fa = "results/csem_mosaics/reference/genome.fa.gz",
    output:
        fa = temp("results/csem_mosaics/viz/genome.fa"),
        fai = "results/csem_mosaics/viz/genome.fa.fai",
    conda:
        "../envs/utils.yaml"
    shell:
        """
        gunzip -c {input.fa} > {output.fa} &&
        samtools faidx {output.fa}
        """

rule wigtobigwig:
    input:
        wig = rules.csem_bam2wig.output.wig,
        chrom_sizes = rules.fa_and_fai.output.fai,
    output:
        bigwig = temp("results/csem_mosaics/csem/{sample}.bw"),
    singularity:
        "docker://quay.io/biocontainers/ucsc-wigtobigwig:377--ha8a8165_3"
    resources:
        time=20,
        mem=24000,
        cpus=1
    shell:
        """
        wigToBigWig {input.wig} {input.chrom_sizes} {output.bigwig}
         """

rule scale_csem_bigwig:
    """
    currently scales the ip to have the same total value as the wce
    """
    input:
        ip = "results/csem_mosaics/csem/{sample}.bw",
        wce = lambda wc: "results/csem_mosaics/csem/{c}.bw".format(c=pep.get_sample(wc.sample).input),
    output:
        bw = temp("results/csem_mosaics/viz/{sample}.scaled.bw"),
    conda:
        "../envs/mosaics.yaml"
    resources:
        time=20,
        mem=24000,
        cpus=1
    script:
        "../scripts/scale_csem_bigwig.R"

rule make_log2f_bw:
    input:
        ip = rules.scale_csem_bigwig.output.bw,
        wce = lambda wc: "results/csem_mosaics/csem/{c}.bw".format(c=pep.get_sample(wc.sample).input),
    output:
        log2f = "results/csem_mosaics/viz/{sample}.log2ratio.bw",
    singularity:
        "docker://quay.io/biocontainers/deeptools:3.5.2--pyhdfd78af_1"
    resources:
        time=20,
        mem=24000,
        cpus=1
    shell:
        """
        bigwigCompare -b1 {input.ip} -b2 {input.wce} --skipZeroOverZero -o {output}
        """