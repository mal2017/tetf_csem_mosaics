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
    wei et al. 2022 https://doi.org/10.1093/molbev/msac080
    """
    input:
        ip = "results/csem_mosaics/csem/{sample}.bw",
        wce = lambda wc: expand("results/csem_mosaics/csem/{c}.bw", c=pep.get_sample(wc.sample).input),
    output:
        bw = temp("results/csem_mosaics/viz/{sample}.scaled.bw"),
        bw_wce = temp("results/csem_mosaics/viz/{sample}-wce.scaled.bw"),
    conda:
        "../envs/mosaics.yaml"
    resources:
        time=20,
        mem=24000,
        cpus=1
    script:
        "../scripts/scale_csem_bigwig.R"

rule make_log2f_bw:
    """
    wei et al. 2022 https://doi.org/10.1093/molbev/msac080
    """
    input:
        ip = rules.scale_csem_bigwig.output.bw,
        wce = rules.scale_csem_bigwig.output.bw_wce,
        #wce = lambda wc: "results/csem_mosaics/csem/{c}.bw".format(c=pep.get_sample(wc.sample).input),
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
        bigwigCompare -b1 {input.ip} -b2 {input.wce} --pseudocount 1 1 --skipZeroOverZero -o {output}
        """

rule plotfingerprint:
    input:
        chip_bam = "results/csem_mosaics/csem-dedup/{sample}.sorted.bam",
        chip_bai = "results/csem_mosaics/csem-dedup/{sample}.sorted.bam.bai",
        wce_bam = lambda wc: expand("results/csem_mosaics/csem-dedup/{c}.sorted.bam", c=pep.get_sample(wc.sample).input),
        wce_bai = lambda wc: expand("results/csem_mosaics/csem-dedup/{c}.sorted.bam.bai", c=pep.get_sample(wc.sample).input),
    output:
        png = "results/csem_mosaics/viz/{sample}.fingerprint.png",
        metrics = "results/csem_mosaics/viz/{sample}.fingerprint.metrics.txt",
    singularity:
        "docker://quay.io/biocontainers/deeptools:3.5.2--pyhdfd78af_1"
    threads: 16
    resources:
        time=20,
        mem=24000,
        cpus=16
    shell:
        """
        plotFingerprint -p {threads} -b {input.chip_bam} {input.wce_bam} --JSDsample {input.wce_bam} --outQualityMetrics {output.metrics} -plot {output.png}
        """
