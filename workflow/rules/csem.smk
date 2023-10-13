rule csem:
    """
    https://deweylab.biostat.wisc.edu/csem/
    """
    input:
        rules.picard_markdup.output.bam,
    output:
        uns = temp("results/csem_mosaics/csem-dedup/{sample}.bam"),
        bam = temp("results/csem_mosaics/csem-dedup/{sample}.sorted.bam"),
        bai = temp("results/csem_mosaics/csem-dedup/{sample}.sorted.bam.bai")
    threads:
        8
    params:
        path = config.get("CSEM_DIR")
    resources:
        time=240,
        mem=48000,
        cpus=48
    shell:
        """
        {params.path}/run-csem --bam \
              -p {threads} \
              {input} \
              200 \
              $(dirname {output.bam})/$(basename -s .sorted.bam {output.bam})
        """

rule csem_generate_bed:
    """
    MOSAiCS expected format
    """
    input:
        bam = rules.csem.output.bam
    output:
        bed = temp("results/csem_mosaics/csem-dedup/{sample}.bed")
    params:
        path = config.get("CSEM_DIR")
    resources:
        time=40,
        mem=24000,
        cpus=1
    shell:
        """
        {params.path}/csem-generate-input --bed {input.bam} $(dirname {output.bed})/$(basename -s .bed {output.bed})
        """
