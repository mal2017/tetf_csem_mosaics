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
        time=40,
        mem=24000,
        cpus=8
    shell:
        """
        {params.path}/run-csem --bam \
              -p {threads} \
              {input} \
              200 \
              $(dirname {output.bam})/$(basename -s .sorted.bam {output.bam})
        """

ruleorder: addrg > csem

rule addrg:
    input:
        rules.csem.output.bam
    output:
        bam = "results/csem_mosaics/csem-dedup/{sample}.sorted.rg.bam"
    resources:
        time=40,
        mem=24000,
        cpus=1
    shell:
        """
        picard AddOrReplaceReadGroups \
            CREATE_INDEX=true \
            I={input} \
            O={output.bam} \
            RGID={wildcards.sample} \
            RGLB={wildcards.sample} \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM={wildcards.sample}
        """

rule csem_generate_bed:
    """
    MOSAiCS expected format
    """
    input:
        bam = rules.addrg.output.bam
    output:
        bed = "results/csem_mosaics/csem-dedup/{sample}.bed"
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