rule multimappers:
    """
    csem reassigns mapq to 100 for unimappers (also ZW tag==1)
    """
    input:
        "results/csem_mosaics/csem-dedup/{sample}.sorted.bam"
    output:
        "results/csem_mosaics/multimapper_stats/{sample}.txt"
    conda:
        "../envs/utils.yaml"
    resources:
        time=120,
        mem=24000,
        cpus=12
    threads: 12
    shell:
        """
        echo -n "sample\tread_type\tn\n{wildcards.sample}\tunimapper\t" > {output}
        printf '%s\n' $(samtools view -@ {threads} -F 0x4 {input} | awk 'int($5) == 100' | cut -f 1 | sort -u | wc -l) >> {output}
        echo -n "{wildcards.sample}\tmultimapper\t" >> {output}
        printf '%s\n' $(samtools view -@ {threads} -F 0x4 {input} | awk 'int($5) < 100' | cut -f 1 | sort -u | wc -l) >> {output}
        """

rule collect_multimappers:
    """
    collect multimapper stats
    """
    input:
        expand("results/csem_mosaics/multimapper_stats/{sample}.txt", sample=SAMPLES)
    output:
        "results/csem_mosaics/multimapper_stats/csem_multimapper_stats.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print $0"\t"FILENAME}}' {input} > {output}
        """

rule pan_te_ranking:
    """
    note that this rule's script is not complete - will need to integrate this to run with inputs, outputs, etc
    passed via snakemake, currently just to test out an analysis with the aid of the cluster
    """
    input:
        bw = "results/csem_mosaics/viz/{sample}.log2ratio.bw"
    output:
        tsv = "results/csem_mosaics/misc/pan_te_ranking/{sample}.pan_te_ranking.tsv"
    conda:
        "../envs/rutils.yaml"
    resources:
        time=40,
        mem=32000,
        cpus=1
    script:
        "../scripts/pan_te_ranking.R"


rule collect_pan_te_ranking:
    """
    """
    input:
        expand("results/csem_mosaics/misc/pan_te_ranking/{s}.pan_te_ranking.tsv", s=CHIPS)
    output:
        "results/csem_mosaics/misc/pan_te_ranking/all.pan_te_ranking.tsv.gz"
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print $0"\t"FILENAME}}' {input} | gzip -c > {output}
        """
