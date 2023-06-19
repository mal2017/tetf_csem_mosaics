rule get_fq:
    """
    assumes provided accession code is a SRA run accession - e.g. a single fastq file
    also assumes SE for now, as it only takes the first fq file

    alternative when an accession: curl $(ffq --ftp {params.acc} | jq '.[0] | .url' | head -n 1 | tr -d '"') > {output}
    """
    output:
        temp("results/csem_mosaics/fastq/{sample}.fastq.gz")
    params:
        acc = lambda wc: list(set(pep.sample_table[pep.sample_table.sample_name == wc.sample].accession))
    retries: 2
    resources:
        time=20,
        mem=24000,
        cpus=1
    shell:
        """
        curl {params.acc} > {output}
        """

rule fastp:
    """
    can test that inputs are phred33 via bbtools: testformat.sh in=results/fastq/embryo_0_8hr_pan_1.fastq.gz
    in general this is true for anything pulled from SRA, I think old submissions are updated to phred33, but still worth checking
    """
    input:
        rules.get_fq.output
    output:
        fq = temp("results/csem_mosaics/trimmed/{sample}.fastp.fastq.gz"),
        h = "results/csem_mosaics/trimmed/{sample}.fastp.html",
        j = "results/csem_mosaics/trimmed/{sample}.fastp.json",
    singularity:
        "docker://quay.io/biocontainers/fastp:0.23.2--h5f740d0_3"
    params:
        max_len = config.get("MAX_READ_LEN"),
    threads:
        4
    resources:
        time=20,
        mem=24000,
        cpus=4
    shell:
        """
        fastp --in1 {input} --out1 {output} -w {threads} -b {params.max_len} -h {output.h} -j {output.j} 
        """

rule get_genome:
    output:
        "results/csem_mosaics/reference/genome.fa.gz"
    params:
        url = config.get("GENOME_URL")
    retries: 2
    shell:
        """
        curl -o {output} {params.url}
        """