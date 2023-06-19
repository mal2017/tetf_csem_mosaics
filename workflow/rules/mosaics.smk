

rule mosaics_bin_level:
    input:
        bed = rules.csem_generate_bed.output,
        sizes = config.get("CHROM_SIZES")
    output:
        dir = directory("results/csem_mosaics/mosaics-bin-level/{sample}/")
    params:
        fraglen = FRAGLEN,
        binsize = BINSIZE,
    conda:
        "../envs/mosaics.yaml"
    resources:
        time=40,
        mem=24000,
        cpus=1
    script:
        "../scripts/mosaics-make-bin-level-file.R"

rule mosaics_main_dedup:
    input:
        chip = "results/csem_mosaics/mosaics-bin-level/{sample}/",
        control = lambda wc: "results/csem_mosaics/mosaics-bin-level/{c}/".format(c=pep.get_sample(wc.sample).input),
        GC = config.get("GC_CONTENT"),
        N= config.get("N_CONTENT"),
        M= config.get("MAPPABILITY"),
    output:
        bed = "results/csem_mosaics/mosaics/{sample}/{sample}.mosaics.bed",
        rds = "results/csem_mosaics/mosaics/{sample}/{sample}.mosaics.rds"
    params:
        chip = lambda wc: expand("results/csem_mosaics/mosaics-bin-level/{c}/{c}.bed_fragL{fl}_bin{bs}.txt",fl=FRAGLEN, bs=BINSIZE, c=wc.sample),
        control = lambda wc: expand("results/csem_mosaics/mosaics-bin-level/{c}/{c}.bed_fragL{fl}_bin{bs}.txt", fl=FRAGLEN, bs=BINSIZE, c=pep.get_sample(wc.sample).input),
    conda:
        "../envs/mosaics.yaml"
    resources:
        time=40,
        mem=24000,
        cpus=1
    script:
        "../scripts/mosaics.R"