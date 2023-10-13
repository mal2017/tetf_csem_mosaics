localrules: unzip_fa_for_motifs, get_repeat_bed

rule get_repeat_bed:
    input:
        config.get("repeat_bed")
    output:
        "results/csem_mosaics/repeats.bed"
    shell:
        """
        cp {input} {output}
        """


rule unzip_fa_for_motifs:
    input:
        fa = rules.get_genome.output
    output:
        fa = "results/csem_mosaics/motifs/genome.fa"
    shell:
        """
        gunzip -c {input.fa} > {output.fa}
        """

rule mask_fasta:
    """
    we want to mask simple repeats for this analysis because I've found it tends to 
    swamp xstreme with abundant low complexity hits, so we use the pattern that matches repeatmasker's simple repeat hits
    """
    input:
        fa = rules.unzip_fa_for_motifs.output.fa,
        bed = rules.get_repeat_bed.output,
    output:
        fa = "results/csem_mosaics/motifs/simple_reps_masked_genome.fa"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2"
    shell:
        """
        grep ')n' {input.bed} | bedtools maskfasta -fi {input.fa} -bed stdin -fo {output.fa}
        """

rule get_fasta_from_peaks:
    """
    I also filter out peaks overlapping more complex repeat repeats entirely because I found a highly enriched motif that was part of the SAR satellite. 
    
    This wasn't included in the original library I used for masking, but is clearly worth removing for this analysis, as I intend for this to be a 'repeat-free' peak set.

    Therefore, I use both my masking and the UCSC repeatmasker track to filter out peaks overlapping repeats.

    UCSC repeat masker track from ./resources/ (see script there for generation of this track, note that this performs chrom renaming too bc UCSC chromnames don't match flybase).
    """
    input:
        fa = rules.mask_fasta.output.fa,
        bed = rules.mosaics_main_dedup.output.bed,
        reps = rules.get_repeat_bed.output,
        ucsc_reps = config.get("ucsc_repeats")
    output:
        fa= "results/csem_mosaics/motifs/peak_seqs/{sample}.non_te_peaks.fasta",
        bed= "results/csem_mosaics/motifs/peak_seqs/{sample}.non_te_peaks.bed"
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2"
    shell:
        """
        cat {input.reps} {input.ucsc_reps} | \
            grep -v ')n' | cut -f 1,2,3 | \
            bedtools intersect -v -a {input.bed} -b stdin > {output.bed} &&
        bedtools getfasta -fi {input.fa} -bed {output.bed} -fo {output.fa}
        """

rule mask_seqs:
    """
    mask simple repeats from the peak sequences
    """
    input:
        tes = rules.get_fasta_from_peaks.output.fa,
    output:
        fa =  "results/csem_mosaics/motifs/peak_seqs/{sample}.non_te_peaks.simple_mask.fasta",
    params:
        entropy = 0.7,
        w = 12,
    singularity:
        "docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    shell:
        """
        bbmask.sh \
            in={input.tes} \
            out={output.fa} \
            w={params.w} \
            entropy={params.entropy}
        """

localrules: get_pan_motifs

rule get_pan_motifs:
    """
    motifs id'd by downloading the meme docker image, starting a shell, and running the following shell command:
    `grep -E -e "FBgn0085432|pan" /opt/meme/share/meme-5.5.0/db/motif_databases/*/*.meme`

    'Archbold' refers to archbold 2014 plos genetics

    cisbp motif, flyreg, flyfactor, etc are all the same or very similar to jaspar1.1, so I don't include them here.
    """
    output:
        jaspar = "results/csem_mosaics/motifs/pan_motifs/jaspar.meme",
        jaspar2 = "results/csem_mosaics/motifs/pan_motifs/jaspar2.meme",
        archbold_hmg = "results/csem_mosaics/motifs/pan_motifs/archbold_hmg.meme",
        archbold_helper = "results/csem_mosaics/motifs/pan_motifs/archbold_helper.meme",
        combined = "results/csem_mosaics/motifs/pan_motifs/pan.meme",
    threads:
        1
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    shell:
        """
        meme-get-motif -id MA0237.1 /opt/meme/share/meme-5.5.0/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar} &&
        meme-get-motif -id MA0237.2 /opt/meme/share/meme-5.5.0/db/motif_databases/JASPAR/JASPAR2022_CORE_insects_redundant_v2.meme > {output.jaspar2} &&
        iupac2meme -dna GCCGCCR > {output.archbold_helper} &&
        iupac2meme -dna SCTTTGWSW > {output.archbold_hmg} &&
        meme2meme {output.archbold_hmg} {output.archbold_helper} {output.jaspar} {output.jaspar2} > {output.combined}
        """

rule xstreme_per_tf:
    input:
        fa = rules.mask_seqs.output.fa,
        known_motifs = rules.get_pan_motifs.output.combined,
    output:
        odir = directory("results/csem_mosaics/motifs/xstreme/{sample}/")
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    threads: 48
    resources:
        time=120,
        mem=96000,
        cpus=48
    shell:
        """
        xstreme --oc '{output.odir}' --m {input.known_motifs} \
            --maxw 10 --p {input.fa} --meme-p {threads} --order 0 --group-thresh 0.4
        """

rule sea_per_tf:
    input:
        fa = rules.mask_seqs.output.fa,
        meme = rules.get_pan_motifs.output.combined,
    output:
        odir = directory("results/csem_mosaics/motifs/sea/{sample}/")
    singularity:
        "docker://memesuite/memesuite:5.5.0"
    threads: 2
    resources:
        time=10,
        mem=12000,
        cpus=2
    shell:
        """
        sea --oc '{output.odir}' \
            --p {input.fa} \
            --m {input.meme} \
            --order 0
        """

