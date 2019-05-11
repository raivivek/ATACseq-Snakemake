##
## Process Samples
##
#:
#: Samples are often run as multiple libraries. We merge libraries and analyze them
#: as single samples. If any libraries are removed from the analysis after QC,
#: the rules below should be be updated to reflect those omissions.
#:

rule make_samples:
    input:
        expand(
            _samples("ataqv", "{sample}.ataqv.json.gz"),
            sample=iterate_all_samples()
        ),
        expand(_samples("gb", "{sample}.bw"), sample=iterate_all_samples()),


rule merge_libraries:
    input:
        lambda wildcards: expand(
            _libraries("pruned", "{library}.pd.bam"),
            library=iterate_sample_libraries(wildcards.sample)
        )
    output:
        bam = _samples("merge_libraries", "{sample}.bam"),
        bai = _samples("merge_libraries", "{sample}.bam.bai")
    resources: io_concurrent = 1
    threads: 2
    shell:
        """
        samtools merge -@{threads} - {input} \
            | samtools sort -m 4G -O bam -o {output.bam} - 2> /dev/null

        samtools index -@{threads} {output.bam}

        samtools flagstat -@{threads} {output.bam} > {output.bam}.flagstat
        """


rule sample_bam2bed:
    input:
        _samples("merge_libraries", "{sample}.bam")
    output:
        _samples("bam2bed", "{sample}.bed")
    shell:
        """bedtools bamtobed -i {input} > {output}"""


rule sample_peaks:
    input:
        _samples("bam2bed", "{sample}.bed")
    output:
        peaks = _samples("macs2", "{sample}_peaks.broadPeak"),
        bdg = _samples("macs2", "{sample}_treat_pileup.bdg.gz"),
    params:
        name = "{sample}",
        outdir = _samples("macs2"),
        genome_size = lambda wildcards: MACS2_GENOME_SIZE[get_sample_genome(wildcards.sample)]
    log:
        _log("{sample}.macs2.out")
    shell:
        """
        macs2 callpeak \
          --outdir {params.outdir} \
          -t {input} \
          -n {params.name} \
          -f BED \
          -g {params.genome_size} \
          --nomodel \
          --shift -100 \
          --extsize 200 \
          --seed 2018 \
          -B \
          --broad \
          --keep-dup all \
          --SPMR \
          &> {log}

        pigz {params.outdir}/{wildcards.sample}_treat_pileup.bdg \
                {params.outdir}/{wildcards.sample}_control_lambda.bdg
        """


rule sample_noblacklist:
    input:
        _samples("macs2", "{sample}_peaks.broadPeak")
    output:
        noblacklist = _samples("macs2", "{sample}_peaks.broadPeak.noblacklist"),
        fdr05 = _samples("macs2", "{sample}_peaks.fdr05.bed")
    params:
        blacklists = lambda wildcards: " ".join(
            get_blacklists(get_sample_genome(wildcards.sample))
        ),
        fdr = config["params"].get("macs2_fdr", 0.05),
    shell:
        """mappabilityFilter -i {input} -b {params.blacklists} \
            | tee {output.noblacklist} \
            | createMasterPeaks --fdr {params.fdr} > {output.fdr05}
        """


rule sample_ataqv:
    """Ataqv-toolkit is a ATAC-seq experiment QC tool developed in Parker Lab.
    The tools provides many useful metrics such as Fragment Length
    Distribution, TSS Enrichment, for quality comparison and downstream
    analysis.

    See: https://github.com/ParkerLab/ataqv
    """
    input:
        bam = _samples("merge_libraries", "{sample}.bam"),
        peaks = _samples("macs2", "{sample}_peaks.fdr05.bed")
    output:
        _samples("ataqv", "{sample}.ataqv.json.gz")
    params:
        name = "{sample}",
        description = "{sample}",
        organism = lambda wildcards: get_organism(get_sample_genome(wildcards.sample)),
        tss = lambda wildcards: get_tss(get_sample_genome(wildcards.sample))
    log:
        _log("{sample}.ataqv.log")
    shell:
        """
        ataqv --peak-file {input.peaks}    \
          --name {params.description}      \
          --metrics-file {output}          \
          --tss-file {params.tss}          \
          --ignore-read-groups             \
          {params.organism}                \
          {input.bam}                      \
          > {log}
        """


rule sample_bigwig:
    input:
        _samples("macs2", "{sample}_treat_pileup.bdg.gz"),
    output:
        _samples("gb", "{sample}.bw")
    params:
        sizes = lambda wildcards: get_chrom_sizes(get_sample_genome(wildcards.sample))
    shell:
        """bdgTobw {input} {params.sizes} {output}"""


rule sample_viewer:
    input:
        expand(_samples("ataqv", "{sample}.ataqv.json.gz"),
                sample=iterate_all_samples)
    output:
        directory(_samples("viewer"))
    params:
        description = config["params"].get("description", "Samples")
    shell:
        """mkarv -d {params.description} {output} {input}"""


# vim:syntax=python
