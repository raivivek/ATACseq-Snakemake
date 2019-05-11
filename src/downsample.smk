##
## (Optional) Perform downsampling; call `rule downsample` directly
##
rule downsample:
    """ Rule that governs downsampling of the processed BAM files.

    Downsampling is desirable when comparing many samples as sequencing-depth
    can easily factor into observed differences.
    """
    input:
        rules.make_samples.input,
        expand(
            _downsampled("macs2", "{sample}_peaks.noblacklist.bed"),
            sample=iterate_all_samples()
        ),
        _downsampled("master-peaks", "all-merged_peaks.noblacklist.bed"),


rule downsample_bam:
  input:
        rules.merge_libraries.output
  output:
        bam = _downsampled('merge_libraries', "{sample}.dwnsmpl.bam"),
        bai = _downsampled('merge_libraries', "{sample}.dwnsmpl.bam.bai"),
        flagstat = _downsampled('merge_libraries', "{sample}.dwnsmpl.bam.flagstat")
  threads: 1
  params:
        seed = config["params"].get("seed", 2018),
        depth = config["params"].get("subsample_depth")
  run:
      # if not supplied; no need to downsample, just print error and exit
      if params.depth is None:
          sys.stderr.write("ERROR: Subsampling depth not supplied; Exiting!\n")
          sys.exit(1)
      else:
          # Compute desired fraction of reads to sub-sample using samtools
          # flagstat.
          #
          # NOTE: __TOTAL_READS computes the total number of # "single-reads" in the
          # given input BAM file, and then we compute the fraction to subsample
          # in __FRACTION (requires `bc` tools, POSIX so almost always
          # available). The __FRACTION value already contains a "." (period) so
          # we omit that in samtools call. `scale=3` in `bc -l` limits the number
          # of decimal places to 3.
          #
          # Further, simply symlink files if desired depth is more then available
          # number of reads in the file.
          #
          # NOTE: (( . )) is not a typo. Using BASH"s arithmetic context.
        shell(
        """
            __TOTAL_READS=$(samtools view -cF 0x100 {input.bam}) \
            && if (( $_TOTAL_READS < {params.depth} )); then
                ln -s {input.bam} {output.bam}
                ln -s {input.bai} {output.bai}
            else
                __FRACTION=$(bc -l <<< "scale=3; {params.depth}/$__TOTAL_READS") \
                && samtools view -@ {threads} -b -h -s {params.seed}$__FRACTION {input.bam} -o {output.bam}
                samtools index {output.bam}
                samtools flagstat {output.bam} > {output.flagstat}
            fi
        """
        )


rule downsample_bam_to_bed:
    input:
        _downsampled("merge_libraries", "{sample}.dwnsmpl.bam"),
    output:
        _downsampled("bam2bed", "{sample}.dwnsmpl.bed"),
    shell:
        """bedtools bamtobed -i {input} > {output}"""


rule downsample_call_peaks:
    input:
        _downsampled("bam2bed", "{sample}.dwnsmpl.bed"),
    output:
        peaks = _downsampled("macs2", "{sample}_peaks.broadPeak"),
        bdg = _downsampled("macs2", "{sample}_treat_pileup.bdg.gz"),
    params:
        name = "{sample}",
        outdir = _downsampled("macs2"),
        genome_size = lambda wildcards: MACS2_GENOME_SIZE[get_sample_genome(wildcards.sample)]
    log:
        _log("{sample}.dwnsmpl.macs2.out")
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


rule downsample_noblacklist:
    input:
        _downsampled("macs2", "{sample}_peaks.broadPeak")
    output:
        _downsampled("macs2", "{sample}_peaks.noblacklist.bed")
    params:
        blacklists = lambda wildcards: " ".join(
            get_blacklists(get_sample_genome(wildcards.sample))
        ),
        fdr = config["params"].get("macs2_fdr", 0.05),
    shell:
        """mappabilityFilter -i {input} -b {params.blacklists} | \
                createMasterPeaks --fdr {params.fdr} > {output}
        """

rule make_master_downsample_bed:
    input:
        expand(
            _downsampled("bam2bed", "{sample}.dwnsmpl.bed"),
            sample=iterate_all_samples()
        )
    output:
        bed = _downsampled("master-peaks", "{prefix}.bed"),
    shell:
        """cat {input} > {output.bed}"""


rule call_master_peaks:
    input:
         _downsampled("master-peaks", "{prefix}.bed"),
    output:
        peaks = _downsampled("master-peaks", "{prefix}_peaks.broadPeak"),
        bdg = _downsampled("master-peaks", "{prefix}_treat_pileup.bdg.gz"),
        bed = _downsampled("master-peaks", "{prefix}_peaks.noblacklist.bed"),
    params:
        name = "{prefix}",
        outdir = _downsampled("master-peaks"),
        fdr = config["params"].get("macs2_fdr", 0.05),
        genome_size = "hs"  # harcoded for now
    log:
        _log("{prefix}.macs2.out")
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

        pigz {params.outdir}/{wildcards.prefix}_treat_pileup.bdg \
                {params.outdir}/{wildcards.prefix}_control_lambda.bdg

        mappabilityFilter -i {output.peaks} -g hg19 | \
                createMasterPeaks --fdr {params.fdr} > {output.bed}
        """


# vim:syntax=python
