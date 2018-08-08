# Snakemake ATAC-seq pipeline

The standard Parker Lab ATAC-seq pipeline in Snakemake (for paired-end data).
Fastq file naming scheme should be '\*.1.fastq.gz' and '\*.2.fastq.gz'. By
default, will work with the following genomes:

1. hg19
2. hg38
3. mm9
4. mm10
5. rn4
6. rn5 

This can be changed by adding the desired genome's information to the `#GENERIC DATA`
section of the Snakefile (although ataqv may fail to run for organisms
besides fly, human, mouse, rat, worm, or yeast -- if you are processing data
from another organism, you will need to edit the pipeline to supply ataqv with
an autosomal reference file).

## Dependencies

Python >=2.7, and the following software packages:

1. fastqc
2. cta (can be downloaded from the Parker Lab github)
3. BWA
4. picard
5. samtools
6. macs2
7. bedtools
8. ataqv

Also, assumes that picard MarkDuplicates can be called using the syntax: `picard
MarkDuplicates ...`.

## Usage:

This Snakemake pipeline requires a config file (YAML format) with the following
information:

```yaml
blacklist:
    hg19:
    - /lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
    - /lab/data/reference/human/hg19/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz
bwa_index:
    hg19: /lab/data/reference/human/hg19/index/bwa/current/hg19
    mm9: /lab/data/reference/mouse/mm9/index/bwa/current/mm9
libraries:
    '100474___2156':
        genome: hg19
        readgroups:
            100474___L1___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L001.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L001.2.fastq.gz
            100474___L2___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L002.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L002.2.fastq.gz
            100474___L3___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L003.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L003.2.fastq.gz
            100474___L4___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L004.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100474_L004.2.fastq.gz
    '100477___2156':
        genome: hg19
        readgroups:
            100477___L1___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L001.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L001.2.fastq.gz
            100477___L2___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L002.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L002.2.fastq.gz
            100477___L3___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L003.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L003.2.fastq.gz
            100477___L4___2156:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L004.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477_L004.2.fastq.gz
results: /lab/work/porchard/atacseq
tss:
    hg19: /lab/data/reference/human/hg19/annot/hg19.tss.refseq.bed
    rn5: /lab/data/reference/rat/rn5/annot/rn5.tss.refseq.bed
whitelist:
    rn5: /lab/data/reference/rat/rn5/annot/rn5.K30.mappable_only.bed.gz
```

### Configuration
The basename for each fastq file must be unique, as must be the readgroup names
(keys). In many cases the only information that will be changing between
ATAC-seq experiments is the library information and the desired output directory
(paths to BWA indices, blacklists, etc. will remain unchanged). It may therefore
by convenient to have a single permanent JSON file with all of the required
information except the library information and the results dir.

If this is the case, you can use the python script at
`src/make_atacseq_config.py` to add library information and the results path to
this unchanging JSON:

```bash
python /path/to/atacseq_snakemake/bin/make_atacseq_config.py \
      -r /path/to/results_dir                                \
      -ref atacseq.generic_data.yaml                         \
      -lib libraries.yaml > complete_atacseq_config.yaml
```

Example files are given given in `examples/`.

### Running Snakemake

Snakemake would then require following files:

1. Snakefile (with your workflow)
2. Configuration (with your library information generated in previous step)
3. (Optionally) Cluster configuration if running on Cluster

One default cluster configuration for the pipeline is `src/cluster.config`. In
this mode, Snakemake outputs all job-logs to `logs/` dir, which must be created
before running Snakemake.

The pipeline can then be run:

```bash
$ mkdir -p logs
$ nohup snakemake -p -j 10 \
        --configfile /path/to/complete_atacseq_config.yaml \
        --snakefile /path/to/Snakefile \
        --cluster-config /path/to/cluster.config \
        --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" \
        &> snakemake.nohup &
```

### Dry-run

It recommended to perform a dry-run of Snakemake to see if it properly generates
all the files in workflow before submitting to cluster. For a dry-run, add `-n`
flag to the command above.
