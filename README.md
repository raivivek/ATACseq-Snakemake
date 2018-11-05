# Snakemake ATAC-seq pipeline

The standard Parker Lab ATAC-seq pipeline in Snakemake (for paired-end data).
Fastq file naming scheme should be `*.1.fastq.gz` and `*.2.fastq.gz`. By
default, will work with the following genomes:

1. hg19
2. hg38
3. mm9
4. mm10
5. rn4
6. rn5 

This can be changed by adding the desired genome's information to
the `#GENERIC DATA` section of the Snakefile (although ataqv may fail to run
for organisms besides fly, human, mouse, rat, worm, or yeast -- if you are
processing data from another organism, you will need to edit the pipeline
to supply ataqv with an autosomal reference file).

## Dependencies

Python >=3.6, and the following software packages:

1. fastqc
2. [cta](https://github.com/ParkerLab/cta)
3. BWA
4. picard (for example, `picard MarkDuplicates`)
5. samtools
6. macs2
7. bedtools
8. [ataqv](https://github.com/ParkerLab/ataqv)

## Usage

The Snakemake file utilizes the `basename` of `fastq` files to organize and
generate files, as must be the readgroup names; thus they **must** be unique.

<details>
Snakemake uses YAML (Yet Another Markup Language) for configuration. It is
human readable, and very often, you can directly edit the file with ease unlike
JSON. Further, YAML is a superset of JSON for the most parts so all YAML
libraries can work with JSON, if needed.
</details>


This Snakemake pipeline requires a YAML config file with the following
information (toggle to see).

<details>

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

</details>

### Configuration

In many cases the only information that will be changing between ATAC-seq
experiments is the **(a) library information, and (b) desired output
directory**; while paths to BWA indices, blacklists, etc. will remain
unchanged.

It may therefore by convenient to have a single permanent file with all of the
required information except the library information and the results dir.

If this is the case, you can use the python script at
`bin/make_atacseq_config.py` to add library information and the results path to
this unchanging YAML:

```bash
python /path/to/atacseq_snakemake/bin/make_atacseq_config.py \
      -r /path/to/results_dir                                \
      -ref atacseq.generic_data.yaml                         \
      -lib libraries.yaml > complete_atacseq_config.yaml
```

#### Generating `libraries.yaml`

Example code is provided in `bin/make_library_config.py` which would generate a
`libraries.yaml` file provided FASTQ files for your ATAC-seq experiment.

#### Generating `reference.yaml`

Example file containing reference data is given in
`config/atacseq.generic_data.yaml`.

### Running Snakemake

Snakemake would then require following files:

1. Snakefile (with your workflow)
2. Configuration (with your library information generated in previous step)
3. (Optionally) Cluster configuration if running on Cluster

A default cluster configuration for the pipeline is `config/cluster.yaml`.

The pipeline can then be run with simple command, which runs atacseq pipeline
with subsampling (if `config['subsample_depth']` is set).

```bash
$ make run_all
```

If subsampling is not desired, one can run the commands one at a time:

```bash
$ make run
$ make downsample
```

### Dry-run

```bash
$ make dry_run
```
