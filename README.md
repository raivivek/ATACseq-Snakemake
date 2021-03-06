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


## Configuration

Snakemake utilizes the `basename` of `fastq` files to organize and
generate files, as must be the readgroup names; thus they **must** be unique.

Snakemake uses YAML (Yet Another Markup Language) for configuration. It is human
readable, and very often, you can directly edit the file with ease unlike JSON.
Further, YAML is a superset of JSON for the most parts so all YAML libraries can
work with JSON, if needed.

In many cases the only information that will be changing between ATAC-seq
experiments is the **(a) library information, and (b) desired output
directory**; while meta-data information such as paths to BWA indices,
blacklists, etc. will remain unchanged.

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

For example, the combined configuration file might look like:

```yaml
results: /lab/work/porchard/atacseq
blacklist:
    hg19:
    - /lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
bwa_index:
    hg19: /lab/data/reference/human/hg19/index/bwa/current/hg19
tss:
    hg19: /lab/data/reference/human/hg19/annot/hg19.tss.refseq.bed
libraries:
    '100474___2156':
        genome: hg19
        readgroups:
            100474___2156_RG1:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2156_RG1.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2156_RG1.2.fastq.gz
            100474_RG2:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2156_RG2.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2156_RG2.2.fastq.gz
    '100477___2157':
        genome: hg19
        readgroups:
            100477___2157_RG1:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2157_RG1.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2157_RG1.2.fastq.gz
            100477___2157_RG2:
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2157_RG2.1.fastq.gz
            - /lab/work/porchard/snakemake_atacseq/data/fastq/100477___2157_RG2.2.fastq.gz
```

### Generating `libraries.yaml`

Example code is provided in `bin/make_library_config.py` which would generate a
`libraries.yaml` file provided FASTQ files for your ATAC-seq experiment.

### Generating `reference.yaml`

Example file containing reference data is given in
`config/atacseq.generic_data.yaml`.

## Usage

### Running Snakemake

Snakemake would then require following files:

1. Snakefile (with your workflow)
2. Configuration (with your library information generated in previous step)
3. (Optionally) Cluster configuration if running on Cluster

A default cluster configuration for the pipeline is `config/cluster.yaml`.

```bash
$ make run
```

### Dry-run

```bash
$ make dry_run
```
