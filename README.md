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

Also, assumes that picard MarkDuplicates can be called using the syntax: `picard MarkDuplicates ...`.

## Usage:

This Snakemake pipeline requires a config file (JSON format) with the following information:

```bash
{
    "blacklist": {  # (Optional) For each genome, a list of blacklisted regions in bed format
                    # (Not required by the pipeline, but should be used if they are available!!).
                    # These are used for peak filtering, and by ataqv.
        "hg19": [
            "/lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz", 
            "/lab/data/reference/human/hg19/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
        ]
    }, 
    "bwa_index": { # (Required) path to BWA indices for each genome needed
        "hg19": "/lab/data/reference/human/hg19/index/bwa/current/hg19",
        "mm9": "/lab/data/reference/mouse/mm9/index/bwa/current/mm9"
    }, 
    "results": "/lab/work/porchard/atacseq", # (Optional) Path to the directory in which results should be placed (default is current working directory is used)
    "libraries": { # (Required) this is where the information for each library is given
        "100474___2156": { # unique ID for first library
            "genome": "hg19", # genome for first library
            "readgroups": { # readgroups for first library
                            # if the library was sequenced across several lanes, multiple readgroups
                            # can be provided, and they will be merged after mapping and before duplicate marking/filtering
                            # in this case, the library was sequenced across four lanes so four readgroups are provided.
                "100474___L1___2156": [ # list of the 2 fastq files for the first lane
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L001.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L001.2.fastq.gz"
                ], 
                "100474___L2___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L002.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L002.2.fastq.gz"
                ], 
                "100474___L3___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L003.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L003.2.fastq.gz"
                ], 
                "100474___L4___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L004.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100474_L004.2.fastq.gz"
                ]
            }
        }, 
        "100477___2156": { # second library begins here
            "genome": "hg19", 
            "readgroups": {
                "100477___L1___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L001.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L001.2.fastq.gz"
                ], 
                "100477___L2___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L002.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L002.2.fastq.gz"
                ], 
                "100477___L3___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L003.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L003.2.fastq.gz"
                ], 
                "100477___L4___2156": [
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L004.1.fastq.gz", 
                    "/lab/work/porchard/snakemake_atacseq/data/fastq/100477_L004.2.fastq.gz"
                ]
            }
        }, 
    }, 
    "tss": { # (Required) for each genome, path to a file of TSSs in bed format (used by ataqv for ATAC-seq quality control)
                # example TSS files can be found in the ataqv GitHub repo
        "hg19": "/lab/data/reference/human/hg19/annot/hg19.tss.refseq.bed", 
        "rn5": "/lab/data/reference/rat/rn5/annot/rn5.tss.refseq.bed"
    }, 
    "whitelist": { # (Optional) if whitelists are present rather than blacklists, they can be provided as well.
                  # These are only used in filtering the peak lists
        "rn5": "/lab/data/reference/rat/rn5/annot/rn5.K30.mappable_only.bed.gz"
    }
}

```

**IMPORTANT**: the basename for each fastq file must be unique, as must be the
readgroup names (keys). In many cases the only information that will be changing
between ATAC-seq experiments is the library information and the desired output
directory (paths to BWA indices, blacklists, etc. will remain unchanged). It may
therefore by convenient to have a single permanent JSON file with all of the
required information except the library information and the results dir. If this
is the case, you can use the python script at `src/make_atacseq_config.py` to
add library information and the results path to this unchanging JSON:

```bash
python bin/make_atacseq_config.py -r /path/to/results_dir /path/to/json_with_everything_except_libraries_and_results.json /path/to/json_with_libraries.json
```

Example files are given given in `examples/`.

In case you are running on a cluster and need a cluster config file for
Snakemake, a template cluster config can be found in `src/` as well.
