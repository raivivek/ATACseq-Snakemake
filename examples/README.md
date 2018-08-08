This directory contains examples of the files needed to employ the Snakemake
pipeline.

Given,

1. a YAML file with library information (see `run_2156.yaml` for an example),
2. a YAML file with paths for the necessary 'generic files' (e.g. BWA indices,
   BED files containing TSS locations, and whitelists/blacklists against which
   peaks should be filtered; see `atacseq.generic_data.yaml`),
3. a path to output directory,

one can create the config file needed for the Snakemake pipeline as follows:

```bash
python /path/to/atacseq_snakemake/bin/make_atacseq_config.py \
      -r /path/to/results_dir                                \
      -ref atacseq.generic_data.yaml                         \
      -lib libraries.yaml > complete_atacseq_config.yaml
```

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
        --snakefile /path/to/atacseq_snakemake/src/Snakefile \
        --cluster-config /path/to/cluster.config \
        --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" \
        &> snakemake.nohup &
```

### Dry-run

It recommended to perform a dry-run of Snakemake to see if it properly generates
all the files in workflow before submitting to cluster. For a dry-run, add `-n`
flag to the command above.
