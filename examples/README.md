This directory contains examples of the files needed to employ the Snakemake
pipeline.

If one has a YAML file with library information (see `run_2156.yaml` for an
example) and a YAML file with paths for the necessary 'generic files' (e.g. BWA
indices, BED files containing TSS locations, and whitelists/blacklists against
which peaks should be filtered; see `atacseq.generic_data.yaml`) one can create
the config file needed for the Snakemake pipeline as follows: 

```bash
python /path/to/atacseq_snakemake/bin/make_atacseq_config.py \
  -r /path/to/results_dir \
  -ref atacseq.generic_data.yaml \
  -lib libraries.yaml > complete_atacseq_config.yaml
```

If running on a cluster, a second config file containing resource allocation
information will be useful as well. A fitting example/template is given in
`cluster.config`. The pipeline can then be run:

```bash
nohup snakemake -p -j 10 \
                --configfile /path/to/complete_atacseq_config.yaml \
                --snakefile /path/to/atacseq_snakemake/src/Snakefile \
                --cluster-config /path/to/cluster.config \
                --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" \
                &> snakemake.nohup &
```
