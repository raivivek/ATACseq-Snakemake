This directory contains examples of the files needed to employ the Snakemake pipeline.

If one has a JSON file with library information (see `run_2156.json` for an example) and 
a JSON file with paths for the necessary 'generic files' (e.g. BWA indices, BED files 
containing TSS locations, and whitelists/blacklists against which peaks should be 
filtered; see `atacseq.generic_data.json`) one can create the config file needed 
for the Snakemake pipeline as follows: 

```bash
python /path/to/atacseq_snakemake/bin/make_atacseq_config.py -r /path/to/results_dir atacseq.generic_data.json libraries.json > complete_atacseq_config.json
```

If running on a cluster, a second config file containing resource allocation
information will be useful as well. A fitting example/template is given in
`cluster.config`. The pipeline can then be run:

```bash
nohup snakemake -p -j 5 \
                --configfile /path/to/complete_atacseq_config.json \
                --snakefile /path/to/atacseq_snakemake/src/Snakefile \
                --cluster-config /path/to/cluster.config
                --cluster "sbatch -t {cluster.time} -n {cluster.N} --mem-per-cpu={cluster.mem}" &> snakemake.nohup &
```
