# Makefile
#
# {{name}}
# {{email}}
# (c) {{affiliation}}
#
# {{date}}
#

.PHONY = dry_run run

# --jn: job name
# -n: dry-run only
# -r: output reason
# -p: commands run
dry_run:
	@snakemake -npr make_libraries make_samples\
		--jn "atacseq.{jobid}" \
		--snakefile src/Snakefile \
		--configfile tmp/config.yaml

unlock:
	@snakemake --unlock \
		--snakefile src/Snakefile \
		--configfile tmp/config.yaml

# nohup: run in background 
# -j: maximum number of jobs to put in queue
#	--keep-going: keep going with independent jobs if some fail
# --rerun-incomplete: re-run any incomplete rules
#

run_all: run downsample

run:
	@nohup snakemake make_libraries make_samples \
		--jn "atacseq.{jobid}" \
		-j 999 \
		--keep-going \
		--rerun-incomplete \
		--snakefile src/Snakefile \
		--configfile config/config.yaml \
		--cluster-config config/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --mail-type=FAIL --mail-user=vivekrai@wolverine.theparkerlab.org" \
		> logs/snakemake.log&

downsample:
	@nohup snakemake downsample\
		--jn "atacseq.{jobid}" \
		-j 999 \
		--keep-going \
		--rerun-incomplete \
		--snakefile src/Snakefile \
		--configfile config/config.yaml \
		--cluster-config config/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus}" \
		> logs/snakemake.log&
