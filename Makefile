# Makefile
#
# {{cookiecutter.name}}
# {{cookiecutter.email}}
# (c) {{cookiecutter.affiliation}}
#
# {{cookiecutter.date}}
#

.PHONY = dry_run run

# --jn: job name
# -n: dry-run only
# -r: output reason
# -p: commands run
dry_run:
	@snakemake -npr \
		--jn "snakejob.{jobid}" \
		--snakefile src/Snakefile \
		--configfile config/config.yaml

# nohup: run in background 
# -j: maximum number of jobs to put in queue
#	--keep-going: keep going with independent jobs if some fail
# --rerun-incomplete: re-run any incomplete rules
run:
	@nohup snakemake \
		--jn "snakejob.{jobid}" \
		-j 999 \
		--keep-going \
		--rerun-incomplete \
		--snakefile src/Snakefile \
		--configfile config/config.yaml \
		--cluster-config config/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus}" \
		> logs/snakemake.log
