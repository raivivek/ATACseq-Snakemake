# Makefile
#
# {{cookiecutter.name}}
# {{cookiecutter.email}}
# (c) {{cookiecutter.affiliation}}
#
# {{cookiecutter.date}}
#

.PHONY = dry_run run

dry_run:
	@snakemake -npr \ # dry-run only (-n); output reason (-r) and commands run (-p)
		--jn "snakejob.{jobid}" \
		--snakefile src/Snakefile \
		--configfile config/config.yaml

run:
	@nohup snakemake \ # run in background
		--jn "snakejob.{jobid}" \
		-j 999 \ # maximum number of jobs to put in queue
		--keep-going \ # keep going with independent jobs if some fail
		--rerun-incomplete \ # re-run any incomplete rules
		--snakefile src/Snakefile \
		--configfile config/config.yaml \
		--cluster-config config/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus}" \ # job-scheduler
		> logs/snakemake.log # pipe log to a file
