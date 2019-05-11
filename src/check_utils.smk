#
# Checks versions for different utilites
#
def check_utils_version(outfile):
    """Check if the tools exists, and if they do, collect their version number.
    Otherwise, exit with an unpleasant error message."""

    utils = [
      "fastqc -v",
      "cta --version",
      "cutadapt --version",
      "samtools --version",
      "bwa",
      "macs2 --version",
      "ataqv --version",
      "picard --version",
      "bedtools --version",
      "pigz --version"
    ]

    all_good = True


    with open(outfile, "w") as f:
      for util in utils:
        binary, args = util.split()
        if not shutil.which(binary):
          all_good = False
          msg = f"{binary} not found!\n"
        else:
          msg = subprocess.check_output(util, shell=True, encoding="utf-8")

        f.write(msg)

    return all_good


##
## debug info
##

rule versions:
    output:
        _versions("version.txt")
    run:
        all_good = check_utils_version(output)
        if not all_good:
            print("Check $PATH. Aborting.", file=sys.stderr)
            sys.exit(1)

# vim:syntax=python
