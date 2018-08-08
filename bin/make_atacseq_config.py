#!/usr/bin/env python

import argparse
import sys
import os
import yaml

__description__ = """
    Given
        i) YAML file with library info,
        ii) a YAML file containing other config information, and
        iii) a path to the desired results directory
    create combined config that can be used for the Snakemake ATAC-seq pipeline.

    Prints to STDOUT.
"""


def parse_arguments():
    parser = argparse.ArgumentParser(prog='python make_atacseq_config.py',
                                     description = __description__)

    parser.add_argument('-ref', action='store', required=True,
        help='Path to the yaml file encoding e.g. BWA index path (see README).')
    parser.add_argument('-lib', action='store', nargs = '+', required=True,
        help='Path to the yaml file encoding library information.')
    parser.add_argument('-r', '--results', action='store',
        help='Base directory in which the ATAC-seq analysis results should be stored (default: current working directory).')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    atacseq_config = {}
    atacseq_config['libraries'] = {}

    with open(args.ref, 'r') as r:
        generic_config = yaml.load(r)
        atacseq_config.update(generic_config)

    for library_yaml in args.lib:
        with open(library_yaml, 'r') as s:
            atacseq_config['libraries'].update(yaml.load(s))

    atacseq_config['results'] = args.results or os.getcwd()

    print(yaml.dump(atacseq_config, indent=4, default_flow_style=False))
