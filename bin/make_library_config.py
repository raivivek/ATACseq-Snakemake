#! /usr/bin/env python

import sys
import re
import yaml
import pathlib

# match files such as atacseq.1.fq.gz or atacseq.2.fastq.gz
FASTQ_RE = '(.*).([12]).(?:fastq|fq).*'

LIBRARIES = {}


def parse_fastq_name(f):
    """ Return FASTQ "library" name and whether it is "first" or "second" set
    of reads (paired-end)."""

    m = re.search(FASTQ_RE, os.path.basename(f))
    library = m.group(1)
    first_or_second = m.group(2)
    return [library, first_or_second]

def create_library_item(fastq):
    """ Parse FASTQ filenames; and return a dictionary with library names
    containing FASTQ files. For example, a library entry looks like
    
    {
        'genome': 'hg19',
        'readgroups': ['ABCD123.1.fastq.gz', 'ABCD123.2.fastq.gz']
    }
    
    """

    pass




if __name__=='__main__':
    #print(yaml.dump(LIBRARIES, indent=4, default_flow_style=False))
    pass

