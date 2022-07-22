#!/usr/bin/env python


import re
import argparse
import logging
import textwrap
from utils import fasta_iterator


def parse_args():
    parser = argparse.ArgumentParser(prog="filter_input.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=textwrap.dedent('''\
                                            filter consensus/query input fasta
                                            ------------------------------------------------------------------

                                            fasta
                                            '''),
                                     argument_default=argparse.SUPPRESS)
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--fasta', metavar='<FILE>', type=str,
                                dest="fasta", help="fasta file"
                                )
    required_group.add_argument('--out', metavar="<FILE>", type=str, dest="out", help="output file to write summary")
    return parser


def filter_input(fasta, out):
    """

    :param fasta:
    
    :return:
    """

    to_exclude = []

    d = dict((re.sub(r'\W+', "_", x[0]), x[1]) for x in fasta_iterator(fasta_name=fasta))
    # d = dict((x[0], x[1]) for x in fasta_iterator(fasta_name=fasta))
    with open(out, 'w') as f_out:
        for seqid, seq in d.items():
            seq = seq.lower()
            pct_ns = float(seq.count('n') / len(seq)) * 100
            pct_ns = "{:.2f}".format(pct_ns)
            f_out.write("seqid {}\n".format(seqid))
            f_out.write("length {}\n".format(len(seq)))
            f_out.write("perc_Ns {}\n".format(pct_ns))
            f_out.write("seq {}\n".format(seq))
    return out


def main():
    """
    run the functions
    :return:
    """
    parser = parse_args()
    args = parser.parse_args()
    n = filter_input(fasta=args.fasta, out=args.out)


if __name__ == '__main__':
    main()
