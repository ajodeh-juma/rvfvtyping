#!/usr/bin/env python


import argparse
import logging
import textwrap
from utils import fasta_iterator

def parse_args():
    parser = argparse.ArgumentParser(prog="check_segments.py",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description=textwrap.dedent('''\
                                            check segments from the classification for sequences
                                            ------------------------------------------------------------------

                                            blast-txt
                                            '''),
                                        argument_default=argparse.SUPPRESS)
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--blast', metavar='<FILE>', type=str,
                                dest="blast", help="DIAMOND BLAST output file (tab-delimited format - 6)"
                                )
    required_group.add_argument('--out', metavar="<FILE>", type=str, dest="out", help="output file to write summary")
    return parser

def check_blast(blast, out):
    """

    :param blast:
    
    :return:
    """

    with open(blast) as fn:
        line = fn.readline() # only read the first line since DIAMOND sorts based on score and evalue
        with open(out, 'w') as f_out:
            line  = line.strip()
            qseqid = line.split('\t')[0]
            sseqid = line.split('\t')[1]
            length = line.split('\t')[3]
            mismatch = line.split('\t')[4]
            gaps = line.split('\t')[5]
            aln_score = (int(length) - (int(mismatch) + int(gaps))) * 3

            f_out.write("qseqid {}\n".format(qseqid))
            f_out.write("sseqid {}\n".format(sseqid))
            f_out.write("length {}\n".format(length))
            f_out.write("mismatches {}\n".format(mismatch))
            f_out.write("gaps {}\n".format(gaps))
            f_out.write("alignment_score {}\n".format(aln_score))
    return out

def main():
    """
    run the functions
    :return:
    """
    parser = parse_args()
    args = parser.parse_args()
    check_blast(blast=args.blast, out=args.out)

if __name__ == '__main__':
    main()
