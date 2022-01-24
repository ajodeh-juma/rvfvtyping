#!/usr/bin/env python

"""
summarize lineages csv files
"""

import os
import sys
import logging
import argparse
from utils import fasta_iterator

import pandas as pd

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    """

    :return:
    """
    parser = argparse.ArgumentParser(prog="snpsTodataframe.py", description=__doc__)
    required = parser.add_argument_group('''input options''')

    required.add_argument('--fasta', action="store", type=str, metavar="<FILE>", dest="fasta",
                          help="path to the masked snps FASTA file (all sequence lengths are equal)")
    required.add_argument('--outfile', action="store", type=str, metavar="<FILE>", dest='outfile',
                          help='filename of the output file'
                          )
    return parser


def fasta2df(fasta, outfile):
    """
    convert fasta sequnec to a dataframe where column names are the sequence identifiers and rows are the nucleotide bases

    :param fasta
    :param outfile
    :return:
    """

    # convert fasta to dict
    fasta_dic = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fasta))
    
    # convert dict to dataframe
    df = pd.DataFrame.from_dict(fasta_dic, orient='index', columns=['sequence']).reset_index().rename(columns={'index': 'id'})

    csv_dict = dict()
    # iterate through the dataframe and create a dict with sequence id as key and a list of nucleotides as values
    for index, row in df.iterrows():
        csv_dict[row['id']] = [base.upper() for base in row['sequence']]
    df2 = pd.DataFrame.from_dict(csv_dict)
    # write to output file
    outfile = os.path.abspath(outfile)
    if not os.path.isdir(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    df2.to_csv(outfile, index=False, sep=",", index_label=True)
    return outfile


def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()
    if args.fasta is None:
        print("required FASTA file missing!")
        sys.exit(1)
    if args.outfile is None:
        print("required output file missing!")
        sys.exit(1)

    fasta2df(fasta=args.fasta, outfile=args.outfile)


if __name__ == '__main__':
    main()