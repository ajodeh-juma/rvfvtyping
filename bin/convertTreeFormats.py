#!/usr/bin/env python

"""
convert tree file from newick format to nexus format
"""

import os
import sys
import logging
import argparse
from Bio import Phylo

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    """

    :return:
    """
    parser = argparse.ArgumentParser(prog="convertTreeFormats.py", description=__doc__)
    parser.add_argument('--in-tree', action="store", type=str, required=True, help="input tree filename",
                        metavar="<FILE>", dest="in_tree")
    parser.add_argument("--in-format", action="store", type=str, required=True,
                        help="format of the input tree file, available choices ['newick', 'nexus']",
                        metavar="str", dest="in_format", choices=['newick', 'nexus'])
    parser.add_argument('--out-tree', action="store", type=str, required=True, help="output tree filename",
                        metavar="<FILE>", dest="out_tree")
    parser.add_argument("--out-format", action="store", type=str, required=True,
                        help="format of the output tree file, available choices ['newick', 'nexus']",
                        metavar="str", dest="out_format", choices=['newick', 'nexus'])
    return parser


def phylo_convert(in_tree, out_tree, in_format, out_format):
    """

    :return:
    """

    return Phylo.convert(in_tree, in_format, out_tree, out_format)


def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()

    if not os.path.isfile(args.in_tree):
        logging.error('cannot find tree file at {}'.format(args.in_tree))
        sys.exit(1)
    if args.in_format == args.out_format:
        logging.info('tree formats are similar: {} & {}'.format(args.in_format, args.out_format))
        sys.exit(1)

    phylo_convert(
        in_tree=args.in_tree, in_format=args.in_format,
        out_tree=args.out_tree, out_format=args.out_format
    )


if __name__ == '__main__':
    main()
