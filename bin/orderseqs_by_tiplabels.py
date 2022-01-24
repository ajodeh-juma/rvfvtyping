#!/usr/bin/env python

"""
order sequences according to the tip labels in a tree
"""

import os
import sys
import logging
import argparse
import textwrap

from dendropy import Tree

from utils import fasta_iterator

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    """

    :return:
    """
    parser = argparse.ArgumentParser(prog="orderseqs_by_tiplabels.py", description=__doc__)
    required = parser.add_argument_group('''input options''')
    
    required.add_argument('--tree', action="store", type=str, metavar="<FILE>", dest='tree',
                          help='path to the tree file in nexus format'
                          )
    required.add_argument('--alignment', action="store", type=str, metavar="<FILE>", dest="alignment",
                          help="path to the alignment FASTA file")
    parser.add_argument('--schema', dest='schema', action='store', default="nexus",
                        choices=['nexus', 'newick', 'nexml'],
                        help='what format is the tree file. This is passed to dendropy. default is \'nexus\'')
    required.add_argument('--outfile', action="store", type=str, metavar="<FILE>", dest='outfile',
                          help='filename of the output ordered FASTA file'
                          )
    return parser


def rearrange_seqs(tree, alignment, outfile, schema):
    """
    order sequences as per the tip labels in a tree (nexus format)

    :param
    :param

    """

    # convert fasta alignment file to dict
    fasta_dic = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=alignment))

    # read the tree
    tree = Tree.get(path=tree, schema=schema.lower(), preserve_underscores=True)

    with open(outfile, 'w') as fh_obj:
        for tip in tree.leaf_node_iter():
            fh_obj.write((">{}\n{}\n".format(tip.taxon.label, textwrap.fill(fasta_dic.get(tip.taxon.label), width=80))))
    return outfile


def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()
    if args.alignment is None:
        print("required alignment FASTA file missing!")
        sys.exit(1)
    if args.tree is None:
        print("required tree file in nexus format missing!")
        sys.exit(1)
    if args.outfile is None:
        print("required output file (.fasta) is missing!")
        sys.exit(1)

    rearrange_seqs(tree=args.tree, alignment=args.alignment, outfile=args.outfile, schema=args.schema)


if __name__ == '__main__':
    main()

