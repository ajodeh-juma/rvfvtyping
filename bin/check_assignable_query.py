#!/usr/bin/env python

"""
get assignable input query FASTA sequences
"""

import os
import sys
import argparse


def parse_args():
    """

    :return:
    """
    parser = argparse.ArgumentParser(prog="check_assignable_query.py", description=__doc__)
    required = parser.add_argument_group('''input options''')
    required.add_argument('--blast', action="store", type=str, metavar="<FILE(s)>", nargs='+', dest='blast',
                          help='one or more tab separated DIAMOND BLASTX output files, paths should be separated by '
                               'space '
                          )
    parser.add_argument('--out', metavar="<FILE>", type=str, dest="out", help="output file to write summary")
    return parser


def check_queries(blast, out):
    """
    
    :param blast
    :param outfile
    :return:
    """

    # dict to get segments
    prot = {'YP_003848704.1': 'L', 'YP_003848705.1': 'M', 'YP_003848706.1': 'S', 'YP_003848707.1': 'S'}

    blast_dict = dict()
    to_drop = []
    to_keep = {}
    for fn in blast:
        with open(fn) as f_obj:
            for line in f_obj:
                line = line.strip()
                qseqid = line.split('\t')[0]
                sseqid = line.split('\t')[1]
                pident = line.split('\t')[2]
                length = int(line.split('\t')[3]) * 3
                mismatch = line.split('\t')[4]
                gaps = line.split('\t')[5]
                evalue = float(line.split('\t')[10])
                bitscore = float(line.split('\t')[11])
                sscinames = line.split('\t')[13]
                stitle = line.split('\t')[14]
                product = stitle.rsplit("[")
                product = ' '.join(product[0].strip().rsplit(" ")[1:])
                aln_score = (int(length) - (int(mismatch) + int(gaps))) * 3

                if not qseqid in to_keep:
                    to_keep[qseqid] = [sseqid]
                else:
                    to_keep[qseqid].append(sseqid)

                if prot.get(sseqid) is None:
                    continue
                else:
                    if not qseqid in blast_dict:
                        blast_dict[qseqid] = [length, sseqid, prot.get(sseqid), product, pident, mismatch, gaps]
                    else:
                        blast_dict[qseqid] += [length, sseqid, prot.get(sseqid), product, pident, mismatch, gaps]

    for k, v in to_keep.items():
        if 'YP_003848705.1' not in v:
            to_drop.append(k)

    for k in to_drop:
        to_keep.pop(k)

    with open(out, 'w') as f_obj:
        if len(to_keep) == 0:
            f_obj.write("to_keep {} true".format(len(to_keep)))
        else:
            f_obj.write("to_drop {} false".format(len(to_drop)))
    # return len(to_keep)
    return out


def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()
    if args.blast is None:
        sys.exit(1)

    m = check_queries(blast=args.blast, out=args.out)
    # print(m, end='')


if __name__ == '__main__':
    main()
