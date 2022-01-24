#!/usr/bin/env python3

"""
summarize lineages, blast and diamond output
"""

import os
import sys
import csv
import logging
import textwrap
import argparse
from Bio import SeqIO

import pandas as pd

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    """

    :return:
    """
    description = textwrap.dedent('''\
                                Summarize lineage assignment results
                                ------------------------------------
                                ''')
    epilog = textwrap.dedent('''\
                        Example usage: python report.py \
                        --query path/to/query/file(s) (fasta) \
                        --report /path/to/assignment/files (csv) \
                        --lineage-csv /path/to/known/lineages (csv) \
                        ''')

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=description,
                                     epilog=epilog
                                     )
    required = parser.add_argument_group('''Mandatory arguments''')

    required.add_argument('--query', action="store", type=str, metavar="<FILE(s)>", nargs='+', dest="query",
                          help="one or more FASTA files, paths should be separated by space")
    required.add_argument('--assignment', action="store", type=str, metavar="<FILE(s)>", nargs='+', dest='assignment',
                          help='one or more csv files from assign lineages function, paths should be separated by space'
                          )
    parser.add_argument('--blast', action="store", type=str, metavar="<FILE(s)>", nargs='+', dest='blast',
                        help='one or more tab separated DIAMOND BLASTX output files, paths should be separated by '
                             'space '
                        )
    required.add_argument("--lineage-csv", action="store", type=str, required=True, metavar="<FILE>",
                          dest="lineage",
                          help="comma-separated values (csv) file having the columns ['accession']")
    return parser


def aggregate_known_lineages(lineage):
    """

    :param lineage:
    :return:
    """
    # read the csv file
    lineages_df = pd.read_csv(lineage)

    # filter out rows with no accessions
    lineages_df = lineages_df[lineages_df['accession'].notnull()]
    # drop duplicates
    lineages_df = lineages_df.drop_duplicates(subset=['accession'], keep='first', inplace=False)

    # group by lineages
    lineage_info_df = lineages_df.groupby('lineage').agg(
        {'year': [min, max], 'lineage': 'count', 'country': lambda x: list(set(list(x)))}).reset_index()
    lineage_info_df.columns = lineage_info_df.columns.droplevel(0)

    lineage_info_df = lineage_info_df.rename(columns={"": "lineage",
                                                      "min": "year_first",
                                                      "max": "year_last",
                                                      "count": "n_seqs",
                                                      "<lambda>": "countries"}
                                                      )
    lineage_info_df['countries'] = lineage_info_df.countries.apply(';'.join)
    # create a dict with accession as key and lineage as value
    lineages_dict = lineage_info_df.set_index('lineage').T.to_dict(orient='list')
    return lineages_dict


def read_query(query):
    """
    read in query files 
    """

    # read in query files
    query_dict = {}
    for fa in query:
        record = SeqIO.read(handle=fa, format="fasta")
        seq = record.seq.lower()
        pct_ns = float(seq.count('n') / len(seq)) * 100
        pct_ns = "{:.2f}".format(pct_ns)
        query_dict[record.id] = [len(seq), pct_ns]
    return query_dict


def read_blast(blast):
    """
    read in blast files
    """

    # dict to get segments
    prot = {'YP_003848704.1': 'L', 'YP_003848705.1': 'M', 'YP_003848706.1': 'S', 'YP_003848707.1': 'S'}

    # read in blast output results
    blast_dict = dict()
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
                if prot.get(sseqid) is None:
                    continue
                else:
                    if not qseqid in blast_dict:
                        blast_dict[qseqid] = [length, sseqid, prot.get(sseqid), product, pident, mismatch, gaps]
                    else:
                        blast_dict[qseqid] += [length, sseqid, prot.get(sseqid), product, pident, mismatch, gaps]
    return blast_dict


def read_assignment(assignment):
    """

    :param assignment:
    :return:
    """

    # read in assignment outputs
    assignment_dict = {}
    for f in assignment:
        header = os.path.splitext(os.path.basename(f))[0].strip('.lineage.csv')
        assignment_dict[header] = f
    return assignment_dict


def summarize_results(query, assignment, lineage, outfile="lineages.csv"):
    """
    
    :param query: 
    :param assignment: 
    :param lineage: 
    :param outfile: 
    :return: 
    """

    query_dict = read_query(query=query)
    assignment_dict = read_assignment(assignment=assignment)
    lineages_dict = aggregate_known_lineages(lineage=lineage)

    columns = [
        'Query', 'Lineage', 'aLRT', 'UFbootstrap', 'Length',
        'Ns(%)', 'Note', 'Year_first', 'Year_last', 'Countries'
    ]

    csv_dict = dict()
    with open(outfile, 'w', newline='') as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(columns)
        for seqid, record in query_dict.items():
            if assignment_dict.get(seqid) is None:
                fasta_lst = record[:2]
                l = record[2:]
                for l in list(chunks(l, 7)):
                    lin = [seqid, 'unassigned', '0', '0'] + fasta_lst + l
                    writer.writerow(lin)
            else:
                f = assignment_dict.get(seqid)
                with open(f) as fin:
                    for line in fin:
                        lin = line.strip().split(',')
                        lin.extend(record)
                        if float(lin[3]) < 70:
                            lin.append("unassigned (bootstrap value < 70)")
                            lin.extend([None,None,None])
                        else:
                            lin.append("assigned (bootstrap value >= 70)")
                            if lineages_dict.get(lin[1]) is not None:
                                f = lineages_dict.get(lin[1])[0]
                                l = lineages_dict.get(lin[1])[1]
                                c = lineages_dict.get(lin[1])[3]
                                lin.extend([f, l, c])
                        writer.writerow(lin)
                        csv_dict[seqid] = lin


def summarize_blast(blast, outfile="diamond_results.csv"):
    """
    report blast results

    :param blast
    :param outfile
    :return:
    """
    blast_dict = read_blast(blast=blast)

    # write to output file
    columns = ['QueryID', 'Length', 'SubjectID', 'Segment', 'Product', 'PercentIdentity', 'Mismatches', 'Gaps']
    csv_dict = dict()
    with open(outfile, 'w', newline='') as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(columns)
        for seqid, record in blast_dict.items():
            record.insert(0, seqid)
            writer.writerow(record)
    return outfile


def chunks(l, n):
    n = max(1, n)
    return (l[i:i + n] for i in range(0, len(l), n))


def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()
    if args.query is None:
        print("required query file(s) missing!")
        sys.exit(1)
    if args.lineage is None:
        print("required lineages file missing!")
        sys.exit(1)
    if args.assignment is None:
        print("required lineage assignment file(s) missing!")
        sys.exit(1)

    if args.blast is not None:
        summarize_results(query=args.query, assignment=args.assignment, lineage=args.lineage)
        summarize_blast(blast=args.blast)
    else:
        summarize_results(query=args.query, assignment=args.assignment, lineage=args.lineage)


if __name__ == '__main__':
    main()
