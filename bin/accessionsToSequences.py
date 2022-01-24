#!/usr/bin/env python3


"""

"""

import os
import re
import sys
import logging
import argparse
import textwrap
import pycountry

import pandas as pd
from Bio import SeqIO
from Bio import Entrez
from itertools import takewhile

from utils import start_time
from utils import end_time
from utils import run_time
from utils import str2bool
from utils import fasta_iterator

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    """
    arguments parser for command line options

    :return:
    """

    parser = argparse.ArgumentParser(prog="accessionsToSequences.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=textwrap.dedent('''\
                                         fetch sequences from NCBI given sequence accession numbers
                                         ------------------------------------------------------------------

                                         lineages.csv input file should have a column labelled 'accession'
                                         isolate,year,country,source,lineage,accession
                                         ZH-548,1977,Egypt,Human,A,NC_014396
                                         '''),
                                     argument_default=argparse.SUPPRESS)
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--accessions', type=str, dest="accessions", metavar="<FILE>",
                                help="comma-separated values (csv) file having a column with accessions labelled "
                                     "as 'accession'"
                                )
    required_group.add_argument('--fasta', metavar='<FILE>', type=str,
                                dest="fasta", help="text file to write the search output"
                                )
    required_group.add_argument('--database', type=str, metavar='<str>', default="nucleotide",
                                choices=["nucleotide", "gene", "protein"], help="NCBI database name"),
    parser.add_argument('--out_dir', metavar='<dir>', dest='out_dir', default="./",
                        help="path to the directory to write the output files"
                        )

    parser.add_argument('--batch-size', type=int, required=False, metavar='batch_size', default=100,
                        help='number of accessions to process per request')
    parser.add_argument('--email', type=str, metavar="<email>", dest="email", required=False,
                        default="someone@gmail.com", help="email address to access Entrez"
                        )
    parser.add_argument('--cleanup', metavar="<boolean>", type=str2bool, default=False,
                        dest='cleanup', choices=[True, False],
                        help='delete the genbank (.gbk) output file once downloaded'
                        )
    return parser


def read_accessions(accessions):
    """

    :param accessions:
    :return:
    """

    # read the csv file
    lineages_df = pd.read_csv(accessions)

    # filter out rows with no accessions
    lineages_df = lineages_df[lineages_df['accession'].notnull()]

    # drop duplicates
    lineages_df = lineages_df.drop_duplicates(subset=['accession'], keep='first', inplace=False)
    # dict to store accession as key, lineage as value
    lineages_dict = lineages_df.set_index("accession")["lineage"].to_dict()
    accessions = list(lineages_df.accession)
    return accessions, lineages_dict


def accessions_to_gb(accessions, database, batch_size, ret_max=10 ** 9):
    """

    :param accessions:
    :param database:
    :param batch_size:
    :param ret_max:
    :return:
    """

    def batch(sequence, size):
        """

        :param sequence:
        :param size:
        :return:
        """
        length = len(accessions)
        for start in range(0, length, size):
            yield sequence[start:min(start + size, length)]

    def extract_records(records_handle):
        """

        :param records_handle:
        :return:
        """
        logging.info("getting records for accession")
        buffer = []
        for line in records_handle:
            if line.startswith("LOCUS") and buffer:
                # yield accession number and record
                yield buffer[0].split()[1], "".join(buffer)
                buffer = [line]
            else:
                buffer.append(line)
        yield buffer[0].split()[1], "".join(buffer)

    def process_batch(accessions_batch):
        """

        :param accessions_batch:
        :return:
        """
        # get GI for query accessions
        logging.info("getting GI accessions")
        query = " ".join(accessions_batch)
        query_handle = Entrez.esearch(db=database, term=query, retmax=ret_max)
        gi_list = Entrez.read(query_handle)['IdList']

        # get GB files
        logging.info("getting the genbank file format")
        search_handle = Entrez.epost(db=database, id=",".join(gi_list))
        search_results = Entrez.read(search_handle)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
        records_handle = Entrez.efetch(db=database, rettype="gb", retmax=batch_size,
                                       webenv=webenv, query_key=query_key)
        yield from extract_records(records_handle)

    accession_batches = batch(accessions, batch_size)
    for acc_batch in accession_batches:
        yield from process_batch(acc_batch)


def write_record(out_dir, accession, record, cleanup):
    """

    :param out_dir:
    :param accession:
    :param record:
    :return:
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    gbk = os.path.join(out_dir, accession + '.gbk')
    fa = os.path.join(out_dir, accession + '.fasta')
    metadata = os.path.join(out_dir, accession + '.tsv')

    with open(gbk, "w") as f_obj:
        print(record, file=f_obj)

    # parse the genbank file
    logging.info("[parsing GenBank file {:>5}]".format(gbk))
    records = SeqIO.parse(gbk, "genbank")  # parse records from the genbank file

    gbk_dict = dict()
    # use module pycountry to create a dictionary of counties and their 3-letter codes/symbols
    countries = dict()
    t = list(pycountry.countries)
    for country in t:
        countries[country.name.split(',')[0]] = country.alpha_3

    with open(fa, 'w') as fh:
        for record in records:
            seq_annotations = record.annotations
            accession = record.id.rsplit('.')[0]
            if len(record.seq) < 490:
                os.remove(fa)
                continue
            else:
                gbk_dict[record.id.rsplit('.')[0]] = [seq_annotations['organism']]
                source = record.features[0]
                for qualifiers in source.qualifiers:
                    if qualifiers == 'strain':
                        strain = ''.join(source.qualifiers['strain'])
                        if accession in gbk_dict:
                            gbk_dict[accession].append(strain)
                    elif qualifiers == 'isolate':
                        isolate = ''.join(source.qualifiers['isolate'])
                        if accession in gbk_dict:
                            gbk_dict[accession].append(isolate)
                    elif qualifiers == 'clone':
                        clone = ''.join(source.qualifiers['clone'])
                        if accession in gbk_dict:
                            gbk_dict[accession].append(clone)
                    elif qualifiers == 'country':
                        country = ''.join(source.qualifiers['country']).rsplit(':')
                        if accession in gbk_dict:
                            gbk_dict[accession].append(country[0])

                    elif qualifiers == 'collection_date':
                        collection_date = ''.join(source.qualifiers['collection_date'])
                        if accession in gbk_dict:
                            gbk_dict[accession].append(collection_date)

                gbk_dict[accession].extend([seq_annotations['references'][0].title,
                                            seq_annotations['references'][0].authors,
                                            seq_annotations['references'][0].journal])
                fh.write((">{}\n{}\n".format(accession, textwrap.fill(str(record.seq), width=80))))

    csv_dict = dict()
    with open(metadata, 'w') as csv_file:
        csv_file.write("Accession\tOrganism\tStrain\tCountry\tYear\tTitle\tAuthors\tJournal")
        csv_file.write("\n")
        for accession, rec in gbk_dict.items():
            if len(rec[1:]) == 0:
                continue
            else:
                rec = ['NA' if v == ' ' else v for v in rec]
                csv_dict[accession] = rec

        try:
            for k, v in csv_dict.items():
                csv_file.write(k + '\t' + "\t".join(v))
                csv_file.write("\n")
        except IOError as error:
            print(f"I/O error: {error}")

    if cleanup is True:
        os.remove(gbk)
    return fa, metadata


def write_to_fasta(fa_files, fasta, metadata_files):
    """

    :param fa_files:
    :param fasta:
    :return:
    """

    fasta = os.path.abspath(fasta)
    metadata = os.path.join(os.path.dirname(fasta), "metadata_references.tsv")
    with open(fasta, 'wt') as fh:
        for fn in fa_files:
            if os.path.isfile(fn):
                seq_dict = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fn))
                for seq_id, seq in seq_dict.items():
                    logging.info("writing {} sequences into {}".format(fn, fasta))
                    fh.write((">{}\n{}\n".format(seq_id, textwrap.fill(seq, width=80))))

    with open(metadata, 'wt') as f_out:
        f_out.write("Accession\tOrganism\tStrain\tCountry\tYear\tTitle\tAuthors\tJournal")
        f_out.write("\n")
        for fn in metadata_files:
            if os.path.isfile(fn):
                for line in open(fn):
                    if line.startswith('Accession'):
                        continue
                    else:
                        f_out.write("\t".join(line.strip().split('\t')))
                        f_out.write("\n")
    return fasta


def main():
    """
    run the functions
    :return:
    """
    parser = parse_args()
    args = parser.parse_args()
    st = start_time()
    input_file = args.accessions
    if not os.path.isfile(input_file):
        logging.error("File {} does not exist")
        sys.exit(1)
    else:
        accessions, lineages_dict = read_accessions(input_file)

    dbase = args.database
    Entrez.email = args.email
    batch_size = args.batch_size
    out_dir = args.out_dir
    fasta = args.fasta
    cleanup = args.cleanup

    fns = list()
    meta_fns = list()
    for acc, record in accessions_to_gb(accessions=accessions, database=dbase, batch_size=batch_size, ret_max=10 ** 9):
        fa, meta = write_record(out_dir=out_dir, accession=acc, record=record, cleanup=cleanup)
        fns.append(fa)
        meta_fns.append(meta)

    write_to_fasta(fa_files=fns, fasta=fasta, metadata_files=meta_fns)

    if cleanup:
        for fn in meta_fns:
            os.remove(fn)
    et = end_time(st)
    run_time(st, et)


if __name__ == '__main__':
    main()
