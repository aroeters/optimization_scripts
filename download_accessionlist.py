#! /usr/bin/env python3

# Eli Korvigo
# https://www.biostars.org/p/66921/
# use biopython environment
# see available databases here: https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly


import argparse
import sys
import os

import Bio.Entrez


RETMAX = 10**9
GB_EXT = ".gbk"


def parse_args(arg_lst):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="A file with accessions to download")
    parser.add_argument("-d", "--database", type=str, required=True,
                        help="NCBI database ID")
    parser.add_argument("-e", "--email", type=str, required=False,
                        default="some_email@somedomain.com",
                        help="An e-mail address")
    parser.add_argument("-b", "--batch", type=int, required=False, default=200,
                        help="The number of accessions to process per request")
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="The directory to write downloaded files to")

    return parser.parse_args(arg_lst)


def read_accessions(fp):
	access_list = list()
	for line in open(fp):
		#~ access_list.append(line.split()[0].strip())
		access_list.append(line.strip())
	return access_list


def accessions_to_gb(accessions, db, batchsize, retmax):
    def batch(sequence, size):
        l = len(accessions)
        for start in range(0, l, size):
            yield sequence[start:min(start + size, l)]

    def extract_records(records_handle):
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
        #~ # get GI for query accessions
        #~ query = " ".join(accessions_batch)
        #~ query_handle = Bio.Entrez.esearch(db=db, term=query, retmax=retmax)
        #~ gi_list = Bio.Entrez.read(query_handle)['IdList']
        gi_list = accessions_batch
        # get GB files
        search_handle = Bio.Entrez.epost(db=db, id=",".join(gi_list))
        search_results = Bio.Entrez.read(search_handle)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
        records_handle = Bio.Entrez.efetch(db=db, rettype="gb", retmax=batchsize, webenv=webenv, query_key=query_key)
        yield from extract_records(records_handle)

    accession_batches = batch(accessions, batchsize)
    for acc_batch in accession_batches:
        yield from process_batch(acc_batch)


def write_record(dir, accession, record):
    with open(os.path.join(dir, accession + GB_EXT), "w") as output:
        print(record, file=output)


def main(argv):
    args = parse_args(argv)
    accessions = read_accessions(os.path.abspath(args.input))
    op_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    dbase = args.database
    Bio.Entrez.email = args.email
    batchsize = args.batch

    for acc, record in accessions_to_gb(accessions, dbase, batchsize, RETMAX):
        write_record(op_dir, acc, record)


if __name__ == "__main__":
    main(sys.argv[1:])
