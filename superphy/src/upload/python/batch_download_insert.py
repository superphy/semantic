#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#usage python batch_download_insert.py

def download_to_insert(accession):
    import subprocess

    from Bio import SeqIO

    SeqIO.write(from_nuccore(accession), 'tmp/' + accession + '.fasta', 'fasta')
    subprocess.call('python insert.py -i ' + 'tmp/' + accession + '.fasta')
        #print 'working on ' + record.id
        #SeqIO.write(record, 'tmp/' + accession)

def test(a):
    print a
    print "boogle"

if __name__ == "__main__":
    import pandas #this is the .csv parser we're using

    from _utils import from_nuccore
    from multiprocessing import Pool

    metadata_table = pandas.read_csv('data/metadata_table.csv')
    accessions = metadata_table['primary_dbxref'].apply(lambda s: s.strip().split(':')[1])

    p = Pool(4) #don't go crazy here
    p.map(test, accessions) #note: you may want to write out the fasta file, but I'm unsure whether it will improve performance as concurrency requires them all to be loaded into memory anyways
    p.map(download_to_insert, accessions)
