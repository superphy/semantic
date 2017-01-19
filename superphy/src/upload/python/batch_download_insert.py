#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#usage python batch_download_insert.py

from Bio import SeqIO

def download_to_insert(accession):
    import subprocess

    r = from_nuccore(accession)

    if r is None:
        print 'OH CRAP'
    else:
        subprocess.call('python insert.py -i ' + from_nuccore(accession))
        
    print 'woogle'

if __name__ == "__main__":
    import pandas #this is the .csv parser we're using

    from _utils import from_nuccore
    from multiprocessing import Pool

    metadata_table = pandas.read_csv('data/metadata_table.csv')
    accessions = metadata_table['primary_dbxref'].apply(lambda s: s.strip().split(':')[1])

    p = Pool(multiprocessing.cpu_count()) #you can use an int instead, just don't go crazy
    #note: you may want to write out the fasta file, but I'm unsure whether it will improve performance as concurrency requires them all to be loaded into memory anyways
    p.map(download_to_insert, accessions)
