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

def downloadFasta_to_insert(url):
    import subprocess, os

    print 'working on ' + url
    r = download_fasta(url)
    print 'done downloading, file at ' + r
    print 'now generating .ttl'
    if r is None:
        print 'OH CRAP'
    else:
        print 'calling subproces'
        subprocess.call(['./insert.py', '-i', r])
        print 'done generating .ttl'
    #os.remove(r) need to add way to check after
    print 'woogle'

if __name__ == "__main__":
    from multiprocessing import Pool, cpu_count

    '''this is testing code using the .csv file
    import pandas #this is the .csv parser we're using
    from _utils import from_nuccore

    metadata_table = pandas.read_csv('data/metadata_table.csv')
    accessions = metadata_table['primary_dbxref'].apply(lambda s: s.strip().split(':')[1])

    p = Pool(multiprocessing.cpu_count()) #you can use an int instead, just don't go crazy
    #note: you may want to write out the fasta file, but I'm unsure whether it will improve performance as concurrency requires them all to be loaded into memory anyways
    p.map(download_to_insert, accessions)
    '''

    #testing using the .txt file as source
    from _utils import download_fasta

    with open('data/download_files.txt') as f:
        lines = f.read().splitlines()
        p = Pool(cpu_count())
        p.map(downloadFasta_to_insert, lines)
