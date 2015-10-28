__author__ = 'Stephen Kan'

from rdflib import Graph
from Bio import SeqIO, Entrez
import gzip
import string
from ftplib import FTP
import os
import inspect
from classes import Sequence, generate_output
from ontology_uploader import upload_data
from sparql import check_NamedIndividual, find_missing_sequences

g = Graph()
currdir = os.path.dirname(inspect.getfile(inspect.currentframe()))

def load_sequences(genome):
    name = genome + '_seq'
    print name

    if check_NamedIndividual(name):
        print name + " already in Blazegraph."

    else:
        try:
            (sequence, bp, contigs) = from_nuccore(genome)
        except ValueError:
            (sequence, bp, contigs) = from_ftp(genome)

        Sequence(g, name, genome, sequence, bp, contigs).rdf()
        upload_data(generate_output(g))


def from_nuccore(accession):
    Entrez.email = "stebokan@gmail.com"
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")

    for record in SeqIO.parse(handle, 'fasta'):
        if str(record.seq) is "":
            raise ValueError("The Genbank file is a master record with no sequence data.")
        else:
            sequence = record.seq
            contigs = 1
            bp = len(sequence)
            return (sequence, bp, contigs)


def from_ftp(accession):
    id = only_abecedarian(accession)
    filetype = 'fsa_nt.gz'

    ftp = FTP('bio-mirror.jp.apan.net')
    ftp.login('anonymous','stebokan@gmail.com')
    ftp.cwd('pub/biomirror/genbank/wgs')
    filename = get_filename(filetype, ftp, id)

    ftp.retrbinary('RETR ' + filename, open(os.path.join(currdir,'tmp/sample.gz'), 'wb').write)

    with gzip.open(os.path.join(currdir,'tmp/sample.gz')) as fasta:
        open(os.path.join(currdir,'tmp/sample.fasta'), 'wb').write(fasta.read())

    handle = open(os.path.join(currdir,'tmp/sample.fasta'), 'rU')
    sequence = ""
    contigs = 0

    for record in SeqIO.parse(handle, 'fasta'):
        sequence += str(record.seq)
        contigs += 1

    bp = len(sequence)

    return (sequence, bp, contigs)


def get_filename(filetype, ftp, id):
    filelist = ftp.nlst()
    for item in filelist:
        if id in str(item) and filetype in str(item):
            return item


def only_abecedarian(str):
    all = string.maketrans('','')
    nodigs = all.translate(all, string.ascii_letters)
    return str.translate(all, nodigs)

def upload_missing_sequences():
    for genome in find_missing_sequences():
        load_sequences(str(genome))

upload_missing_sequences()