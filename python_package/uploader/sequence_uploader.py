__author__ = 'Stephen Kan'

from rdflib import Graph
from Bio import SeqIO, Entrez
import gzip
import string
from ftplib import FTP
import os
import gc
import inspect
from classes import Sequence, generate_output
from ontology_uploader import upload_data
from sparql import check_NamedIndividual, find_missing_sequences

class SequenceUploader(object):
    def __init__(self):
        pass

    def load_sequences(self, genome, accession):
        name = accession + '_seq'
        print name

        if check_NamedIndividual(name):
            print name + " already in Blazegraph."

        else:
            try:
                g = Graph()
                try:
                    (sequence, bp, contigs, is_from) = self.from_nuccore(accession)
                    sequence = Sequence(g, name, genome, sequence, bp, contigs)
                    sequence.rdf()
                    sequence.add_is_from(is_from)
                except ValueError:
                    (sequence, bp, contigs) = self.from_ftp(accession)
                    sequence = Sequence(g, name, genome, sequence, bp, contigs)
                    sequence.rdf()
                    sequence.add_is_from("WGS")

                upload_data(generate_output(g))
            except TypeError:
                f = open(os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), "outputs/seq_errors.txt"), "a")
                f.write("%s - %s: The records for this sequence are not retrievable." %(genome, accession))
                print "%s - %s: The records for this sequence are not retrievable." %(genome, accession)


    def from_nuccore(self, accession):
        Entrez.email = "stebokan@gmail.com"
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")

        for record in SeqIO.parse(handle, 'fasta'):
            if str(record.seq) is "":
                raise ValueError("The Genbank file is a master record with no sequence data.")
            else:
                sequence = record.seq
                contigs = 1
                bp = len(sequence)

                if "plasmid" in record.description.lower():
                    is_from = "PLASMID"
                else:
                    is_from = "CORE"

                return (sequence, bp, contigs, is_from)


    def from_ftp(self, accession):
        currdir = os.path.dirname(inspect.getfile(inspect.currentframe()))

        ftp = FTP('bio-mirror.jp.apan.net')
        ftp.login('anonymous','stebokan@gmail.com')
        ftp.cwd('pub/biomirror/genbank/wgs')
        filename = self.get_filename('fsa_nt.gz', ftp, self.only_abecedarian(accession))
        ftp.retrbinary('RETR ' + filename, open(os.path.join(currdir,'tmp/sample.gz'), 'wb').write)

        with gzip.open(os.path.join(currdir,'tmp/sample.gz')) as fasta:
            open(os.path.join(currdir,'tmp/sample.fasta'), 'wb').write(fasta.read())

        handle = open(os.path.join(currdir,'tmp/sample.fasta'), 'rb')

        sequence = ""
        contigs = 0

        for record in SeqIO.parse(handle, 'fasta'):
            sequence += str(record.seq)
            contigs += 1

        bp = len(sequence)

        return (sequence, bp, contigs)


    def get_filename(self, filetype, ftp, id):
        filelist = ftp.nlst()
        for item in filelist:
            if id in str(item) and filetype in str(item):
                return item


    def only_abecedarian(self, str):
        all = string.maketrans('','')
        nodigs = all.translate(all, string.ascii_letters)
        return str.translate(all, nodigs)


    def upload_missing_sequences(self):
        for (genome, accession) in find_missing_sequences():
            self.load_sequences(str(genome), str(accession))
            gc.collect()

if __name__ == "__main__":
    SequenceUploader().upload_missing_sequences()