#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from rdflib import Graph
from Bio import SeqIO, Entrez
import gzip
from _utils import strip_non_alphabetic
from ftplib import FTP
import gc
from _utils import generate_path, generate_output
from classes import Sequence
from blazegraph_upload import BlazegraphUploader
from _sparql import check_NamedIndividual, find_missing_sequences
import hashlib
from sequence_validation import SequenceValidator
import traceback

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class SequenceUploader(object):
    def __init__(self):
        pass

    def upload_missing_sequences(self):
        for (genome, accession) in find_missing_sequences():
            self.load_sequences(str(genome), str(accession))
            gc.collect()


    def load_sequences(self, genome, accession):
        name = accession + '_seq'
        g = Graph()
        print name

        if check_NamedIndividual(name):
            print name + " already in Blazegraph."

        else:
            try:
                bp, contigs, is_from, sequences = self.get_seqdata(accession)
                checksum = self.generate_checksum(sequences)
                seq_RDF = Sequence(g, name, genome, sequences, bp, contigs, checksum, is_from)

                if is_from is not "PLASMID":
                    (isValid, hits) = SequenceValidator(accession, sequences, bp, contigs, checksum).validate()
                    if isValid:
                        seq_RDF.rdf()
                        seq_RDF.add_hits(hits)
                    seq_RDF.add_seq_validation(isValid)
                else:
                    seq_RDF.rdf()
                    if genome == accession:
                        seq_RDF.add_seq_validation(True)

                BlazegraphUploader().upload_data(generate_output(g))

            except TypeError:
                self.error_logging(accession, genome)

    def error_logging(self, accession, genome):
        with open(generate_path("outputs/seq_errors.txt"), "a") as f:
            f.write("%s - %s: The records for this sequence are not retrievable.\n" % (genome, accession))
            f.write(traceback.format_exc() + "\n" + "================================" + "\n" + "\n")
        print "%s - %s: The records for this sequence are not retrievable." % (genome, accession)

    def generate_checksum(self, sequences):
        seqhash = hashlib.md5()
        for contig in sequences:
            seqhash.update(str(contig))
        checksum = seqhash.hexdigest()
        return checksum

    def get_seqdata(self, accession):
        try:
            (sequences, bp, contigs, is_from) = self.from_nuccore(accession)
        except ValueError:
            (sequences, bp, contigs, is_from) = self.from_ftp(accession)
        return bp, contigs, is_from, sequences

    def from_nuccore(self, accession):
        Entrez.email = "stebokan@gmail.com"
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")

        for record in SeqIO.parse(handle, 'fasta'):
            if str(record.seq) is "":
                raise ValueError("The Genbank file is a master record with no sequence data.")
            else:
                sequences = [record.seq]
                contigs = 1
                bp = sum(len(contig) for contig in sequences)

                if "plasmid" in record.description.lower():
                    is_from = "PLASMID"
                else:
                    is_from = "CORE"

                return (sequences, bp, contigs, is_from)


    def from_ftp(self, accession):
        ftp = FTP('bio-mirror.jp.apan.net')
        ftp.login('anonymous','stebokan@gmail.com')
        ftp.cwd('pub/biomirror/genbank/wgs')
        filename = self.get_filename('fsa_nt.gz', ftp, strip_non_alphabetic(accession))

        if filename:
            ftp.retrbinary('RETR ' + filename, open(generate_path('tmp/loading.gz'), 'wb').write)

            with gzip.open(generate_path('tmp/loading.gz')) as fasta:
                open(generate_path('tmp/loading.fasta'), 'wb').write(fasta.read())

            handle = open(generate_path('tmp/loading.fasta'), 'rb')

            sequences = []
            contigs = 0

            for record in SeqIO.parse(handle, 'fasta'):
                sequences.append(str(record.seq))
                contigs += 1

            sequences = sorted(sequences, key=len)
            bp = sum(len(contig) for contig in sequences)
            is_from = "WGS"

            return (sequences, bp, contigs, is_from)


    def get_filename(self, filetype, ftp, id):
        filelist = ftp.nlst()
        for item in filelist:
            if id in str(item) and filetype in str(item):
                return item

    def save_files(self):
        pass

if __name__ == "__main__":
    SequenceUploader().upload_missing_sequences()