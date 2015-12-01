#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import gc
import gzip
import hashlib
import sys
import traceback
from Bio import SeqIO, Entrez
from ftplib import FTP
from rdflib import Graph
from _sparql import check_NamedIndividual, find_missing_sequences
from _utils import generate_path, generate_output, strip_non_alphabetic
from classes import Sequence
from blazegraph_upload import BlazegraphUploader
from sequence_validation import SequenceValidator

reload(sys)
sys.setdefaultencoding("utf-8")

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"


class SequenceUploader(object):
    def upload_missing_sequences(self):
        for (genome, accession) in find_missing_sequences():
            seqdata = SequenceMetadata(genome, accession)
            self.load_sequence(seqdata)
            gc.collect()

    def load_sequence(self, metadata):
        print metadata.name
        if check_NamedIndividual(metadata.name):
            print "%s already in Blazegraph." % metadata.name
        else:
            try:
                self.validated_upload(metadata)
            except TypeError:
                self.error_logging(metadata)

    def validated_upload(self, metadata):
        self.get_seqdata(metadata)
        g = Graph()
        seq_RDF = Sequence(g, **metadata.generate_kwargs())

        if metadata.is_from is not "PLASMID":
            (isValid, hits) = SequenceValidator(metadata).validate()
            if isValid:
                seq_RDF.rdf()
                seq_RDF.add_hits(hits)
            seq_RDF.add_seq_validation(isValid)
        else:
            seq_RDF.rdf()
            if metadata.accession is metadata.genome:
                seq_RDF.add_seq_validation(True)

        BlazegraphUploader().upload_data(generate_output(g))

    def error_logging(self, metadata):
        with open(generate_path("outputs/seq_errors.txt"), "a") as f:
            f.write("Genome: %s - Accession: %s.\n" % (metadata.genome, metadata.accession))
            f.write("%s \n ================================ \n\n" % traceback.format_exc())
        print "%s - %s: The records for this sequence are not retrievable." % (metadata.genome, metadata.accession)

    def get_seqdata(self, metadata):
        try:
            self.from_nuccore(metadata)
        except ValueError:
            self.from_ftp(metadata)

    def from_nuccore(self, metadata):
        Entrez.email = "superphy.info@gmail.com"
        handle = Entrez.efetch(db="nuccore", id=metadata.accession, rettype="fasta", retmode="text")

        for record in SeqIO.parse(handle, 'fasta'):
            if str(record.seq) is "":
                raise ValueError("The Genbank file is a master record with no sequence data.")
            else:
                metadata.add_sequences([record.seq])
                metadata.contigs = 1

                if "plasmid" in record.description.lower():
                    metadata.is_from = "PLASMID"
                else:
                    metadata.is_from = "CORE"

    def from_ftp(self, metadata):
        ftp = FTP('bio-mirror.jp.apan.net')
        ftp.login('anonymous', 'stebokan@gmail.com')
        ftp.cwd('pub/biomirror/genbank/wgs')
        filename = self.get_filename('fsa_nt.gz', ftp, strip_non_alphabetic(str(metadata.accession)))

        if filename:
            ftp.retrbinary('RETR ' + filename, open(generate_path('tmp/loading.gz'), 'wb').write)

            with gzip.open(generate_path('tmp/loading.gz')) as fasta:
                open(generate_path('tmp/loading.fasta'), 'wb').write(fasta.read())

            handle = open(generate_path('tmp/loading.fasta'), 'rb')

            sequences = []

            for record in SeqIO.parse(handle, 'fasta'):
                sequences.append(str(record.seq))
                metadata.contigs += 1

            metadata.add_sequences(sorted(sequences, key=len))
            metadata.is_from = "WGS"

    def get_filename(self, filetype, ftp, id):
        filelist = ftp.nlst()
        for item in filelist:
            if id in str(item) and filetype in str(item):
                return item


class SequenceMetadata(object):
    def __init__(self, genome, accession):
        self.name = accession + "_seq"
        self.accession = accession
        self.genome = genome
        self.sequences = None
        self.bp = 0
        self.contigs = 0
        self.checksum = None
        self.is_from = None

    def add_sequences(self, sequences):
        self.sequences = sequences
        self.bp = sum(len(contig) for contig in sequences)

    def generate_checksum(self):
        seqhash = hashlib.md5()
        for contig in self.sequences:
            seqhash.update(str(contig))
        self.checksum = seqhash.hexdigest()

    def generate_kwargs(self):
        self.generate_checksum()
        kwargs = {"name": self.name,
                  "genome": self.genome,
                  "sequences": self.sequences,
                  "bp": self.bp,
                  "contigs": self.contigs,
                  "checksum": self.checksum,
                  "is_from": self.is_from}
        return kwargs


if __name__ == "__main__":
    SequenceUploader().upload_missing_sequences()
