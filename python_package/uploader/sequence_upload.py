#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from ftplib import FTP
import gc
import gzip
import hashlib
import sys
import traceback

from Bio import Entrez, SeqIO
from rdflib import Graph

from _sparql import check_NamedIndividual, find_missing_sequences
from _utils import generate_output, generate_path, strip_non_alphabetic
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

    def load_sequence(self, seqdata):
        print seqdata.name
        if check_NamedIndividual(seqdata.name):
            print "%s already in Blazegraph." % seqdata.name
        else:
            try:
                self.validated_upload(seqdata)
            except TypeError:
                self.error_logging(seqdata)

    def validated_upload(self, seqdata):
        self.get_seqdata(seqdata)
        g = Graph()
        seq_rdf = Sequence(g, **seqdata.generate_kwargs())

        if seqdata.dict["is_from"] is not "PLASMID":
            (isValid, hits) = SequenceValidator(seqdata).validate()
            if isValid:
                seq_rdf.rdf()
                seq_rdf.add_hits(hits)
            seq_rdf.add_seq_validation(isValid)
        else:
            seq_rdf.rdf()
            if seqdata.accession is seqdata.genome:
                seq_rdf.add_seq_validation(True)

        BlazegraphUploader().upload_data(generate_output(g))

    def error_logging(self, seqdata):
        with open(generate_path("outputs/seq_errors.txt"), "a") as f:
            f.write("Genome: %s - Accession: %s.\n" % (seqdata.genome, seqdata.accession))
            f.write("%s \n ================================ \n\n" % traceback.format_exc())
        print "%s - %s: The records for this sequence are not retrievable." % (seqdata.genome, seqdata.accession)

    def get_seqdata(self, seqdata):
        try:
            self.from_nuccore(seqdata)
        except ValueError:
            self.from_ftp(seqdata)

    def from_nuccore(self, seqdata):
        Entrez.email = "superphy.info@gmail.com"
        handle = Entrez.efetch(db="nuccore", id=seqdata.accession, rettype="fasta", retmode="text")

        self.read_fasta(handle, seqdata)

        if seqdata.dict["sequences"] == ['']:
            raise ValueError("The Genbank file is a master record with no sequence data.")

        if seqdata.dict["is_from"] is None:
            seqdata.dict["is_from"] = "CORE"

    def from_ftp(self, seqdata):
        seq_id = strip_non_alphabetic(str(seqdata.accession))
        self.download_file(seq_id, 'fsa_nt.gz')

        with open(generate_path('tmp/loading.fasta'), 'rb') as handle:
            self.read_fasta(handle, seqdata)
        seqdata.dict["is_from"] = "WGS"

    def read_fasta(self, handle, seqdata):
        sequences = []

        for record in SeqIO.parse(handle, 'fasta'):
            sequences.append(str(record.seq))

            if "plasmid" in record.description.lower():
                seqdata.dict["is_from"] = 'PLASMID'

        seqdata.add_sequences(sorted(sequences, key=len))

    def download_file(self, id, filetype):
        ftp = FTP('bio-mirror.jp.apan.net')
        ftp.login('anonymous', 'superphy.info@gmail.com')
        ftp.cwd('pub/biomirror/genbank/wgs')

        filenames = ftp.nlst()
        filename = [s for s in filenames if id in s and filetype in s]

        if len(filename) is not 1:
            raise TypeError("No files could be found for download.")
        else:
            ftp.retrbinary('RETR ' + filename[0], open(generate_path('tmp/loading.gz'), 'wb').write)
            with gzip.open(generate_path('tmp/loading.gz')) as fasta, \
                    open(generate_path('tmp/loading.fasta'), 'wb') as output:
                output.write(fasta.read())


class SequenceMetadata(object):
    def __init__(self, genome, accession):
        self.name = accession + "_seq"
        self.accession = accession
        self.genome = genome
        self.hits = None

        keys = ["name", "genome", "sequences", "bp", "contigs", "checksum", "is_from"]
        self.dict = {key: None for key in keys}
        self.dict["name"] = accession + "_seq"
        self.dict["genome"] = genome

    def add_sequences(self, sequences):
        self.dict["sequences"] = sequences
        self.dict["contigs"] = len(sequences)
        self.dict["bp"] = sum(len(contig) for contig in sequences)
        self.generate_checksum()

    def generate_checksum(self):
        seqhash = hashlib.md5()
        for contig in self.dict["sequences"]:
            seqhash.update(str(contig))
        self.dict["checksum"] = seqhash.hexdigest()

    def generate_kwargs(self):
        for key, value in self.dict.iteritems():
            if not value:
                raise TypeError("Missing sequence metadata: %s" % key)
        return self.dict


if __name__ == "__main__":
    SequenceUploader().upload_missing_sequences()
