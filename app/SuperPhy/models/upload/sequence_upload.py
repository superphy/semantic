#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""This module looks for eligible genomes and checks NCBI Genbank and WGS
databases for relevant sequence metadata, such as FASTA files and descriptions.
Genomes are eligible if they do not already have a valid core or WGS sequence
associated with them, or if they have already been processed and found to have
an invalid genome.

TODO: parallelize this along with the rest of the upload/validation scripts.

Classes:
    SequenceUploader: retrieves and uploads sequence metadata of eligible
    genomes.
    SequenceMetadata: stores sequence metadata for uploading
"""

from ftplib import FTP
import gc
import gzip
import hashlib
import sys
import traceback
from urllib2 import HTTPError

from Bio import Entrez, SeqIO
from rdflib import Graph

from SuperPhy.models.upload._sparql import check_named_individual, \
    find_missing_sequences
from SuperPhy.models.upload._utils import generate_output, generate_path, \
    strip_non_alphabetic
from SuperPhy.models.upload.classes import Sequence
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader
from SuperPhy.models.upload.sequence_validation import SequenceValidator

reload(sys)
sys.setdefaultencoding("utf-8")

class SequenceUploader(object):
    """A class for retrieving and uploading sequence metadata of eligible
    genomes
    """
    def upload_missing_sequences(self):
        """Compiles a list of genomes with missing and unvalidated sequences
        and uploads them to Blazegraph.
        """
        for (genome, accession) in find_missing_sequences():
            seqdata = SequenceMetadata(genome, accession)
            try:
                self.load_sequence(seqdata)
                if seqdata.dict["is_from"] == "PLASMID":
                    self.upload(seqdata, self.plasmid_rdf)
                else:
                    SequenceValidator(seqdata).validate()
                    self.upload(seqdata, self.nonplasmid_rdf)
                gc.collect()
            except TypeError:
                self.error_logging(seqdata)

    def load_sequence(self, seqdata):
        """Checks to see if the sequence is already uploaded onto Blazegraph,
        and if not, try loading the sequence.
        Logs all TypeErrors (assumption: all other errors are mistakes in the
        process and not something that requires
        manual curation).

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        print seqdata.name
        if check_named_individual(seqdata.name):
            print "%s already in Blazegraph." % seqdata.name
            raise TypeError
        else:
            self.get_seqdata(seqdata)
    @classmethod
    def upload(cls, seqdata, func):
        """Uploads sequence data to Blazegraph based on the inputted function
        argument

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump:
            func: function that handles RDF conversion based on the type of
            genome
        """
        graph = Graph()
        seq_rdf = Sequence(graph, **seqdata.generate_kwargs())

        func(seqdata, seq_rdf)

        BlazegraphUploader().upload_data(generate_output(graph))

    @classmethod
    def plasmid_rdf(cls, seqdata, seq_rdf):
        """Sets up RDF triples for a plasmid sequence and its metadata. Marks
        genome as possessing a valid sequence if
        the accession id matches the genome id (i.e. metadata validation has
        not merged genomes with the same biosample
        id together yet)

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
            seq_rdf: an initalized RDF converter for Sequence data
        """
        seq_rdf.rdf()

        if seqdata.accession == seqdata.genome:
            seq_rdf.add_seq_validation(True)

    @classmethod
    def nonplasmid_rdf(cls, seqdata, seq_rdf):
        """Sets up RDF triples for a nonplasmid sequence and its metadata in
        accordance with its sequence validation
         results

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
            seq_rdf: an initalized RDF converter for Sequence data
        """
        seq_rdf.add_seq_validation(seqdata.valid)

        if seqdata.valid:
            seq_rdf.rdf()
            seq_rdf.add_hits(seqdata.hits)

    #should replace with logger module
    @classmethod
    def error_logging(cls, seqdata):
        """Logs errors regarding sequence uploading to a file, for manual
        curation.

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        with open(generate_path("outputs/seq_errors.txt"), "a") as file_:
            file_.write(
                "Genome: %s - Accession: %s.\n" % (
                    seqdata.genome, seqdata.accession
                )
            )
            file_.write("%s \n ================================ \n\n" % traceback.format_exc())
        print "%s - %s: The records for this sequence are not retrievable." % (
            seqdata.genome, seqdata.accession
            )

    def get_seqdata(self, seqdata):
        """Tries to get the FASTA sequence, first through the NCBI Genbank
        Nucleotide database, and then in the WGS
        pipeline directory of the NCBI FTP server.

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        try:
            self.from_nuccore(seqdata)
        except ValueError:
            self.from_ftp(seqdata)

    def from_nuccore(self, seqdata):
        """Obtains the FASTA sequence via the NCBI Genbank Nucleotide database
        using Entrez EUtils. If the sequence is not already annotated as a
        plasmid, annotate on the assumption that it is a core genome. If there
        is nothing found for the sequence, raise a ValueError.

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        Entrez.email = "superphy.info@gmail.com"
        handle = None
        i = 0

        while i < 3:
            try:
                handle = Entrez.efetch(
                    db="nuccore",
                    id=seqdata.accession,
                    rettype="fasta",
                    retmode="text"
                )
                self.read_fasta(handle, seqdata)

                if seqdata.dict["sequences"] == ['']:
                    raise ValueError(
                        'The Genbank file is a master record withno sequence '
                        'data.'
                        )

                if seqdata.dict["is_from"] is None:
                    seqdata.dict["is_from"] = "CORE"
                break
            except HTTPError:
                i += 1
                continue
        try:
            handle is None
        except NameError:
            raise TypeError("Could not retrieve file for analysis")

    def from_ftp(self, seqdata):
        """Obtains the FASTA sequence via the NCBI FTP server in the WGS genome
        pipeline and labels the sequence as being from the WGS piepline.

        Args:
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        seq_id = strip_non_alphabetic(str(seqdata.accession))
        self.download_file(seq_id, 'fsa_nt.gz')

        with open(generate_path('tmp/loading.fasta'), 'rb') as handle:
            self.read_fasta(handle, seqdata)
        seqdata.dict["is_from"] = "WGS"
    @classmethod
    def read_fasta(cls, handle, seqdata):
        """Reads the file object containing a properly formatted FASTA document
        and retrieves data to store in a SequenceMetadata instance.

        Args:
            handle: a file instance containing a properly formatted FASTA
            document
            seqdata: a SequenceMetadata instance storing sequence-related data
            that would otherwise be a data clump
        """
        sequences = []

        for record in SeqIO.parse(handle, 'fasta'):
            sequences.append(str(record.seq))

            if "plasmid" in record.description.lower():
                seqdata.dict["is_from"] = 'PLASMID'

        seqdata.add_sequences(sorted(sequences, key=len))
    @classmethod
    def download_file(cls, id_, filetype):
        """Downloads the gzip file with the correct id and filetype and unzips
        it and transfers its contents into a temporary FASTA file for further
        processing. If no files on the server match, returns a TypeError.

        Args:
            id(str): a WGS project ID, composed of only alphabetics
            filetype(str): the type of file to be found. 'fsa_nt.gz' is the
            default, but there are other options for amino acids and other
            formats
        """
        ftp = FTP('bio-mirror.jp.apan.net')
        ftp.login('anonymous', 'superphy.info@gmail.com')
        ftp.cwd('pub/biomirror/genbank/wgs')

        filenames = ftp.nlst()
        filename = [s for s in filenames if id_ in s and filetype in s]

        if len(filename) is not 1:
            raise TypeError("No files could be found for download.")
        else:
            ftp.retrbinary(
                'RETR ' + filename[0],
                open(generate_path('tmp/loading.gz'), 'wb').write
            )
            with gzip.open(generate_path('tmp/loading.gz')) as fasta, \
                    open(generate_path('tmp/loading.fasta'), 'wb') as output:
                output.write(fasta.read())


class SequenceMetadata(object):
    """A class for storing sequence metadata for uploading into Blazegraph.
    """
    def __init__(self, genome, accession):
        """Initializes the class with the necessary fields. The dict is used to
        construct the kwargs to pass into
        classes.Sequences for uploading into Blazegraph

        Args:
            genome(str): the accession of the genome that the sequence is
            associated with
            accession(str): the accession that the sequence is identified by
        """
        self.name = str(accession) + "_seq"
        self.accession = str(accession)
        self.genome = str(genome)
        self.hits = None
        self.valid = None

        keys = ["name", "genome", "sequences", "bp", "contigs", "checksum",
                "is_from"]
        self.dict = {key: None for key in keys}
        self.dict["name"] = self.name
        self.dict["genome"] = self.genome

    def add_sequences(self, sequences):
        """Add the sequence FASTA to the class, and calculate metadata
        attributes for the sequence such as the number of base pairs, the
        number of contigs, and the MD5 checksum used to check identity

        Args:
            sequences: a list of sequences(str) to add to the class
        """
        self.dict["sequences"] = sequences
        self.dict["contigs"] = len(sequences)
        self.dict["bp"] = sum(len(contig) for contig in sequences)
        self.generate_checksum()

    def generate_checksum(self):
        """Generates an MD5 checksum from the sorted contigs in the sequence
        (sorted from smallest to largest; order matters).
        """
        seqhash = hashlib.md5()
        for contig in self.dict["sequences"]:
            seqhash.update(str(contig))
        self.dict["checksum"] = seqhash.hexdigest()

    def generate_kwargs(self):
        """Returns a dict of kwargs to pass into the classes.Sequence
        constructor. If any dict entries are NoneType, raise an error
        indicating the first found that is missing.
        """
        for key, value in self.dict.iteritems():
            if not value:
                raise TypeError("Missing sequence metadata: %s" % key)
        return self.dict

if __name__ == "__main__":
    SequenceUploader().upload_missing_sequences()
