#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""This module checks to see if a retrieved sequence is valid, by checking its
base pair count, number of contigs from sequencing, and the validity of its
characters. The sequence is also ran against a reference database containing
sequences for regions that would uniquely identify it as E. coli.  As well, it
is checked for uniqueness in the database by comparing to stored md5 checksums.

Hits and validity are recorded on the SequenceMetadata object used to
instantiate the validator

"""
import string
import subprocess
import sys

from Bio.Blast import NCBIXML

from superphy.upload._sparql import check_checksum
from superphy.upload._utils import generate_path

reload(sys)

__author__ = "Stephen Kan"
__copyright__ = """
    © Copyright Government of Canada 2012-2015. Funded by the Government of
    Canada Genomics Research and Development Initiative"""
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class SequenceValidator(object):
    """ A class for evaluating if a given E. coli sequence is valid.
    """
    def __init__(self, seqdata):
        """Initializes the class with reference values for validation

        Args:
            seqdata: SequenceMetadata object containing relevant data for analysis
        """

        self.min_bp = 3500000
        self.max_bp = 7500000
        self.max_contigs = 10000

        self.seqdata = seqdata

    def validate(self):
        """Handles the whole sequence validation process. After obtaining the
        results for each check, it determines how the sequence should be
        handled in sequence uploading by modifying the associated
        SequenceMetadata object.

        TODO: refactor this more for clarity and ease of testing?
        """
        self.filter_passing_hits()

        checks = {"number of hits":self.check_hits(),
                  "base pair count":self.check_bp(),
                  "contig count":self.check_contigs(),
                  "characters": self.check_chars(),
                  "checksum":not check_checksum(self.seqdata.checksum)}

        failed_checks = {(k, v) for k, v in checks.iteritems() if v is False}

        if failed_checks:
            """
            replace this with logger, break would be replaced by a raised
            Exception where the Exception would be caught by the
            Sequence_Upload code
            """
            for k, v in failed_checks:
                with open(generate_path("outputs/seq_errors.txt"), "a") as file_:
                    file_.write(
                        '%s failed validation:'
                        'the %s was not valid\n' %(self.seqdata.accession, k)
                    )
            self.seqdata.valid = False
        else:
            self.seqdata.valid = True

    def check_hits(self):
        """
        Checks if there are sufficient qualifying hits against the test
        database to qualify as a valid genome. If there are more hits than
        regions in the database, the sequence is also disqualified
        (this is probably more indicative of a programming issue)

        Returns: a boolean indicating if a sequence passes this check
        """
        return 3 <= len(self.seqdata.hits) <= 10

    def check_bp(self):
        """
        Checks if the number of base pairs in the sequence is valid
        (enough to encompass a mostly intact E. coli sequence, not so much that
        it would indicate contamination of the sample with foreign data or
        improper alignment construction generating artifact sequences)

        Returns: a boolean indicating if a sequence passes this check
        """
        return self.min_bp <= self.seqdata.bp <= self.max_bp

    def check_contigs(self):
        """
        Checks if the number of contigs composing the sequence is valid
        (for WGS samples in particular). If the contiq number is too high, the
        quality of the alignment is too poor to be of use for the project.

        Returns: a boolean indicating if a sequence passes this check
        """
        return 0 < self.seqdata.numcontigs <= self.max_contigs

    def check_chars(self):
        """Checks if the characters composing the sequence are valid: they only
        contain IUPAC codes for nucleotides, including uncertain nucleotides.

        Returns: a boolean indicating if a sequence passes this check
        """
        allowed_chars = r"[^ACGTUNXRYSWKMBDHVacgtunxryswkmbdhv\.-]"
        s = "".join(str(seq) for (accession_name, seq) in self.seqdata.contigs)
        trans_table = string.maketrans('', '')
        return not s.translate(trans_table, allowed_chars)

    def filter_passing_hits(self):
        """
        Reads the result from the command line BLAST using fileIO and parses it
        to look for the top scoring hits at 90% and above. If there are
        multiple hits, select the highest scoring one.
        """
        self.create_fasta()
        self.blastn_commandline()

        hits = {}
        result_handle = open(generate_path("tmp/validate.xml"))
        for record in NCBIXML.parse(result_handle):
            for entry in record.alignments:
                hit = entry.hit_def
                seqlen = entry.length
                hsp = entry.hsps[0]
                percent_ident = (float(hsp.positives) / float(seqlen)) * 100

                if 90 <= percent_ident <= 100:
                    if hit in hits:
                        if percent_ident > hits[hit]:
                            hits[hit] = percent_ident
                    else:
                        hits[hit] = percent_ident
        del result_handle
        self.seqdata.hits = hits
    @classmethod
    def blastn_commandline(cls):
        """Runs a command line BLAST on the generated FASTA sequence using the
        database composed of 10 E. coli species-specific genomic regions and
        outputs the results into XML format into another file.
        """
        command = generate_path("../../blast/ncbi-blast*/bin/blastn")
        fasta = generate_path("tmp/validate.fasta")
        db = generate_path("data/blast/ValidationDB")
        results = generate_path("tmp/validate.xml")

        subprocess.call(
            '%s -query %s -db %s -outfmt 5 -out %s -best_hit_score_edge 0.05 '
            '-best_hit_overhang 0.1' % (
                command, fasta, db, results
            ), shell=True
        )

    def create_fasta(self):
        """Writes a FASTA sequence to a file for use by the command line version of BLAST. Obtains nucleotide data from
        the sequence data object used to initialize the validator and writes each entry as a separate FASTA object.
        Contigs from WGS samples must be kept separate to avoid false matches based on misaligned sequences.
        """
        with open(generate_path("tmp/validate.fasta"), "w") as file_:
            for (accession_name, seq) in self.seqdata.contigs:
                file_.write(">%s\n%s\n" %(self.seqdata.accession, seq))