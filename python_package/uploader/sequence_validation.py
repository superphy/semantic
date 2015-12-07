#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""

"""
import string
import subprocess
import sys

from Bio.Blast import NCBIXML

from _sparql import check_checksum
from _utils import generate_path

reload(sys)
sys.setdefaultencoding("utf-8")

__author__ = "Stephen Kan"
__copyright__ = "© Copyright Government of Canada 2012-2015. Funded by the Government of Canada Genomics Research and Development Initiative"
__license__ = "ASL"
__version__ = "2.0"
__maintainer__ = "Stephen Kan"
__email__ = "stebokan@gmail.com"

class SequenceValidator(object):
    def __init__(self, seqdata):

        self.min_bp = 3500000
        self.max_bp = 7500000
        self.max_contigs = 10000

        self.seqdata = seqdata

    def validate(self):
        self.create_fasta()
        self.blastn_commandline()
        self.filter_passing_hits()

        checks = {"number of hits":self.check_hits(),
                  "base pair count":self.check_bp(),
                  "contig count":self.check_contigs(),
                  "characters": self.check_chars(),
                  "checksum":not check_checksum(self.seqdata.dict["checksum"])}

        failed_checks = {(k, v) for k, v in checks.iteritems() if v is False}

        if failed_checks:
            # replace this with logger, break would be replaced by a raised Exception where the Exception would be
            # caught by the Sequence_Upload code
            for k, v in failed_checks:
                with open(generate_path("outputs/seq_errors.txt"), "a") as f:
                    f.write("%s failed validation: the %s was not valid\n" %(self.seqdata.accession, k))
            self.seqdata.valid = False
        else:
            self.seqdata.valid = True

    def check_hits(self):
        return 3 <= len(self.seqdata.hits) <= 10

    def check_bp(self):
        return self.min_bp <= self.seqdata.dict["bp"] <= self.max_bp

    def check_contigs(self):
        return 0 < self.seqdata.dict["contigs"] <= self.max_contigs

    def check_chars(self):
        allowed_chars = "[^ACGTUNXRYSWKMBDHVacgtunxryswkmbdhv\.-]"
        s = "".join(str(contig) for contig in self.seqdata.dict["sequences"])
        trans_table = string.maketrans('','')
        return not s.translate(trans_table, allowed_chars)

    def filter_passing_hits(self):
        result_handle = open(generate_path("tmp/validate.xml"))
        hits = {}

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


    def blastn_commandline(self):
        command = generate_path("../../blast/ncbi-blast*/bin/blastn")
        fasta = generate_path("tmp/validate.fasta")
        db = generate_path("data/blast/ValidationDB")
        results = generate_path("tmp/validate.xml")

        subprocess.call("%s -query %s -db %s -outfmt 5 -out %s -best_hit_score_edge 0.05 -best_hit_overhang 0.1"
                         % (command, fasta, db, results),
                        shell=True)


    def create_fasta(self):
        with open(generate_path("tmp/validate.fasta"), "w") as f:
            for contig in self.seqdata.dict["sequences"]:
                f.write(">%s\n%s\n" %(self.seqdata.accession, contig))