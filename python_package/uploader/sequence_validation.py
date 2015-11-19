__author__ = 'Stephen Kan'

from Bio.Blast import NCBIXML
from _utils import path
import subprocess

"""
SearchIO is still experimental, but there are no expected API changes. The file format I am using,
BLAST+ tabular format, is not expected to change.

See http://comments.gmane.org/gmane.comp.python.bio.general/8598
"""


class SequenceValidator(object):
    def __init__(self, genome, sequence):
        pass

    def validator(self):
        self.create_fasta()
        self.blastn_commandline()
        hits = self.filter_passing_hits()

        if len(hits) < 3:
            pass
        else:
            pass


    def filter_passing_hits(self):
        result_handle = open(path("tmp/results.xml"))

        hits = {}

        for record in NCBIXML.parse(result_handle):
            for entry in record.alignments:
                hit = entry.hit_def
                seqlen = entry.length
                hsp = entry.hsps[0]
                percent_ident = (float(hsp.positives) / float(seqlen)) * 100

                if percent_ident >= 90:
                    if hit in hits:
                        if percent_ident > hits[hit]:
                            hits[hit] = percent_ident
                    else:
                        hits[hit] = percent_ident

        del result_handle
        return hits

    def blastn_commandline(self):
        command = path("../../blast/ncbi-blast*/bin/blastn")
        fasta = path("tmp/validate.fasta")
        db = path("data/blast/ValidationDB")
        results = path("tmp/validate.xml")

        subprocess.call("%s -query %s -db %s -outfmt 6 -out %s -best_hit_score_edge 0.05 -best_hit_overhang 0.1"
                         % (command, fasta, db, results),
                        shell=True)

    def create_fasta(self):
        with open(path("tmp/validate.fasta")) as f:
            for contig in self.sequence:
                f.write(">%s\n%s" %(self.genome, contig))

if __name__ == "__main__":
    print SequenceValidator().filter_passing_hits()