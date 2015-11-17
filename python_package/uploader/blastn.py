__author__ = 'Stephen Kan'

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import sparql
from rdflib import Graph, Namespace, Literal, XSD
from ontology_uploader import upload_data
from caller_path_gen import path
from classes import generate_output

class SequenceValidator(object):
    def __init__(self):
        self.g = Graph()
        self.n =Namespace("https://github.com/superphy#")


    def validator(self):
        data = sparql.get_unvalidated_sequences()

        while data:
            for (name, sequence) in data:
                with open(path('tmp/validate.fasta'), 'w') as f:
                    f.write(sequence)

                self.blastn_commandline()
                above_threshold = self.filter_passing_hits()

                if len(above_threshold) >= 3:
                    self.g.add((self.n[name], self.n.validated, Literal("PASSED", datatype=XSD.string)))
                    print "PASSED"

                else:
                    self.g.add((self.n[name], self.n.validated, Literal("FAILED", datatype=XSD.string)))
                    print "FAILED"

                for (hit, identity) in above_threshold:
                    self.g.add((self.n[name], self.n.has_hit, Literal(str(hit), datatype=XSD.string)))

            upload_data(generate_output(self.g))
            data = sparql.get_unvalidated_sequences()


    def filter_passing_hits(self):
        result_handle = open(path("tmp/validate.xml"))
        blast_record = NCBIXML.read(result_handle)

        gene_identities = []

        for alignment in blast_record.alignments:
            hit = alignment.hit_def
            gene_length = alignment.length
            hsp = alignment.hsps[0]
            identity = (float(hsp.positives) / float(gene_length)) * 100
            gene_identities.append((hit, identity))

        above_threshold = [(x, y) for (x,y) in gene_identities if y >= 90.0]
        return above_threshold


    def blastn_commandline(self):
        blastn_cline = NcbiblastnCommandline(cmd=path("../../blast/ncbi*/bin/blastn"),
                                             query=path("tmp/validate.fasta"),
                                             db=path("data/blast/ValidationDB"),
                                             evalue=0.001, outfmt=5,
                                             out=path("tmp/validate.xml"))
        blastn_cline()


if __name__ == "__main__":
    SequenceValidator().validator()