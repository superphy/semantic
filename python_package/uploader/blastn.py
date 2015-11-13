__author__ = 'Stephen Kan'

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import sparql
from rdflib import Graph, Namespace, Literal, XSD
from ontology_uploader import upload_data
import os
import inspect

class SequenceValidator(object):
    def __init__(self):
        self.g = Graph()
        self.n =Namespace("https://github.com/superphy#")
        self.currdir = os.path.dirname(inspect.getfile(inspect.currentframe()))


    def validator(self):
        data = sparql.get_unvalidated_sequences()

        while data:
            for (name, sequence) in data:
                with open("tmp/validate.fasta", "w") as f:
                    f.write(sequence)

                blastn_cline = NcbiblastnCommandline(cmd=self.path( "../../blast/ncbi*/bin/blastn"),
                                                     query=self.path("tmp/validate.fasta"), db=self.path("data/blast/ValidationDB"),
                                                     evalue = 0.001, outfmt=5, out=self.path("tmp/validate.xml"))
                blastn_cline()

                result_handle = open(self.path("tmp/validate.xml"))
                blast_record =  NCBIXML.read(result_handle)

                gene_identities = []

                for alignment in blast_record.alignments:
                    hsp = alignment.hsps[0]
                    if hsp.expect < 0.1:
                        gene_identities.append((float(hsp.positives) / float(hsp.align_length)) * 100)

                above_threshold = [ x for x in gene_identities if x >= 90.0]

                print name
                if len(above_threshold) >= 3:
                    self.g.add((self.n[name], self.n.validated, Literal("PASSED", datatype=XSD.string)))
                    print "PASSED"
                    upload_data(self.generate_output())
                else:
                    self.g.add((self.n[name], self.n.validated, Literal("FAILED", datatype=XSD.string)))
                    print "FAILED"
                    upload_data(self.generate_output())
            data = sparql.get_unvalidated_sequences()

    def generate_output(self):
        """
        Returns RDF Graph data in the turtle format and clears the Graph
        """

        output = self.g.serialize(format="turtle")
        self.g.remove( (None, None, None) )
        return output

    def path(self, filepath):
        return os.path.join(self.currdir, filepath)
