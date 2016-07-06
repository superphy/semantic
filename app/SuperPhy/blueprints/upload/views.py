"""
This view is for uploading related endpoints
"""
import subprocess
import os
import rdflib
from Bio import SeqIO

from SuperPhy.models import Response
from SuperPhy.models.upload.classes.sequence import Sequence
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader

from SuperPhy.blueprints.upload import upload

def graph_to_json(g):
    """
    Pass in a rdflib.Graph and get back a chunk of JSON using
    the Talis JSON serialization for RDF:
    Source: https://gist.github.com/edsu/76729
    """
    json = {}

    # go through all the triples in the graph
    for s, p, o in g:

        # initialize property dictionary if we've got a new subject
        if not json.has_key(s):
            json[s] = {}

        # initialize object list if we've got a new subject-property combo
        if not json[s].has_key(p):
            json[s][p] = []

        # determine the value dictionary for the object
        v = {'value': unicode(o)}
        if isinstance(o, rdflib.URIRef):
            v['type'] = 'uri'
        elif isinstance(o, rdflib.BNode):
            v['type'] = 'bnode'
        elif isinstance(o, rdflib.Literal):
            v['type'] = 'literal'
            if o.language:
                v['lang'] = o.language
            if o.datatype:
                v['datatype'] = unicode(o.datatype)

        # add the triple
        json[s][p].append(v)
    return json

@upload.route('/foo', methods=['GET', 'POST'])
def foobar():
    """
    Endpoint for testing & devloping upload functionality
    """
    u = Uploader("foo.fasta", None)
    data = u.send_fasta_data()
    return Response.default(data)

class Uploader(object):
    """
    Class for moving fasta files and meta-data to blazegraph.
    """
    def __init__(self, fasta_file, meta_file):
        basedir = os.path.realpath(os.path.dirname(__file__)).rsplit("SuperPhy", 1)[0]
        filepath = os.path.join(basedir, 'uploads', fasta_file)

        self.fasta_path = filepath
        self.meta_file = meta_file

        self.accession = ""
        self.data = []
        self.base_pairs = 0
        self.contig_numbers = 0
        self.md5sum = ""
        self.get_fasta_data()

    def fake_fasta_data(self):
        """
        Fake fasta data in the format provided.
        """
        self.accession = "newSequence"
        self.accession = "ATCCnewGenome"
        self.data = (">contig1", "ATGC", ">contig2", "GGGG")
        self.base_pairs = 42
        self.contig_numbers = 1
        self.md5sum = "fakeCheckSum"

    def get_fasta_data(self):
        """
        Open the fasta file and get the fasta data.
        """
        base_pairs = 0
        contigs = 0
        data = []
        for record in SeqIO.parse(self.fasta_path, "fasta"): #Contigs
            contigs += 1
            base_pairs += len(record)
            try:
                id_ = record.id.split(">", 1)[1].split(" ", 1)[0]
            except IndexError:
                id_ = record.id
            data.append(id_)
            data.append(repr(record.seq))


            #print record.seq
            #print repr(record.seq)


        self.md5sum = subprocess.Popen("md5sum {}".format(self.fasta_path),\
            shell=True, stdout=subprocess.PIPE).stdout.read(32)
        self.accession = "???"
        self.data = data
        self.base_pairs = base_pairs
        self.contig_numbers = contigs

    def send_fasta_data(self):
        """
        Uploads and sends the fasta data to Blazegraph
        """
        graph = rdflib.Graph()
        seq = Sequence(graph, "{}{}".format(self.accession, "_seq"),\
            self.accession, self.data, self.base_pairs, self.contig_numbers,\
            self.md5sum, "WGS")
        seq.rdf()

        output = graph.serialize(format='turtle')

        uploader = BlazegraphUploader
        uploader.upload_data(output)

        return graph_to_json(graph)

'''
@upload.route('/', methods=['GET', 'POST'])
def genome_example():
    """
    Example of making a database insertion of genome data.
    """
    graph = rdflib.Graph()
    seq = Sequence(graph, "newSequence{}".format("_seq"), "ATCCnewGenome", data, 42, 1, "fakeCheckSum", "WGS")
    seq.rdf()

    output = graph.serialize(format='turtle')

    uploader = BlazegraphUploader
    uploader.upload_data(output)

    return Response.default(graph_to_json(graph))

'''
