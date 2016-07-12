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
from SuperPhy.models.upload.guelph_pipeline import Pipeline


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

@upload.route('/fasta', methods=['GET', 'POST'])
def fasta():
    basedir = os.path.realpath(os.path.dirname(__file__)).rsplit("SuperPhy", 1)[0]
    uploads = os.path.join(basedir, 'uploads')

    fasta_files = ["AYQH01000001.fasta", "JCM5491.fasta", "JEMI01000001.fasta", "KI929742.fasta", "MG1655.fasta"]
    for pos, item in enumerate(fasta_files):
        fasta_files[pos] = os.path.join(uploads, item)
        print fasta_files[pos]

    meta_file = os.path.join(uploads, "genomes.csv")
    pipeline = Pipeline(fasta_files, meta_file)
    return Response.default({"triples": pipeline.process()})




@upload.route('/foo', methods=['GET', 'POST'])
def foobar():
    """
    Endpoint for testing & devloping upload functionality
    """
    u = Uploader("JCM5491.fasta", None)
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
        self.md5sum = subprocess.Popen("md5sum {}".format(self.fasta_path),\
            shell=True, stdout=subprocess.PIPE).stdout.read(32)
        self.base_pairs = 0
        self.contig_numbers = 0
        for record in SeqIO.parse(self.fasta_path, "fasta"): #Contigs
            print record.name.rsplit('.')[0]
            #if not self.accession:
            #    self.accession = record.name.split("|")[3].split(".")[0]
            self.contig_numbers += 1
            self.base_pairs += len(record.seq)
            self.data.append((">" + self.accession, str(record.seq)))

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
