"""
This view is for uploading related endpoints
"""
import subprocess
import os
import rdflib
from Bio import SeqIO
import json
import tablib

from flask import request, url_for

from SuperPhy.models import Response
from SuperPhy.models.upload.classes.sequence import Sequence
from SuperPhy.models.upload.blazegraph_upload import BlazegraphUploader
from SuperPhy.models.upload.guelph_pipeline import Pipeline
from SuperPhy.models.upload.metadata_upload import GenomeMetadataUploader

from SuperPhy.models.genome import Genomes

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

    fasta_files = ['AYQH00000000.fasta', 'JEMI00000000.fasta', 'MG1655.fasta', 'JCM5491.fasta', 'KI929742.fasta']
    for pos, item in enumerate(fasta_files):
        print os.path.join(uploads, item)
        fasta_files[pos] = os.path.join(uploads, item)
        print fasta_files[pos]

    meta_file = os.path.join(uploads, "genomes.csv")
    #pipeline = Pipeline(fasta_files, meta_file)
    
    m = GenomeMetadataUploader(os.path.join(uploads, "2_genome.json"), "Eschericia Coli")
    m.upload() #Import, not upload
    #return Response.default({"triples": pipeline.process()})
    return Response.default({})


@upload.route('/foo', methods=['GET', 'POST'])
def foobar():
    """
    Endpoint for testing & devloping upload functionality
    """
    folderpath = os.path.join(os.path.realpath(os.path.dirname(__file__)).rsplit("SuperPhy", 1)[0], 'uploads')
    """
    Download a genome with multiple contigs from:
    ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/all_assembly_versions/
    """
    fastapath = os.path.join(folderpath, "KI929759.fasta")
    
    #This is where we make rdf.
    genomes = Genomes()
    genomes.add_sequence(fastapath)
    #Fake data
    genomes.add_metadata({
        "Accession":"KI929759",
        "Bioproject_Id":"15578",
        "Biosample_Id":"2435896",
        "Genome_Uri":"https://github.com/superphy#KI929759",
        "Serotype_H":"9",
        "Serotype_O":"111",
        "Strain":"E110019"
    })

    genomes.upload()

    return genomes.data.serialize(format='turtle')

@upload.route('/post', methods=['POST', 'PUT'])
def uploading():
    """
    
    """
    #This needs to have protection!!!!!!!!!!!!!!!
    
    uploads = os.path.join(os.path.realpath(os.path.dirname(__file__)).rsplit("SuperPhy", 1)[0], 'uploads')

    files = [request.files['fasta'], request.files['meta']]
    for file_ in files:
        file_.save(os.path.join(uploads, file_.filename))
        #filename = secure_filename(file_.filename)
        #file_.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    genomes = Genomes()
    genomes.add_sequence(os.path.join(uploads, os.path.join(uploads, files[0].filename)))

    metadata = tablib.Dataset()
    metadata.csv = open(os.path.join(uploads, files[1].filename)).read()
    for ordered_dict in metadata.dict:
        genomes.add_metadata(dict(ordered_dict))

    return Response.default({"turtle":[item for item in genomes.data.serialize(format="turtle").split('\n')]})

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
