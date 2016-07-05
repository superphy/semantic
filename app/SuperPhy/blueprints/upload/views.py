"""
.
"""
import rdflib

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

@upload.route('/', methods=['GET', 'POST'])
def genome_example():
    """
    Example of making a database insertion of genome data.
    """

    graph = rdflib.Graph()
    seq = Sequence(graph, "newSequence_seq", "ATCCnewGenome", (">contig1", "ATGC",\
        ">contig2", "GGGG"), 42, 1, "fakeCheckSum", "WGS")
    seq.rdf()

    output = graph.serialize(format='turtle')

    uploader = BlazegraphUploader
    uploader.upload_data(output)

    return Response.default(graph_to_json(graph))
