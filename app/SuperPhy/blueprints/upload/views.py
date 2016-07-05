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


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch



def getFile(filename):
    from Bio import SeqIO
    for seq_record in SeqIO.parse(filename, "fasta"): #Contigs
        print seq_record.id
        print repr(seq_record.seq)
        print len(seq_record)

    #output sum of all seq_records as bp
    #Accession #
    #Data = [seq_record.id, repr(seq_record.seq), seq_record.id, repr(seq_record.seq)]
    #contigs


def seqdata(accession, data, bp):
    graph = rdflib.Graph()
    seq = Sequence(graph, "{}{}".format(accession, "_seq"), accession, data, bp, contigs, md5sum, "WGS")
    seq.rdf()

    output = graph.serialize(format='turtle')

    uploader = BlazegraphUploader
    uploader.upload_data(output)

    return Response.default(graph_to_json(graph))



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
