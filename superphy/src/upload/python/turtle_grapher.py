# a collection of turtle generation functions
# not meant to be run directly
# only contains functions related to rdflib.Graph manipulation
# only nonspecific stuff: shouldn't contain any functions directly related
# to data structure(rdf triple organization) of the modules you're dev'ing

from turtle_utils import generate_uri as gu

'''
a guide to all those predicates:
    dc:description      -human readable version for display on the web
'''


def generate_graph():
    '''
    Parses all the Namespaces defined in the config file and returns a graph
    with them bound.

    Return:
        (rdflib.Graph): a graph with all the defined Namespaces bound to it.
    '''
    import settings

    from rdflib import Namespace, Graph

    graph = Graph()

    for key in settings.namespaces.keys():
        if key is 'root':
            graph.bind('', settings.namespaces['root'])
        else:
            graph.bind(key, settings.namespaces[key])

    return graph


def generate_output(graph):
    """
    Returns RDF Graph data in the turtle format and clears the Graph

    Args:
        graph (rdflib.Graph): container object to store RDF triples
    """

    output = graph.serialize(format="turtle")
    return output


def generate_turtle_skeleton(graph, fasta_file, uriIsolate, uriGenome):
    '''
    Handles the main generation of a turtle object.

    NAMING CONVENTIONS:
    uriIsolate: this is the top-most entry, a uniq. id per file is allocated by checking our DB for the greatest most entry (not in this file)
        ex. :spfy234
    uriAssembly: aka. the genome ID, this is a sha1 hash of the file contents
        ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba
    uriContig: indiv contig ids; from SeqIO.record.id - this should be uniq to a contig (at least within a given file)
        ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba/contigs/FLOF01006689.1
        note: the record.id is what RGI uses as a prefix for ORF_ID (ORF_ID has additional _314 or w/e #s)

    Args:
        graph(rdflib.Graph): the graph instance that is 1:1 with a .fasta file
        fasta_file(str): path to the .fasta file (this should already incl the directory)
        spfyID(hash): currently a hash value generated from the name of the fasta file
    Returns:
        graph: the graph with all the triples generated from the .fasta file
    '''

    from Bio import SeqIO
    from rdflib import Literal
    from turtle_utils import uri_to_basename
    from os.path import basename

    # ex. :spfy234
    graph.add((uriIsolate, gu('rdf:type'), gu('ncbi:562')))
    graph.add((uriIsolate, gu('ge:0001567'), Literal("bacterium")))
    graph.add((uriIsolate, gu('dc:description'), Literal(uri_to_basename(uriIsolate))))

    # ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba
    # associatting isolate URI with assembly URI
    graph.add((uriIsolate, gu('g:Genome'), uriGenome))

    # this is used as the human readable display of Genome
    graph.add((uriGenome, gu('dc:description'), Literal(basename(fasta_file))))

    # uri for bag of contigs
    # ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba/contigs/
    uriContigs = gu(uriGenome, "/contigs")
    graph.add((uriGenome, gu('so:0001462'), uriContigs))

    for record in SeqIO.parse(open(fasta_file), "fasta"):

        # ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba/contigs/FLOF01006689.1
        uriContig = gu(':' + record.id)
        # linking the spec contig and the bag of contigs
        graph.add((uriContigs, gu('g:Contig'), uriContig))
        graph.add((uriContig, gu('g:DNASequence'), Literal(record.seq)))
        graph.add((uriContig, gu('g:Description'),
                   Literal(record.description)))
        graph.add((uriContig, gu('dc:description'), Literal(record.description)))

    return graph
