#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#use: python insert.py -i samples/ANLJ01.1.fsa_nt

from _utils import generate_output, parse_nih_name, generate_uri, upload_data
gu = generate_uri

def generate_graph():
    '''
    Parses all the Namespaces defined in the config file and returns a graph
    with them bound.

    Return:
        (rdflib.Graph): a graph with all the defined Namespaces bound to it.
    '''

    from ConfigParser import SafeConfigParser
    from rdflib import Namespace, Graph

    graph = Graph()

    parser = SafeConfigParser()
    parser.read('config.cfg')

    for name in parser.items('Namespaces'):
        if name[0] is 'root':
            graph.bind('', Namespace(name[1]))
        else:
            graph.bind(name[0], Namespace(name[1]))

    return graph

def generate_turtle(graph, fasta_file):
    '''
    '''

    from Bio import SeqIO
    from rdflib import Literal

    for record in SeqIO.parse(open(fasta_file), "fasta"):
        identifiers = parse_nih_name(record.description)

        #creating :spfy1/ANLJ01
        #this is repetitive for the same assembly
        uriAssembly = gu(uriIsolate, '/' + identifiers['assembly'])
        #associatting isolate URI with assembly URI
        graph.add((uriIsolate, g.Genome, uriAssembly))
        graph.add((uriIsolate, gu('g:Name'), Literal(identifiers['accession_id'][0:4])))
        graph.add((uriIsolate, gu('ge:0001567'), Literal("bacterium"))) #rdflib.Namespace seems to not like numbers hence ge + '0001567'

        #no longer using blank node, instead uri for bag of contigs
        uriContigs = gu(uriAssembly + "/contigs")
        graph.add((uriAssembly, gu('so:0001462'), uriContigs))

        #creating :spfy1/ANLJ01/00001.1 ie. the contig uri
        uriContig = gu(uriAssembly, '/' + identifiers['contig'])
        graph.add((uriContigs, g.Contig, uriContig))
        graph.add((uriContig, g.DNASequence, Literal(record.seq)))

    return graph

def get_tracking_id():
    #we do this outside of record as we want same uri for all isolates
    #todo: add some check if same fasta files represents same isolate
    #grabs current id #
    parser = SafeConfigParser()
    parser.read('config.cfg')
    i = parser.getint('Database', 'id_tracking')
    return i

if __name__ == "__main__":

    import argparse
    import os #for batch cleanup

    from Bio import SeqIO
    from rdflib import Namespace, BNode, Graph, URIRef, Literal
    from ConfigParser import SafeConfigParser

    #setting up graph
    graph = generate_graph()

    #parsing cli-input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        help = "FASTA file",
    )
    args = parser.parse_args()

    print("Importing FASTA from: " + args.i)

    #we do this outside of record as we want same uri for all isolates
    #todo: add some check if same fasta files represents same isolate
    #grabs current id #
    parser = SafeConfigParser()
    parser.read('config.cfg')
    i = parser.getint('Database', 'id_tracking')
    #creating :id1
    #note: these adds become repetitive as the fasta file references the same species (will need it or a check for importing directories)
    uriIsolate = gu(':spfy' + str(i))
    i += 1
    parser.set('Database', 'id_tracking', str(i)) #saves the id so we don't overlap
    with open(r'config.cfg', 'wb') as configfile:
        parser.write(configfile)

    for record in SeqIO.parse(open(args.i), "fasta"):
        identifiers = parse_nih_name(record.description)

        #creating :spfy1/ANLJ01
        #this is repetitive for the same assembly
        uriAssembly = gu(uriIsolate, '/' + identifiers['assembly'])
        #associatting isolate URI with assembly URI
        graph.add((uriIsolate, g.Genome, uriAssembly))
        graph.add((uriIsolate, gu('g:Name'), Literal(identifiers['accession_id'][0:4])))
        graph.add((uriIsolate, gu('ge:0001567'), Literal("bacterium"))) #rdflib.Namespace seems to not like numbers hence ge + '0001567'

        #no longer using blank node, instead uri for bag of contigs
        uriContigs = gu(uriAssembly + "/contigs")
        graph.add((uriAssembly, gu('so:0001462'), uriContigs))

        #creating :spfy1/ANLJ01/00001.1 ie. the contig uri
        uriContig = gu(uriAssembly, '/' + identifiers['contig'])
        graph.add((uriContigs, g.Contig, uriContig))
        graph.add((uriContig, g.DNASequence, Literal(record.seq)))

    '''
    #for testing
    print("Writing out...")
    #we use the i value for when we're testing batches
    graph.serialize(destination='outputs/newFormat' + str(i) + '.ttl', format='turtle')
    '''

    print("Uploading to Blazegraph")
    print upload_data(generate_output(graph))
    print 'uploaded wooot!'

    os.remove(args.i)
