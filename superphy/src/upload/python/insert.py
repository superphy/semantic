#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#use: python insert.py -i samples/ANLJ01.1.fsa_nt

#importing Stephen's work for uploading to blazegraph
from blazegraph_upload import BlazegraphUploader
from _utils import generate_output, parse_nih_name, generate_uri
gu = generate_uri #shorthand to make it easier to code

if __name__ == "__main__":

    import argparse

    from Bio import SeqIO
    from rdflib import Namespace, BNode, Graph, URIRef, Literal
    from ConfigParser import SafeConfigParser

    #setting up graph
    graph = Graph()
    superphy = Namespace('https://www.github.com/superphy#')
    g = Namespace('http://www.biointerchange.org/gfvo#')
    so = Namespace('http://purl.obolibrary.org/obo/SO_')
    ge = Namespace('http://purl.obolibrary.org/obo/GENEPIO_')
    graph.bind('so', so)
    graph.bind('ge', ge)
    graph.bind('g', g)
    graph.bind('', superphy)

    #parsing cli-input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        help = "FASTA file",
    )
    args = parser.parse_args()

    print("Importing FASTA from: " + args.i)

    ##### will have to wrap entire thing for multiple files


    #we do this outside of record as we want same uri for all isolates
    #todo: add some check if same fasta files represents same isolate
    #grabs current id #
    parser = SafeConfigParser()
    parser.read('defaults.cfg')
    i = parser.getint('Database', 'id_tracking')
    #creating :id1
    #note: these adds become repetitive as the fasta file references the same species (will need it or a check for importing directories)
    uriIsolate = gu(':spfy' + str(i))
    i += 1
    parser.set('Database', 'id_tracking', str(i)) #saves the id so we don't overlap
    with open(r'defaults.cfg', 'wb') as configfile:
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

    #for testing
    print("Writing out...")
    graph.serialize(destination='outputs/newFormat.ttl', format='turtle')

    print("Uploading to Blazegraph")
    BlazegraphUploader().upload_data(generate_output(graph))

    #####will have to wrap entire thing for multiple files
