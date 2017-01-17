#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#use: python insert.py -i samples/ANLJ01.1.fsa_nt

#importing Stephen's work for uploading to blazegraph
from blazegraph_upload import BlazegraphUploader
from _utils import generate_output, parse_nih_name

if __name__ == "__main__":

    import argparse

    from Bio import SeqIO
    from rdflib import Namespace, BNode, Graph, URIRef, Literal

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
    for record in SeqIO.parse(open(args.i), "fasta"):
        accession_id = record.id.split("|")[3].split(".")[0]
        identifiers = parse_nih_name(record.description)
        print d
        species = record.description.split("|")[4].split(" ")[3]
        assembly = accession_id[0:6]
        contig = accession_id[6:12]

        #creating :id1
        #note: these adds become repetitive as the fasta file references the same species (will need it or a check for importing directories)
        uriSpecies = superphy["idEcoli" + species]
        graph.add((uriSpecies, g.Name, Literal(accession_id[0:4])))
        graph.add((uriSpecies, ge['0001567'], Literal("bacterium"))) #rdflib.Namespace seems to not like numbers hence ge + '0001567'

        #creating :genomeID1
        #this is repetitive for the same assembly
        uriAssembly = superphy["genomeID" + assembly]
        graph.add((superphy["idEcoli" + species], g.Genome, uriAssembly))

        #no longer using blank node, instead uri
        uriContigs = superphy["genomeID" + assembly + "/contigs"]
        graph.add((uriAssembly, so['0001462'], uriContigs))

        #creating :genomeID1/contigID1
        uriContig = superphy["genomeID" + assembly + '/contigID' + contig]
        graph.add((uriContigs, g.Contig, uriContig))
        graph.add((uriContig, g.DNASequence, Literal(record.seq)))

    print("Writing out...")
    graph.serialize(destination='outputs/newFormat.ttl', format='turtle')

    print("Uploading to Blazegraph")
    BlazegraphUploader().upload_data(generate_output(graph))
