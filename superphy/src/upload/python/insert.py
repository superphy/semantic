#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#use: python insert.py -i samples/ANLJ01.1.fsa_nt

if __name__ == "__main__":

    import argparse

    from Bio import SeqIO
    from rdflib import Namespace, BNode, Graph, URIRef, Literal

    #setting up graph
    graph = Graph()
    superphy = 'https://github.com/superphy#'
    g = Namespace("http://www.biointerchange.org/gfvo#")
    so = 'http://purl.obolibrary.org/obo/SO_'
    ge = 'http://purl.obolibrary.org/obo/GENEPIO_'

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
        species = record.description.split("|")[4].split(" ")[3]
        assembly = accession_id[0:6]
        contig = record.description.split("|")[4].split(" ")[4].split(".")[2]

        #creating :id1
        #note: these adds become repetitive as the fasta file references the same species (will need it or a check for importing directories)
        graph.add((URIRef(superphy + "idEcoli" + species), g.Name, Literal(accession_id[0:4])))
        graph.add((URIRef(superphy + "idEcoli" + species), URIRef(ge + '0001567'), Literal("bacterium"))) #rdflib.Namespace seems to not like numbers hence ge + '0001567'

        #creating :genomeID1
        #this is repetitive for the same assembly
        graph.add((URIRef(superphy + "idEcoli" + species), g.Genome, URIRef(superphy + "genomeID" + assembly)))

        #no longer using blank node, instead uri
        graph.add((URIRef(superphy + "genomeID" + assembly), URIRef(so + '0001462'), URIRef(superphy + "genomeID" + assembly + "contigs")))

        #creating :genomeID1/contigID1
        graph.add((URIRef(superphy + "genomeID" + assembly + "contigs"), g.Contig, URIRef(superphy + "genomeID" + assembly + '/contigID' + contig)))
        graph.add((URIRef(superphy + "genomeID" + assembly + '/contigID' + contig), g.DNASequence, Literal(record.seq)))

    print("Writing out...")
    graph.serialize(destination='outputs/newFormat.txt', format='turtle')
