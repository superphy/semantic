#!/usr/bin/env python2
# -*- coding: UTF-8 -*-

# use: python insert.py -i samples/ANLJ01.1.fsa_nt

import logging

from _utils import generate_output, generate_uri as gu, upload_data

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

def generate_turtle(graph, fasta_file, uriIsolate):
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

    TODO:
    -make a check against the db so spfyID is unique to particular isolates
    '''

    from Bio import SeqIO
    from rdflib import Literal

    # ex. :spfy234
    # rdflib.Namespace seems to not like numbers hence ge + '0001567'
    graph.add((uriIsolate, gu('rdf:type'), gu('ncbi:562')))
    # rdflib.Namespace seems to not like numbers hence ge + '0001567'
    graph.add((uriIsolate, gu('ge:0001567'), Literal("bacterium")))

    # ex. :spfy234/GCA_900089785.1_CQ10_genomic.fna
    uriAssembly = gu(uriIsolate, '/' + fasta_file.split('/')
                     [-1])  # done to ensure 1:1 for now
    # associatting isolate URI with assembly URI
    graph.add((uriIsolate, gu('g:Genome'), uriAssembly))

    # uri for bag of contigs
    # ex. :spfy234/GCA_900089785.1_CQ10_genomic.fna/contigs
    uriContigs = gu(uriAssembly, "/contigs")
    graph.add((uriAssembly, gu('so:0001462'), uriContigs))

    for record in SeqIO.parse(open(fasta_file), "fasta"):

        # ex. :spfy234/GCA_900089785.1_CQ10_genomic.fna/FLOF01006689.1
        uriContig = gu(uriAssembly, '/contigs/' + record.id)
        # linking the spec contig and the bag of contigs
        graph.add((uriContigs, gu('g:Contig'), uriContig))
        graph.add((uriContig, gu('g:DNASequence'), Literal(record.seq)))

    return graph


def call_ectyper(graph, fasta_file, uriIsolate):
    # i don't intend to import anything from ECTyper (there are a lot of
    # imports in it - not sure if we'll use them all)
    import subprocess

    from rdflib import Literal
    from ast import literal_eval
    from os.path import splitext

    logging.info('calling ectyper from fun call_ectyper')
    # concurrency is handled at the batch level, not here (note: this might change)
    # we only use ectyper for serotyping, amr is handled by rgi directly
    ectyper_dict = subprocess.check_output(['./ecoli_serotyping/src/Tools_Controller/tools_controller.py',
                                            '-in', fasta_file,
                                            '-s', '1',
                                            '-vf', '1'
                                            ])
    logging.info('inner call completed')

    # because we are using check_output, this catches any print messages from tools_controller
    # TODO: switch to pipes
    if 'error' in ectyper_dict.lower():
        logging.error('ectyper failed for' + fasta_file)
        print 'ECTyper failed for: ', fasta_file
        print 'returning graph w/o serotype'
        return graph

    logging.info('evalulating ectyper output')
    # generating the dict
    ectyper_dict = literal_eval(ectyper_dict)
    logging.info('evaluation okay')

    # we are calling tools_controller on only one file, so grab that dict
    ectyper_dict = ectyper_dict[splitext(fasta_file)[0].split('/')[-1]]

    # serotype parsing
    graph = parse_serotype(graph, ectyper_dict['Serotype'], uriIsolate)
    logging.info('serotype parsed okay')

    # amr
    graph = generate_amr(graph, uriIsolate, fasta_file)

    # vf
    graph = parse_gene_dict(
        graph, ectyper_dict['Virulence Factors'], uriIsolate, fasta_file)

    return graph


def parse_serotype(graph, serotyper_dict, uriIsolate):
    if 'O type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001076'),
                   Literal(serotyper_dict['O type'])))
    if 'H type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001077'),
                   Literal(serotyper_dict['H type'])))
    if 'K type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001684'),
                   Literal(serotyper_dict['K type'])))

    return graph


def parse_gene_dict(graph, gene_dict, uriIsolate, fasta_file):
    '''
    My intention is to eventually use ECTyper for all of the calls it was meant for.
    Just need to update ECTyper dict format to ref. AMR / VF by contig. as opposed to genome directly.

    These are the common gene related triples to both AMR / VF.
    Note: we are working from uriIsolate and assume that the calling functions (
    generate_amr() and generate_vf() are doing the transformations to the
    gene_dict.keys so that they are contig ids (as they differ in return value
    between VF & AMR from ECTyper)
    )

    TODO: offshore rgi calls to ectyper and make it return a dict in the format we need
    -currently, we'll handle ORF_ID to contig id transform in generate_amr()

    Args:
    graph(rdflib.Graph): the running graph with all our triples
    gene_dict({{}}): a dictionary of genes with a assoc info
        ex. {'Some_Contig_ID':[{'START','STOP','ORIENTATION','GENE_NAME'}]}
    uriIsolate(rdflib.URIRef): the base uri of the isolate
        ex. :spfy324

    TODO: merge common components with generate_amr()
    '''

    for contig_id in gene_dict.keys():
        for gene_record in gene_dict[contig_id]:
            # recreating the contig uri
            uriGenome = gu(':' + fasta_file.split('/')
                           [-1])
            uriContig = gu(uriGenome, '/contigs/' +
                           contig_id)  # now at contig uri
            graph.add((uriGenome, gu('so:0001462'), uriContig))

            # after this point we switch perspective to the gene and build down to
            # relink the gene with the contig

            bnode_start = BNode()
            bnode_end = BNode()

            gene_name = gene_record['GENE_NAME'].replace(' ', '_')

            graph.add((gu(':' + gene_name), gu('faldo:Begin'), bnode_start))
            graph.add((gu(':' + gene_name), gu('faldo:End'), bnode_end))

            graph.add((bnode_start, gu('rdf:type'), gu('faldo:Position')))
            graph.add((bnode_start, gu('rdf:type'), gu('faldo:ExactPosition')))
            graph.add((bnode_end, gu('rdf:type'), gu('faldo:Position')))
            graph.add((bnode_end, gu('rdf:type'), gu('faldo:ExactPosition')))

            if gene_record['ORIENTATION'] is '+':
                graph.add((bnode_start, gu('rdf:type'), gu(
                    'faldo:ForwardStrandPosition')))
                graph.add((bnode_end, gu('rdf:type'), gu(
                    'faldo:ForwardStrandPosition')))
            else:
                graph.add((bnode_start, gu('rdf:type'), gu(
                    'faldo:ReverseStrandPosition')))
                graph.add((bnode_end, gu('rdf:type'), gu(
                    'faldo:ReverseStrandPosition')))

            graph.add((bnode_start, gu('faldo:Position'),
                       Literal(gene_record['START'])))
            graph.add((bnode_start, gu('faldo:Reference'), uriContig))

            graph.add((bnode_end, gu('faldo:Position'),
                       Literal(gene_record['STOP'])))
            graph.add((bnode_end, gu('faldo:Reference'), uriContig))

            ####

    return graph


def generate_amr(graph, uriIsolate, fasta_file):
    import subprocess
    import pandas

    from os import rename
    from rdflib import BNode, Literal

    if '/' in fasta_file:
        outputname = fasta_file.split('/')[-1]
    else:
        outputname = fasta_file

    # differs from ectyper as we dont care about the temp results, just the final .tsv
    # direct (the main) call
    subprocess.call(['rgi',
                     '-i', fasta_file,
                     '-o', 'outputs/' + outputname])

    # the rgi_json call in rgitool.py isn't needed for us
    # this generates the .tsv we want
    subprocess.call(['rgi_jsontab',
                     '-i', 'outputs/' + outputname + '.json',
                     '-o', 'outputs/' + outputname])

    rename('outputs/' + outputname + '.txt', 'outputs/' + outputname + '.tsv')

    amr_results = pandas.read_table('outputs/' + outputname + '.tsv')
    amr_results = amr_results[
        ['ORF_ID', 'START', 'STOP', 'ORIENTATION', 'CUT_OFF', 'Best_Hit_ARO']]

    # triple generation
    for i in amr_results.index:

        orf_id = amr_results['ORF_ID'][i].strip()
        contig_id = orf_id.split(orf_id.split('_')[-1])[0].split('_')[0]

        # recreating the contig uri
        uriContig = gu(uriIsolate, '/' + fasta_file.split('/')
                       [-1])  # now at assembly id
        uriContig = gu(uriContig, '/contigs/' + contig_id)  # now at contig uri

        # after this point we switch perspective to the gene and build down to
        # relink the gene with the contig

        bnode_start = BNode()
        bnode_end = BNode()

        gene_name = amr_results['Best_Hit_ARO'][i].replace(' ', '_')

        graph.add((gu(':' + gene_name), gu('faldo:Begin'), bnode_start))
        graph.add((gu(':' + gene_name), gu('faldo:End'), bnode_end))

        graph.add((bnode_start, gu('dc:Description'),
                   Literal(amr_results['CUT_OFF'][i])))
        graph.add((bnode_end, gu('dc:Description'),
                   Literal(amr_results['CUT_OFF'][i])))

        graph.add((bnode_start, gu('rdf:type'), gu('faldo:Position')))
        graph.add((bnode_start, gu('rdf:type'), gu('faldo:ExactPosition')))
        graph.add((bnode_end, gu('rdf:type'), gu('faldo:Position')))
        graph.add((bnode_end, gu('rdf:type'), gu('faldo:ExactPosition')))

        if amr_results['ORIENTATION'][i] is '+':
            graph.add((bnode_start, gu('rdf:type'), gu(
                'faldo:ForwardStrandPosition')))
            graph.add((bnode_end, gu('rdf:type'), gu(
                'faldo:ForwardStrandPosition')))
        else:
            graph.add((bnode_start, gu('rdf:type'), gu(
                'faldo:ReverseStrandPosition')))
            graph.add((bnode_end, gu('rdf:type'), gu(
                'faldo:ReverseStrandPosition')))

        graph.add((bnode_start, gu('faldo:Position'),
                   Literal(amr_results['START'][i])))
        graph.add((bnode_start, gu('faldo:Reference'), uriContig))

        graph.add((bnode_end, gu('faldo:Position'),
                   Literal(amr_results['STOP'][i])))
        graph.add((bnode_end, gu('faldo:Reference'), uriContig))

        ####

    return graph

if __name__ == "__main__":

    import argparse
    import os  # for batch cleanup

    from Bio import SeqIO
    from rdflib import Namespace, BNode, Graph, URIRef, Literal
    from ConfigParser import SafeConfigParser
    from _utils import generate_uri_hash

    # setting up graph
    graph = generate_graph()

    # parsing cli-input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        help="FASTA file",
        required=True
    )
    parser.add_argument(
        "--disable-serotype",
        help="Disables use of the Serotyper. Serotyper is triggered by default.",
        action="store_true"
    )
    parser.add_argument(
        "--disable-vf",
        help="Disables use of ECTyper to get associated Virulence Factors. VFs are computed by default.",
        action="store_true"
    )
    parser.add_argument(
        "--disable-amr",
        help="Disables use of RGI to get Antimicrobial Resistance Factors.  AMR genes are computed by default.",
        action="store_true"
    )

    # note: by in large, we expect uri to be given as just the unique string value  (be it the hash or the integer) without any prefixes, the actual rdflib.URIRef object will be generated in this script

    # this is mainly for batch computation
    parser.add_argument(
        "--uri-genome",
        help="Allows the specification of the Genome URI separately. Expect just the hash (not an actual uri).",
    )
    # This is both for batch computation and for future extensions where there
    # are multiple sequencings per isolate (Campy)
    parser.add_argument(
        "--uri-isolate",
        help="Allows the specification of the Isolate URI separately. Expect just the integer (not the full :spfyID)",
        type=int
    )
    args = parser.parse_args()

    # starting logging
    logging.basicConfig(
        filename='outputs/' + __name__ + args.i.split('/')[-1] + '.log',
        level=logging.INFO
    )

    print("Importing FASTA from: " + args.i)
    logging.info('importing from' + args.i)

    # check if a genome uri isn't set yet
    if args.uri_isolate is None:
        # this is temporary, TODO: include a spqarql query to the db
        uriIsolate = gu(':spfy' + str(hash(args.i.split('/')[-1])))
    else:
        uriIsolate = gu(':spfy' + args.uri_isolate)

    # if the fasta_file hash was not precomputed (batch scripts should
    # precompute it), we compute that now
    if args.uri_genome is None:
        uriGenome = gu(':' + generate_uri_hash(args.i))
    else:
        uriGenome = gu(':' + args.uri_genome)

    # we make a dictionary from the cli-inputs and add are uris to it
    args_dict = vars(args)
    args_dict['uriIsolate'] = uriIsolate
    args_dict['uriGenome'] = uriGenome

    logging.info('generating barebones ttl from file')
    graph = generate_turtle(graph, args.i, uriIsolate)
    logging.info('barebones ttl generated')

    logging.info('calling ectyper')
    graph = call_ectyper(graph, args.i, uriIsolate)
    logging.info('ectyper call completed')

    print "Uploading to Blazegraph"
    logging.info('uploading to blazegraph')
    confirm = upload_data(generate_output(graph))
    print confirm
    logging.info(confirm)
    print 'uploaded wooot!'

    # removing fasta
    # os.remove(args.i)

    # individual fasta logs are wiped on completion (or you'd have several
    # thousand of these)
    os.remove('outputs/' + __name__ + args.i.split('/')[-1] + '.log')
