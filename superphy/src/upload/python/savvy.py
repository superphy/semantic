#!/usr/bin/env python2
# -*- coding: UTF-8 -*-

# use: python savvy.py -i samples/ANLJ01.1.fsa_nt

# S:erotype
# A:ntimicrobial Resistance
# V:irulence Factors
# -vy

import logging
# long function calls, cause datastruct_savvy is important
import datastruct_savvy

from rdflib import Graph
from turtle_utils import generate_uri as gu
from turtle_grapher import generate_output, generate_graph, generate_turtle_skeleton

from os.path import basename

# bruteforce
from insert import upload_graph


def call_ectyper(graph, args_dict):
    # i don't intend to import anything from ECTyper (there are a lot of
    # imports in it - not sure if we'll use them all)
    import subprocess

    from rdflib import Literal
    from ast import literal_eval
    from os.path import splitext

    #logging.info('calling ectyper from fun call_ectyper')
    # concurrency is handled at the batch level, not here (note: this might change)
    # we only use ectyper for serotyping and vf, amr is handled by rgi directly
    if not args_dict['disable_serotype'] or not args_dict['disable_vf']:
        ectyper_dict = subprocess.check_output(['./ecoli_serotyping/src/Tools_Controller/tools_controller.py',
                                                '-in', args_dict['i'],
                                                '-s', str(
                                                    int(not args_dict['disable_serotype'])),
                                                '-vf', str(
                                                    int(not args_dict['disable_vf']))
                                                ])
        #logging.info('inner call completed')

        # because we are using check_output, this catches any print messages from tools_controller
        # TODO: switch to pipes
        print ectyper_dict.lower()
        if 'error' in ectyper_dict.lower():
            #logging.error('ectyper failed for' + args_dict['i'])
            print 'ECTyper failed for: ', args_dict['i']
            print 'returning graph w/o serotype'
            return graph

        #logging.info('evalulating ectyper output')
        # generating the dict
        ectyper_dict = literal_eval(ectyper_dict)
        # logging.info(ectyper_dict)
        #logging.info('evaluation okay')

        # TODO: edit ectyper sure were not using this ducktape approach
        # we are calling tools_controller on only one file, so grab that dict
        key, ectyper_dict = ectyper_dict.popitem()

        if not args_dict['disable_serotype']:
            # serotype parsing
            #logging.info('parsing Serotype')
            graph = datastruct_savvy.parse_serotype(
                graph, ectyper_dict['Serotype'], args_dict['uriIsolate'])
            #logging.info('serotype parsed okay')

        if not args_dict['disable_vf']:
            # vf
            #logging.info('parsing vf')
            graph = datastruct_savvy.parse_gene_dict(
                graph, ectyper_dict['Virulence Factors'], args_dict['uriGenome'])
            #logging.info('vf parsed okay')

    if not args_dict['disable_amr']:
        # amr
        #logging.info('generating amr')
        graph = generate_amr(graph, args_dict['uriGenome'], args_dict['i'])
        #logging.info('amr generation okay')

    return graph


def generate_amr(graph, uriGenome, fasta_file):
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

    amr_results.rename(
        columns={'ORF_ID': 'contig_id', 'Best_Hit_ARO': 'GENE_NAME'}, inplace=True)

    # sometimes there are spaces at the end of the contig id, also we remove
    # the additional occurance tag that RGI adds to contig ids
    amr_results['contig_id'] = amr_results['contig_id'].apply(
        lambda n: n.strip().rsplit('_', 1)[0])

    # note: you might be tempted to prefix a set_index('contig_id') but
    # remember, the same contig might have multiple genes
    amr_results = amr_results.to_dict(orient='index')

    # we have to manually check for contigs with multiple genes
    # TODO: write something less horrendously slow and memory consuming
    amr_dict = {}
    for i in amr_results.keys():
        contig_id = amr_results[i]['contig_id']
        if contig_id not in amr_dict.keys():
            amr_dict[contig_id] = []
        amr_dict[contig_id].append(dict((keys, amr_results[i][keys]) for keys in (
            'START', 'STOP', 'GENE_NAME', 'ORIENTATION', 'CUT_OFF', 'GENE_NAME')))
    # wipe the amr_results early
    amr_results = None

    graph = datastruct_savvy.parse_gene_dict(graph, amr_dict, uriGenome)

    return graph


def savvy(args_dict):
    '''
    Args:
        args_dict(dict): i prefer working with args in a dictionary, rather than a namespace is all

    Returns:
        (rdflib.Graph): a graph object with the VF/AMR/Serotype added to it via ECTyper/RGI
    '''
    from os import remove  # for batch cleanup
    # starting #logging
    '''
    logging.basicConfig(
        filename='outputs/' + __name__ +
        args_dict['i'].split('/')[-1] + '.log',
        level=logging.INFO
    )
    '''

    print("Importing FASTA from: " + args_dict['i'])
    #logging.info('importing from' + args_dict['i'])

    # setting up graph
    graph = generate_graph()

    #logging.info('generating barebones ttl from file')
    graph = generate_turtle_skeleton(
        graph, args_dict['i'], args_dict['uriIsolate'], args_dict['uriGenome'])
    #logging.info('barebones ttl generated')

    if not (args_dict['disable_serotype'] and args_dict['disable_vf'] and args_dict['disable_amr']):
        #logging.info('calling ectyper')
        graph = call_ectyper(graph, args_dict)
        #logging.info('ectyper call completed')

    # individual fasta logs are wiped on completion (or you'd have several
    # thousand of these)
    #remove('outputs/' + __name__ + args_dict['i'].split('/')[-1] + '.log')
    print upload_graph(graph)
    return graph

if __name__ == "__main__":
    import argparse
    import os  # for batch cleanup

    from ConfigParser import SafeConfigParser
    from turtle_utils import generate_hash

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
    # note: by in large, we expect uri to be given as just the unique string
    # value  (be it the hash or the integer) without any prefixes, the actual
    # rdflib.URIRef object will be generated in this script
    # this is mainly for batch computation
    parser.add_argument(
        "--uriGenome",
        help="Allows the specification of the Genome URI separately. Expect just the hash (not an actual uri).",
    )
    # This is both for batch computation and for future extensions where there
    # are multiple sequencings per isolate (Campy)
    parser.add_argument(
        "--uriIsolate",
        help="Allows the specification of the Isolate URI separately. Expect just the integer (not the full :spfyID)",
        type=int
    )
    args = parser.parse_args()
    # we make a dictionary from the cli-inputs and add are uris to it
    # mainly used for when a func needs a lot of the args
    args_dict = vars(args)

    # starting#logging
    '''
    logging.basicConfig(
        filename='outputs/' + __name__ +
        args_dict['i'].split('/')[-1] + '.log',
        level=logging.INFO
    )
    '''

    # check if a genome uri isn't set yet
    if args_dict['uriIsolate'] is None:
        # this is temporary, TODO: include a spqarql query to the db
        uriIsolate = gu(':spfy' + str(hash(args_dict['i'].split('/')[-1])))
    else:
        uriIsolate = gu(':spfy' + args_dict['uriIsolate'])

    # if the fasta_file hash was not precomputed (batch scripts should
    # precompute it), we compute that now
    if args_dict['uriGenome'] is None:
        uriGenome = gu(':' + generate_hash(args_dict['i']))
    else:
        uriGenome = gu(':' + args_dict['uriGenome'])

    args_dict['uriIsolate'] = uriIsolate
    args_dict['uriGenome'] = uriGenome

    savvy(args_dict)
