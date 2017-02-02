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

def call_ectyper(graph, args_dict):
    # i don't intend to import anything from ECTyper (there are a lot of
    # imports in it - not sure if we'll use them all)
    import subprocess

    from rdflib import Literal
    from ast import literal_eval
    from os.path import splitext

    fasta_file = args_dict['i']
    uriIsolate = args_dict['uriIsolate']

    logging.info('calling ectyper from fun call_ectyper')
    # concurrency is handled at the batch level, not here (note: this might change)
    # we only use ectyper for serotyping and vf, amr is handled by rgi directly
    ectyper_dict = subprocess.check_output(['./ecoli_serotyping/src/Tools_Controller/tools_controller.py',
                                            '-in', fasta_file,
                                            '-s', str(
                                                int(not args_dict['disable_serotype'])),
                                            '-vf', str(
                                                int(not args_dict['disable_vf']))
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
    ectyper_dict = ectyper_dict[splitext(args_dict['i'])[0].split('/')[-1]]

    if not args_dict['disable_serotype']:
        # serotype parsing
        logging.info('parsing Serotype')
        graph = datastruct_savvy.parse_serotype(graph, ectyper_dict['Serotype'], uriIsolate)
        logging.info('serotype parsed okay')

    if not args_dict['disable_vf']:
        # vf
        logging.info('parsing vf')
        graph = datastruct_savvy.parse_gene_dict(
            graph, ectyper_dict['Virulence Factors'], uriGenome)
        logging.info('vf parsed okay')

    if not args_dict['disable_amr']:
        # amr
        logging.info('generating amr')
        graph = generate_amr(graph, uriGenome, fasta_file)
        logging.info('amr generation okay')

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

def savvy(args):
    '''
    Args:
        args(dict): i prefer working with args in a dictionary, rather than a namespace is all

    Returns:
        (rdflib.Graph): a graph object with the VF/AMR/Serotype added to it via ECTyper/RGI
    '''

    # starting logging
    logging.basicConfig(
        filename='outputs/' + __name__ + args['i'].split('/')[-1] + '.log',
        level=logging.INFO
    )

    print("Importing FASTA from: " + args['i'])
    logging.info('importing from' + args['i'])

    # setting up graph
    graph = generate_graph()

    logging.info('generating barebones ttl from file')
    graph = generate_turtle_skeleton(graph, args['i'], uriIsolate, uriGenome)
    logging.info('barebones ttl generated')

    if not (args['disable_serotype'] and args['disable_vf'] and args['disable_amr']):
        logging.info('calling ectyper')
        graph = call_ectyper(graph, args_dict)
        logging.info('ectyper call completed')

    # individual fasta logs are wiped on completion (or you'd have several
    # thousand of these)
    os.remove('outputs/' + __name__ + args.i.split('/')[-1] + '.log')
    return graph
