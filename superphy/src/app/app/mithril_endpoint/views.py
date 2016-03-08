#views.py
from flask import render_template, redirect, url_for, abort, flash, request,\
    current_app, make_response, Response, jsonify
from superphy.shared import sparql

from . import mithril

@mithril.after_request
def add_header(response):
    #Append after request the nessesary headers.
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response

@mithril.route('/query', methods = ['GET'])
def query():
    #General query to test if the requests are working.
    data = (sparql.get_x_triples(10))
    return jsonify(data)

@mithril.route('/meta', methods = ['GET', 'POST'])
def meta():
    results = (sparql.get_all_genome_metadata())
    return jsonify(results)

@mithril.route('/genomes', methods = ['GET', 'POST'])
def genomes():
    results = (sparql.get_all_genome_metadata())
    return jsonify(results)

@mithril.route('/genome/<genomeid>', methods = ['GET', 'POST'])
def genome(genomeid):
    results = (sparql.get_genome_metadata(genomeid))
    return jsonify(results)

@mithril.route('/genes', methods = ['GET', 'POST'])
def genes():
    """
    General query that returns all genes and their metadata.
    """
    results = (sparql.get_all_genes())
    return jsonify(results)