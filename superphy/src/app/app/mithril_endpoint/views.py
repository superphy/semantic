#views.py
"""
views.py
provides the endpoints for this particular blueprint.
"""
from flask import jsonify
from superphy.shared import sparql

from . import mithril

@mithril.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response

@mithril.route('/genomes', methods=['GET', 'POST'])
def genomes():
    """
    General query that returns all genomes and their metadata.
    """
    results = (sparql.get_all_genome_metadata())
    return jsonify(results)

@mithril.route('/genome/<genomeid>', methods=['GET', 'POST'])
def genome(genomeid):
    """
    Returns the metadata of a particular genome in json format.
    """
    results = (sparql.get_genome_metadata(genomeid))
    return jsonify(results)

@mithril.route('/genes', methods=['GET', 'POST'])
def genes():
    """
    General query that returns all genes and their metadata.
    """
    results = (sparql.get_all_genes())
    return jsonify(results)

@mithril.route('/gene/<geneid>', methods=['GET', 'POST'])
def gene(geneid):
    """
    Returns the metadata of a particular gene in json format.
    """
    results = (sparql.get_gene(geneid))
    return jsonify(results)
