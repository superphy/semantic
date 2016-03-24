#views.py
"""
views.py
provides the endpoints for this particular blueprint.
"""
from flask import jsonify
from superphy.shared import sparql

from . import data


@data.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response

@data.route('/query', methods=['GET'])
def query():
    """
    #General query to test if the requests are working.
    """
    data = (sparql.get_x_triples(10))
    return jsonify(data)

@data.route('/meta', methods=['GET', 'POST'])
def meta():
    """
    General query that returns all genomes and their metadata.
    """
    results = (sparql.get_all_genome_metadata())
    return jsonify(results)

@data.route('/genomes', methods=['GET', 'POST'])
def genomes():
    """
    General query that returns all genomes and their metadata.
    """
    results = (sparql.get_all_genome_metadata())
    return jsonify(results)

@data.route('/genome/<genomeid>', methods=['GET', 'POST'])
def genome(genomeid):
    """
    Returns the metadata of a particular genome in json format.
    """
    results = (sparql.get_genome_metadata(genomeid))
    return jsonify(results)

@data.route('/genes', methods=['GET', 'POST'])
def genes():
    """
    General query that returns all genes and their metadata.
    """
    results = (sparql.get_all_genes())
    return jsonify(results)

@data.route('/vf', methods=['GET', 'POST'])
def vfs():
    """
    General query that returns all virulence factors.
    """
    results = (sparql.get_all_genes('vf'))
    return jsonify(results)

@data.route('/amr', methods=['GET', 'POST'])
def amrs():
    """
    General query that returns all antimicrobial resistance genes.
    """
    results = (sparql.get_all_genes('amr'))
    return jsonify(results)


@data.route('/gene/<geneid>', methods=['GET', 'POST'])
def gene(geneid):
    """
    Returns the metadata of a particular gene in json format.
    """
    results = (sparql.get_gene(geneid))
    return jsonify(results)


@data.route('/region/<geneid>/<genomeid>', methods=['GET', 'POST'])
def region(geneid, genomeid):
    """
    Queries for the instances of geneid in genomeid.
    """
    results = (sparql.find_regions(geneid, genomeid))
    return jsonify(results)

