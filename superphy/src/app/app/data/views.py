#views.py
"""
views.py
provides the endpoints for this particular blueprint.
"""
import datetime
from flask import jsonify, request

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
    results = (sparql.get_x_triples(10))
    return jsonify(results)

@data.route('/meta', methods=['GET', 'POST'])
def meta():
    """
    General query that returns all genomes and their metadata.
    """
    results = (sparql.get_all_genome_metadata())
    results['date'] = (datetime.datetime.now() + datetime.timedelta(minutes=30)).isoformat()
    return jsonify(results)

@data.route('/meta2', methods=['GET'])
def meta2():
    """
    General query that returns all genomes and metadata in a nicer format.
    """
    results = (sparql.get_all_genome_metadata())
    bindings = results['results']['bindings'][:5]
    rows = []
    for binding in bindings:
        row = {}
        for item in results['head']['vars']:
            try:
                row[item] = binding[item]['value']
            except KeyError:
                row[item] = ''
        rows.append(row)
    return jsonify({'data':rows})

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

@data.route('/regions/<genomeid>', methods=['GET', 'POST'])
def regions(genomeid):
    """
    Returns all the genes inside a particular genome
    """
    results = (sparql.find_alleles(genomeid))
    return jsonify(results)


@data.route('/region/<geneid>/<genomeid>', methods=['GET', 'POST'])
def region(geneid, genomeid):
    """
    Queries for the instances of geneid in genomeid.
    """
    results = (sparql.find_regions(geneid, genomeid))
    return jsonify(results)

@data.route('/genesearchresults', methods=['POST'])
def genesearchresults():
    """
    Endpoint for returning gene search results
    """
    data = request.get_json()
    print data
    genome = data["genome"]
    genes = data["genes"]
    results = sparql.get_regions(genome, genes)
    return jsonify(results)

