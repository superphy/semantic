#views.py
"""
views.py
provides the endpoints for this particular blueprint.
"""
import datetime
import os
from flask import jsonify, request
from flask_wtf import Form
from wtforms import StringField
from wtforms.validators import DataRequired

from werkzeug import secure_filename

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
    data_ = request.get_json()
    print "data", data_
    genomes_ = data_["genome"]
    genes_ = data_["genes"]

    ## Prepping dictionary to be returned
    genomeDict = {}
    for genome in genomes_:
        genomeDict[genome] = {}
        for gene in genes_:
            genomeDict[genome][gene] = 0
    
    results = sparql.get_regions(genomes_, genes_)
    bindings = results['results']['bindings']
    for binding in bindings:
        accession = binding['Genome']['value'].split("#")[1]
        gene_name = binding['Gene_Name']['value']
        try:
            genomeDict[accession][gene_name] += 1
        except KeyError:
            "Genome or gene doesn't exist in dictionary"

    return jsonify(genomeDict)

@data.route('/categories/<type>', methods=['GET', 'POST'])
def getcategories(type):
    """
    Endpoint for returning categories of genes.
    """
    categoryDict = {}
    results = (sparql.get_categories(type))
    bindings = results['results']['bindings']
    for binding in bindings:
        category = binding["Category"]
        subcategory = binding["Subcategory"]
    return jsonify(results)

# results = (sparql.get_all_genome_metadata())
#     bindings = results['results']['bindings'][:5]
#     rows = []
#     for binding in bindings:
#         row = {}
#         for item in results['head']['vars']:
#             try:
#                 row[item] = binding[item]['value']
#             except KeyError:
#                 row[item] = ''
#         rows.append(row)
#     return jsonify({'data':rows})

# Uploading stuff (actually belongs in uploader)

class MyForm(Form):
    """
    This belongs in models. This is dummy  function.
    """
    name = StringField('name', validators=[DataRequired()])
#Change this later to only allow '.faldo'?
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif']) 
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@data.route('/', methods=['GET', 'POST'])
def upload_file():
    """
    We also want to be changing the filename somehow, and keep track of who 
    is uploading what
    """
    if request.method == 'POST':
        file_ = request.files['file0']
        path = os.path.join(os.path.abspath(__file__)
            .split('superphy')[:2][:1].pop(), 'superphy/superphy/posted_files')
        print path
        if file and allowed_file(file_.filename):
            filename = secure_filename(file_.filename)
            file_.save(os.path.join(path, filename))
            return jsonify({})
    return jsonify({})
