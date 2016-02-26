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
    data = (sparql.get_x_tripples(3))
    return jsonify(data)

@mithril.route('/meta', methods = ['GET', 'POST'])
def meta():
    results = (sparql.get_genome_meta_data())
    return jsonify(results)