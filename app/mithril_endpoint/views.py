#views.py
from flask import render_template, redirect, url_for, abort, flash, request,\
    current_app, make_response, Response, jsonify
from superphy import sparql

from . import mithril

#Which one?
@mithril.after_app_request
@mithril.after_request
def add_header(response):
    #Append after request the nessesary headers.
    response.headers['Access-Control-Allow-Origin'] = '*'
    return response

@mithril.route('/query', methods = ['GET'])
def query():
    #General query to test if the requests are working.
    data = (sparql.get_x_tripples(3))
    return jsonify(data)