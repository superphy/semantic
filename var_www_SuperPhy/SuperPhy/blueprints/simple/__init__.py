import os
from flask import request, jsonify, url_for, Blueprint

simple = Blueprint('simple', __name__)

@simple.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response

@simple.route('/foo')
def foo():
    return jsonify({'foo':url_for('simple.bar')})

@simple.route('/bar')
def bar():
    return jsonify({'bar':'bar'})
