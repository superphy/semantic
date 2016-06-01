
import os

from flask import Blueprint, request, jsonify, url_for
from werkzeug import secure_filename

from SuperPhy import auth

example = Blueprint('example', __name__)

@example.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response

@example.route('/foo')
def foo():
    return jsonify({'foo':url_for('example.bar')})

@example.route('/bar')
def bar():
    return jsonify({'bar':'bar'})

@example.route('/alpha')
def beta():
    return jsonify({'route':'/alpha', 'function':'beta'})

@example.route('/bing')
@auth.login_required
def test_for_login():
    """
    This is using the auth extenstion. Sign up a new user (/api), and test 
    it here.
    """
    return jsonify({"YOU DID":"IT"})


@example.route('/echo', methods=['POST'])
def echo():
    return jsonify(request.json)

@example.route('/file', methods=['POST'])
def file():
    file_ = request.files['file']
    if file_:

        filename = secure_filename(file_.filename)
        file.save(os.path.join(os.path.dirname(os.path.realpath(__file__)), filename))

