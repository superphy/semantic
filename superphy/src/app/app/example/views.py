from flask import request, jsonify, url_for

from . import example

from .. import auth

@example.route('/foo')
def foo():
    return jsonify({'foo':url_for('example.bar')})

@example.route('/bar')
def bar():
    return jsonify({'bar':'bar'})

@example.route('/bing')
@auth.login_required
def test_for_login():
    """
    This is using the auth extenstion. Sign up a new user (/api), and test 
    it here.
    """
    return jsonify({"YOU DID":"IT"})

