"""
api for users
"""

from flask import abort, request, jsonify, g, url_for

from . import api
from ..models import User
from .. import auth, db, cors


@auth.verify_password
def verify_password(username_or_token, password):
    """
    takes a plain password as argument and returns True if the password is
    correct or False if not
    """
    # first try to authenticate by token
    user = User.verify_auth_token(username_or_token)
    if not user:
        # try to authenticate with username/password
        user = User.query.filter_by(username=username_or_token).first()
        if not user or not user.verify_password(password):
            return False
    g.user = user
    return True

@api.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type, authorization"
    return response

@api.route('/users', methods=['POST'])
def new_user():
    """
    a client can register a new user with a POST request to /api/users. The
    body of the request needs to be a JSON object that has username and
    password fields
    """
    username = request.json.get('username')
    password = request.json.get('password')
    if username is None or password is None:
        abort(400)    # missing arguments
    if User.query.filter_by(username=username).first() is not None:
        abort(400)    # existing user
    user = User(username=username)
    user.hash_password(password)
    db.session.add(user)
    db.session.commit()
    return (jsonify({'username': user.username}), 201,
            {'Location': url_for('api.get_user', id=user.id, _external=True)})

@api.route('/users/<int:id>')
def get_user(id):
    user = User.query.get(id)
    if not user:
        abort(400)
    return jsonify({'username': user.username})

@api.route('/token')
@auth.login_required
def get_auth_token():
    """
    If the token can be decoded then the id encoded in it is used to load the
    user, and that user is returned.
    """
    token = g.user.generate_auth_token(600)
    return jsonify({'token': token.decode('ascii'), 'duration': 600})

@api.route('/resource')
@auth.login_required
def get_resource():
    return jsonify({'data': 'Hello, %s!' % g.user.username})