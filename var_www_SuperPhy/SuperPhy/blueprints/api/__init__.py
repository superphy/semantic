from flask import Blueprint

api = Blueprint('api', __name__)

from SuperPhy.blueprints.api import users
from SuperPhy.blueprints.api import routes

@api.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type, authorization"
    return response