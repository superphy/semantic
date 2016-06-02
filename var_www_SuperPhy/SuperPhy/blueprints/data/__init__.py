from flask import Blueprint

data = Blueprint('data', __name__)

from SuperPhy.blueprints.data import views

@data.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response
