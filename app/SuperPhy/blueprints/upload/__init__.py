from flask import Blueprint

upload = Blueprint('upload', __name__)

from SuperPhy.blueprints.upload import old_views
from SuperPhy.blueprints.upload import views

'''
@upload.after_request
def add_header(response):
    """
    Append after request the nessesary headers.
    """
    response.headers['Access-Control-Allow-Headers'] = "accept, content-type"
    return response
'''
