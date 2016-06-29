from flask import Blueprint

upload = Blueprint('upload', __name__)

from SuperPhy.blueprints.upload import old_views