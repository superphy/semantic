from flask import Blueprint

api = Blueprint('api', __name__)

from SuperPhy.blueprints.api import users
