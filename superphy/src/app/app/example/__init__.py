from flask import Blueprint

example = Blueprint('example', __name__)

from . import views