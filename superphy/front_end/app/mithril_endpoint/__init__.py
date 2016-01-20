from flask import Blueprint

mithril = Blueprint('mithril', __name__)

from . import views