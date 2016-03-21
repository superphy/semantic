"""
__init__.py

app for the flask server

"""

import os

from flask import Flask

from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.httpauth import HTTPBasicAuth

from .config import config



db = SQLAlchemy()
auth = HTTPBasicAuth()


from .models import User

def create_app(config_name):
    """
    creates the application from the aggregation of blueprints
    """
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    db.init_app(app)

    with app.app_context():
        db.create_all()

    from.example import example as example_blueprint
    app.register_blueprint(example_blueprint, url_prefix='/example')

    from .data import data as data_blueprint
    app.register_blueprint(data_blueprint, url_prefix='/data')

    from .api import api as api_blueprint
    app.register_blueprint(api_blueprint, url_prefix='/api')

    from .upload import upload as upload_blueprint
    app.register_blueprint(upload_blueprint, url_prefix='/upload')

    return app
