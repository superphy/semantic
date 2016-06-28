
import os
from flask import Flask, render_template, jsonify

from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.httpauth import HTTPBasicAuth
from flask.ext.cors import CORS

db = SQLAlchemy()
auth = HTTPBasicAuth()
cors = CORS

app = Flask(__name__)

from SuperPhy.config import config

config_name = os.getenv('FLASK_CONFIG') or 'default'
app.config.from_object(config[config_name])
config[config_name].init_app(app)

cors(app)

db.init_app(app)
with app.app_context():
    db.create_all()


from SuperPhy.models import User

from SuperPhy.blueprints.data import data as data_blueprint
from SuperPhy.blueprints.api import api as api_blueprint
from SuperPhy.blueprints.upload import upload as upload_blueprint

@app.route("/")
def index():
	return render_template('index.html')

#@app.errorhandler(Exception)
#def exception_handler(error):
#    return jsonify({"ERROR": repr(error), "FOO": error})

app.register_blueprint(data_blueprint, url_prefix='/data')
app.register_blueprint(api_blueprint, url_prefix='/api')
app.register_blueprint(upload_blueprint, url_prefix='/upload')

def create_app(config_name):
    """
    #creates the application from the aggregation of blueprints
    """
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)
    cors(app)
    db.init_app(app)

    with app.app_context():
        db.create_all()

    app.register_blueprint(data_blueprint, url_prefix='/data')
    app.register_blueprint(api_blueprint, url_prefix='/api')
    app.register_blueprint(upload_blueprint, url_prefix='/upload')

    return app

if __name__ == "__main__":
    app.run()
