from flask import Flask

from flask.ext.sqlalchemy import SQLAlchemy

from config import config

db = SQLAlchemy()

def create_app(config_name):
    app = Flask(__name__)
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    db.init_app(app)

    #if not app.debug and not app.testing and not app.config['SSL_DISABLE']:
    #    from flask.ext.sslify import SSLify
    #    sslify = SSLify(app)
        
    from .data import data as data_blueprint
    app.register_blueprint(data_blueprint, url_prefix='/data')

    from .auth import auth as auth_blueprint
    app.register_blueprint(auth_blueprint, url_prefix='/auth')
     
    return app