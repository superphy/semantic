import os
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.login import LoginManager
#from flask.ext.openid import OpenID
from config import basedir

app = Flask(__name__)
app.config.from_object('config')
db = SQLAlchemy(app)

from app import views, models