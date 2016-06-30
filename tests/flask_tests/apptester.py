import os
import sys
import unittest

#NETWORK_FACEING_SUBDIRECTORY = 'app'
#TOPLEVEL = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
#sys.path.append(os.path.join(TOPLEVEL, NETWORK_FACEING_SUBDIRECTORY))

from flask import json, jsonify, url_for
from base64 import b64encode

import SuperPhy

from SuperPhy.config import basedir
from SuperPhy import db
from SuperPhy import User
from SuperPhy import create_app

class AppTester(unittest.TestCase):
    """
    The base testing functionality of flask apps is this superclass.
    Inherit this from another class to test REST endpoints.
    """
    def setUp(self):
        self.app = SuperPhy.create_app('testing')
        self.app.config['SERVER_NAME'] = 'localhost:5000'
        self.client = self.app.test_client()
        self.app_context = self.app.app_context()
        self.app_context.push()


    def tearDown(self):
        db.session.remove()
        db.drop_all()
        self.app_context.pop()

    def get_api_headers(self, username, password):
        return {
            'Authorization': 'Basic ' + b64encode(
                (username + ':' + password).encode('utf-8')).decode('utf-8'),
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        }

    def get_headers(self):
        return {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        }
