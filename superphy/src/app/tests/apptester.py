import os
import unittest


from flask import json, jsonify, url_for
from base64 import b64encode

from superphy.app.config import basedir
from superphy.app import create_app, db
from superphy.app import User

class AppTester(unittest.TestCase):
    """
    The base testing functionality of flask apps is this superclass.
    Inherit this from another class to test REST endpoints.
    """
    def setUp(self):
        self.app = create_app('testing')
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