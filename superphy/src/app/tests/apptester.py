import os
import unittest


from flask import json, jsonify, url_for

from superphy.app.config import basedir
from superphy.app import create_app, db
from superphy.app import User

class AppTester(unittest.TestCase):
    def setUp(self):
        self.app = create_app('testing')
        self.client = self.app.test_client()

    def tearDown(self):
        with self.app.app_context():
            db.session.remove()
            db.drop_all()