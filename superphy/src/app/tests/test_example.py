import os
import unittest


from flask import json, jsonify, url_for

from superphy.app.config import basedir
from superphy.app import create_app, db
from superphy.app import User

class BasicTestCase(unittest.TestCase):
    def setUp(self):
        self.app = create_app('testing')
        self.client = self.app.test_client()

    def tearDown(self):
        with self.app.app_context():
            db.session.remove()
            db.drop_all()

    def test_json(self):
        """
        tests the example function

        This is designed for the first test. If it breaks, it means the entire
        directory structure is broken. You have been moving things around.
        """
        test_result = json.loads('{"foo": "/example/bar"}')

        url = '/example/foo'
        headers = {"Content-Type": "application/json"}
        data = None

        with self.app.app_context():
            resp = self.client.get(url, data=data, headers=headers)
        self.assertEqual(test_result, json.loads(resp.data))

    def test_add_user(self):
        """
        tests adding a user

        this should be moved to it's own class!
        """
        test_result = json.loads('{"username": "foo"}')

        url = '/api/users'
        headers = {"Content-Type": "application/json"}
        data1 = '{"username":"foo", "password":"bar"}'
        data2 = '{"username":"bing", "password":"bang"}'

        #add a user to fresh database
        with self.app.app_context():
            resp = self.client.post(url, data=data1, headers=headers)
        self.assertEqual(test_result, json.loads(resp.data))
        self.assertEqual(201, resp.status_code)

        #add a user that already exists
        with self.app.app_context():
            resp = self.client.post(url, data=data1, headers=headers)
        self.assertEqual(400, resp.status_code)

        #add another user to database
        with self.app.app_context():
            resp = self.client.post(url, data=data2, headers=headers)
        self.assertEqual(201, resp.status_code)

if __name__ == '__main__':
    unittest.main()
