import os
import unittest

from flask import json, jsonify, url_for

from apptester import AppTester
from superphy.app import create_app, db
from superphy.app.models import User

class UserTestCase(AppTester):
    """
    tests the user namespace
    """
    def test_add_user(self):
        """
        @api.route('/users', methods=['POST'])
        """
        test_result = json.loads('{"username": "foo"}')

        #add a user to fresh database
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(201, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))


        #add a user that already exists
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(400, resp.status_code)

        #add another user to database
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"bing", "password":"bang"}'
        )
        self.assertEqual(201, resp.status_code)

    def test_get_user(self):
        """
        test if we can get data from behind protected endpoint
        """
        #without credentials
        test_result = json.loads('{"username": "foo"}')
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_headers()
        )
        self.assertEqual(401, resp.status_code)

        #without proper username & password
        test_result = json.loads('{"username": "foo"}')
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_api_headers('bing', 'bang')
        )
        self.assertEqual(401, resp.status_code)

        #add a user to fresh database
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(201, resp.status_code)

        #test if we can get data from behind protected endpoint
        test_result = json.loads('{"data": "Hello, foo!"}')
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_api_headers('foo', 'bar')
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

    def test_tokens(self):
        #init
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(201, resp.status_code)

        resp = self.client.get(
            url_for('api.get_auth_token', _external=False),
            headers=self.get_api_headers('foo', 'bar'),
            data=None
        )
        self.assertEqual(200, resp.status_code)
        old_token = json.loads(resp.data)['token']

        resp = self.client.get(
            url_for('api.get_auth_token', _external=False),
            headers=self.get_api_headers('foo', 'bar'),
            data=None
        )
        self.assertEqual(200, resp.status_code)

        #Check if token works for @auth.loginrequired
        test_result = json.loads('{"data": "Hello, foo!"}')
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_api_headers(json.loads(resp.data)['token'], 'x')
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

        #Check if old_token works for @auth.loginrequired
        test_result = json.loads('{"data": "Hello, foo!"}')
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_api_headers(old_token, 'x')
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

if __name__ == '__main__':
    unittest.main()
