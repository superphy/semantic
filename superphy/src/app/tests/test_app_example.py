"""
    testing module for /example/*
"""

import unittest

from flask import json, url_for

from apptester import AppTester
from superphy.app import create_app, db

class ExampleTestCase(AppTester):
    """
    tests the example namespace
    """
    def test_foo(self):
        """
        @example.route('/foo')
        """
        test_result = json.loads('{"foo": "/example/bar"}')
        resp = self.client.get(
            url_for('example.foo', _external=False),
            data=None,
            headers=self.get_headers()
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

    def test_bar(self):
        """
        @example.route('/bar')
        """
        test_result = {'bar': 'bar'}
        resp = self.client.get(
            url_for('example.bar', _external=False),
            headers=self.get_headers(),
            data=None
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

    def test_apha(self):
        """
        @example.route('/alpha')

        This test is to show that url_for is calling the function, not the
        route. This means that when you use it, you need to use the correct
        one.
        """
        test_result = {'route':'/alpha', 'function':'beta'}
        self.assertEqual(url_for('example.beta', _external=False), '/example/alpha')
        resp = self.client.get(
            url_for('example.beta', _external=False),
            headers=self.get_headers(),
            data=None
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

    def test_bing(self):
        """
        @example.route('/bing')
        @api.route('/users', methods=['POST'])
        """
        
        resp = self.client.post(
            url_for('api.new_user', _external=False),
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(201, resp.status_code)

        test_result = {"YOU DID":"IT"}
        resp = self.client.get(
            url_for('example.test_for_login', _external=False),
            headers=self.get_api_headers('foo', 'bar'),
            data=None
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

        #test accessing non-protected data with a username
        resp = self.client.get(
            url_for('example.bar', _external=False),
            headers=self.get_api_headers('foo', 'bar'),
            data=None
        )

if __name__ == '__main__':
    unittest.main()
