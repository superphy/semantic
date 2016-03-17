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
            '/example/foo',
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
            '/example/bar',
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
        test_result = {"YOU DID":"IT"}
        resp = self.client.post(
            '/api/users',
            headers=self.get_headers(),
            data='{"username":"foo", "password":"bar"}'
        )
        self.assertEqual(201, resp.status_code)

        resp = self.client.get(
            '/example/bing',
            headers=self.get_api_headers('foo', 'bar'),
            data=None
        )
        self.assertEqual(200, resp.status_code)
        self.assertEqual(test_result, json.loads(resp.data))

        #self.assertEqual(url_for('example.foo', _external=False), '/example/foo')


if __name__ == '__main__':
    unittest.main()
