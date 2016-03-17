import os
import unittest

from flask import json, jsonify, url_for

from apptester import AppTester

class BasicTestCase(AppTester):
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
        self.assertEqual(test_result, json.loads(resp.data))
        self.assertEqual(201, resp.status_code)

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

        #test if we can get data from behind protected endpoint
        resp = self.client.get(
            url_for('api.get_resource', _external=False),
            headers=self.get_api_headers('foo', 'bar')
        )
        self.assertEqual(200, resp.status_code)

if __name__ == '__main__':
    unittest.main()
