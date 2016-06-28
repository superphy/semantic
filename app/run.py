#!/usr/bin/python
import os
import sys
import logging
import subprocess

activate_this = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'venv/bin/activate_this.py')
execfile(activate_this, dict(__file__=activate_this))

#activate = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'venv/bin/activate')
#os.system("bash %s" % activate)

#logging.basicConfig(stream=sys.stderr)
sys.path.insert(0,"/var/www/SuperPhy/")

from SuperPhy import app as application
application.secret_key = 'Add your secret key'


application.run(host='0.0.0.0', debug=True, use_reloader=False, port=5000)