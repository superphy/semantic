#!/usr/bin/python
import os
import sys
import logging

activate_this = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'venv/bin/activate_this.py')
execfile(activate_this, dict(__file__=activate_this))

logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
sys.path.insert(0,"/var/www/SuperPhy/")
#sys.path.insert(0,"/var/www/html/SuperPhy/")

from SuperPhy import app as application
application.secret_key = 'Add your secret key'
