#!/bin/bash

source="$(pwd)"/venv/bin/activate
xvfb-run -a nosetests && cd tests && xvfb-run -a jasmine-ci && cd ../
deactivate