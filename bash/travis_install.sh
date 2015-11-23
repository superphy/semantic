#!/bin/bash

chmod a+x *.py

if ! find venv/bin/pip | read v; then
	echo Setting up Virtual Environment
	virtualenv --no-site-packages venv
fi

source venv/bin/activate
pip install --upgrade pip
pip install -r venv/requirements.txt

if ! find venv/bin/node | read v; then
    nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt
fi

deactivate