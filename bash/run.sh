#!/bin/bash
source bash/config
source venv/bin/activate
cd "$(pwd)"${HTML_DIRECTORY}
gulp
cd ../..
./manage.py runserver --host 0.0.0.0