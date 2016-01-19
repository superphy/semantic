#!/bin/bash

GULPFILE_LOCATION="/app/static/"
source ../venv/bin/activate
#Point symlink to index html file.
bash bash/deploy_apache.sh
root="$(pwd)"
index="${root}""${GULPFILE_LOCATION}"
echo $index
#Compile Coffeescript
cd $index
gulp
cd $root

#Start Flask server
./manage.py runserver --host 0.0.0.0