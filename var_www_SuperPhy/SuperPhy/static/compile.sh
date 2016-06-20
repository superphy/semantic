#!/bin/bash
#Runs gulp in app/static

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

source ../../../venv/bin/activate
gulp