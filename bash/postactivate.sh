#!/bin/bash

################################################################
## postactivate.sh
##
##   Scripts / code in this file will be run when
##   venv/bin/activate is called.
##  
##   This script is called at the end of the venv/bin/activate
##   bash script using "source". All changes will impact the
##   current shell and will remain after the venv is deactivated
##
##   Variable CWD will contain the full path to the bash/
##   directory.
##
##
## author = Matt Whiteside
## copyright = Copyright 2015, Public Health Agency of Canada
## license = CC-SY
## version = 4.0
## maintainer = Matt Whiteside
## email = matthew.whiteside@canada.ca
##
################################################################

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Initialize Superphy ENV variables
source "$CWD/init_env.sh"