#!/bin/bash
#Use this when you imported requirements and added another requirement to the environment.
source venv/bin/activate

pip freeze > "$(pwd)"/venv/requirements.txt
freeze > "$(pwd)"/venv/npm-requirements.txt
deactivate