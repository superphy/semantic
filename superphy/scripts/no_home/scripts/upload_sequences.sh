#!/bin/bash
#makes a superphy namespace
source venv/bin/activate
python "$(pwd)"/python_package/uploader/sequence_upload.py
deactivate

exit 0