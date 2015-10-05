#!/bin/bash
./start.db
python -m app.__init__
find . -name "*.pyc" -exec rm -rf {} \;