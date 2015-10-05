#!/bin/bash
./start.db
python -m app/__init__.py
find . -name "*.pyc" -exec rm -rf {} \;