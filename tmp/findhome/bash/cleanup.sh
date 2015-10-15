#!/bin/bash
#Run this script from project root
(sleep 1 ; find . -name "*.pyc" -exec rm -rf {} \;) & sleep 0