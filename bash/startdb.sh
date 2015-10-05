#!/bin/bash
cd app/db
java -jar bigdata-bundled.jar &> /dev/null
cd ..
cd ..