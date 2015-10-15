#!/bin/bash
cd db
(java -jar bigdata-bundled.jar &> /dev/null) & sleep 0
cd ..