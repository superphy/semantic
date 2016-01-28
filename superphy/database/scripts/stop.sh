#!bin/bash
source scripts/config
kill -9 `lsof -t -i TCP:${PORT} -c java -a`

echo Blazegraph has been terminated!