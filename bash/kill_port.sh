#!/bin/bash
source bash/config
port=${1:-${NULL}}
if [ "$port" == ${TRIPLESTORE_PORT} ]
#script to kill all programs using port 9999 (Blazegraph in particular). Use with caution.
	then kill -9 `lsof -t -i TCP:${TRIPLESTORE_PORT} -c java -a`
elif [ "$port" == "5000" ]
#script to kill a python server on port 5000. Use with caution.
then
	var1=$(lsof -i TCP:5000 | grep "python" | cut -d ' ' -f 3);
	kill -9 $var1 &> /dev/null;
else
	echo "Error: Choose a hardcoded port in the script to kill"
	.
fi