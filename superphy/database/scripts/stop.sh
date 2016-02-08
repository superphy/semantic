#!bin/bash
source scripts/config
result="$(ps ax | grep ${PORT} | grep java | awk '{print $1}') $(lsof -t -i:$PORT)" && echo $result
while [ "$result" != " " ]; do
	result="$(ps ax | grep ${PORT} | grep java | awk '{print $1}') $(lsof -t -i:$PORT)" && echo $result
	if [ "$result" != " " ]; then
		{
			kill -9 $(echo $result | awk '{print $1}')
		} || {
			echo "Sudo killing existing blazegraph..."
			sudo kill -9 $(echo $result | awk '{print $1}')
		}
		echo "Blazegraph has been terminated!"
	else
		echo "No Blazeraph running on port $PORT..."
	fi
done