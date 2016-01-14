#!/bin/bash
#This removes the old symlink from apache,
#and adds a new one; the new one points to the working 
#index.html that is in the folder that the script is run on.
source bash/config

folder=${1:-${default}}
if [ "$folder" == "$NULL" ]
	then
	read -r -p "What folder in apache? " folder
fi

echo symlink = ${folder} ...

origin="$(pwd)"
static=${origin}${HTML_DIRECTORY}
apachefolder=$APACHE$folder


[ -L $apachefolder ] && sudo rm $apachefolder && echo "Removing existing symlink at ${apachefolder} ..."
echo Making symlink from ${static}, to ${apachefolder} ...


sudo ln  $static $apachefolder -s

echo Apache Deployed!