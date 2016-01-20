#!bin/bash
#bash venv_linker.sh superphy/uploader/python_modules/ uploader

#target="superphy/front_end/python_modules"
#name="front_end"

target=${1:-"null"}
name=${2:-"null"}

num=$(grep -o "/" <<< ${symlink::-1} | wc -l)
for ((i=0;i<5;i++));
do
	target="../"$target
done
echo $target
ln $target $name -s