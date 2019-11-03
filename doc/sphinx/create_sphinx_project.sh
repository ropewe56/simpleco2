#!/bin/bash

#sphinx-quickstart
# => project1
# => split source and build

project1=$1

# first invocation
echo "sphinx-build -b html $project1/source $project1/build"
sphinx-build -b html $project1/source $project1/build

#
cd $project1

# set SPHINXBUILD="python3 -msphinx" in makefile
make SPHINXBUILD="python3 -msphinx" html
