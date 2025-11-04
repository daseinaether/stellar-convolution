#!/bin/bash

# Script to install syntheticstellarpopconvolve in the current venv

VERSION_NUMBER=$(grep -oP '__version__ = "\K[^"]+' syntheticstellarpopconvolve/_version.py | awk '{print $1}')

echo "installing syntheticstellarpopconvolve version $VERSION_NUMBER"

# we can only use python3 and python3, but allow
# the user to set these in environment variables
# PYTHON and PIP.
: "${PYTHON:="python3"}"
: "${PIP:="pip3"}"

# do stuff...
$PYTHON setup.py clean
cd docs
$PIP uninstall -y syntheticstellarpopconvolve
cd ../
$PYTHON setup.py build --force
$PYTHON setup.py sdist
$PIP install -v dist/syntheticstellarpopconvolve-$VERSION_NUMBER.tar.gz
