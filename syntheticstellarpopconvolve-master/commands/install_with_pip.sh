#!/bin/bash
# Script to install syntheticstellarpopconvolve in the current venv with editable flags. with this any change in this repo will be available in the venv

VERSION_NUMBER=$(grep -oP '__version__ = "\K[^"]+' ../syntheticstellarpopconvolve/_version.py | awk '{print $1}')
echo "installing syntheticstellarpopconvolve version $VERSION_NUMBER"

# we can only use python3 and python3, but allow
# the user to set these in environment variables
# PYTHON and PIP.
: "${PYTHON:="python3"}"
: "${PIP:="pip3"}"

# do stuff...
cd ..
$PYTHON setup.py clean
cd docs
$PIP uninstall -y syntheticstellarpopconvolve
cd ../
$PIP install -v  .
