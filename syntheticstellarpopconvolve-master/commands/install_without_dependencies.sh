#!/bin/bash
# NOTE: this script needs revision. Does not work properly currently

# Script to install syntheticstellarpopconvolve in the current venv

VERSION_NUMBER=$(grep -oP '__version__ = "\K[^"]+' ../syntheticstellarpopconvolve/_version.py | awk '{print $1}')
echo "installing syntheticstellarpopconvolve version $VERSION_NUMBER"

# Clean up all the stuff from before
python setup.py clean --all

# Go into a directory that doesnt contain 'syntheticstellarpopconvolve' so pip will uninstall the one in the venv, not the local one.
cd src
pip uninstall -y syntheticstellarpopconvolve
cd ../

# Create build, sdist and install it into the venv
python setup.py build --force
python setup.py sdist
pip install --ignore-installed --no-dependencies -v dist/syntheticstellarpopconvolve-$VERSION_NUMBER.tar.gz
