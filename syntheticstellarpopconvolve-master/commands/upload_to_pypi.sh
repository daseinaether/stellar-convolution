#!/bin/bash
# NOTE: this script needs revision. Does not work properly currently

#
NAME_CURRENT_FILE="`realpath \"$0\"`"
DIRNAME_CURRENT_FILE=$(dirname $NAME_CURRENT_FILE)
DIRNAME_PROJECT_ROOT=$(dirname $DIRNAME_CURRENT_FILE)
cd $DIRNAME_PROJECT_ROOT

# Get current version
VERSION_NUMBER=$(grep -oP '__version__ = "\K[^"]+' syntheticstellarpopconvolve/_version.py | awk '{print $1}')

# Create dist
echo "Creating source distribution for syntheticstellarpopconvolve-$VERSION_NUMBER"
python setup.py sdist

# Checking validity
echo ""
echo "Checking validity of for syntheticstellarpopconvolve-$VERSION_NUMBER source distribution:"
twine check dist/syntheticstellarpopconvolve-$VERSION_NUMBER.tar.gz

echo ""
echo "Uploading syntheticstellarpopconvolve version $VERSION_NUMBER to pypi."
read -p "Continue? y/n " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    twine upload dist/syntheticstellarpopconvolve-$VERSION_NUMBER.tar.gz
fi
