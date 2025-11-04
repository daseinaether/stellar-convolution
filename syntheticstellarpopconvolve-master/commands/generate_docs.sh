# !/bin/bash
# Script to generate the docs

#
NAME_CURRENT_FILE="`realpath \"$0\"`"
DIRNAME_CURRENT_FILE=$(dirname $NAME_CURRENT_FILE)
DIRNAME_PROJECT_ROOT=$(dirname $DIRNAME_CURRENT_FILE)

#
DOCS_DIR="$DIRNAME_PROJECT_ROOT/docs/"


pandoc ../README.md -t rst -o ../docs/source/_includes/readme.rst

#
# echo "$NAME_CURRENT_FILE"
# echo "$DIRNAME_CURRENT_FILE"
# echo "$DIRNAME_PROJECT_ROOT"
# echo "$DOCS_DIR"

# Generate docs
echo "Generating documentation"
cd $DOCS_DIR
make clean
make html
echo "Done"
