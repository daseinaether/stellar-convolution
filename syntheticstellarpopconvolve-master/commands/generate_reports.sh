#!/bin/bash
# Script to generate the docstring coverage
# https://github.com/HunterMcGushion/docstr_coverage/tree/5092bb5898305c51a7a970af79e9b4d302b7ed1f

#
NAME_CURRENT_FILE="`realpath \"$0\"`"
DIRNAME_CURRENT_FILE=$(dirname $NAME_CURRENT_FILE)
DIRNAME_PROJECT_ROOT=$(dirname $DIRNAME_CURRENT_FILE)

#
TESTS_DIR="$DIRNAME_PROJECT_ROOT/syntheticstellarpopconvolve/tests"
REPORTS_DIR="$DIRNAME_PROJECT_ROOT/reports"
BADGE_DIR="$DIRNAME_PROJECT_ROOT/badges"


# echo "$NAME_CURRENT_FILE"
# echo "$DIRNAME_CURRENT_FILE"
# echo "$DIRNAME_PROJECT_ROOT"
# echo "$TESTS_DIR"
# echo "$REPORTS_DIR"
# echo "$BADGE_DIR"

# Create main reports directory
mkdir -p "$REPORTS_DIR"

## Docstring coverage:
command -v docstr-coverage >/dev/null 2>&1 || { echo >&2 "docstr-coverage is not installed.  Aborting."; exit 1; }

#
echo "Generating docstring report"
DOCSTRING_COV_DIR="$REPORTS_DIR/docstring_coverage"
mkdir -p "$DOCSTRING_COV_DIR/"
docstr-coverage syntheticstellarpopconvolve --exclude="$TESTS_DIR/*" -e ".*/tests" -e "$TESTS_DIR" -v 3 --badge "$DOCSTRING_COV_DIR/docstring_coverage.svg" > "$DOCSTRING_COV_DIR/docstring_coverage.txt" 2>&1
cp "$DOCSTRING_COV_DIR/docstring_coverage.svg" "$BADGE_DIR/docstring_coverage.svg"
echo "Done"

## test coverage
command -v coverage >/dev/null 2>&1 || { echo >&2 "coverage is not installed. Aborting."; exit 1; }
command -v coverage-badge >/dev/null 2>&1 || { echo >&2 "coverage-badge is not installed. Aborting."; exit 1; }

echo "Generating test coverage html report"
TEST_COV_DIR="$REPORTS_DIR/test_coverage"
mkdir -p "$TEST_COV_DIR/"
cd $TEST_COV_DIR
coverage run --source=syntheticstellarpopconvolve "$TESTS_DIR/main.py"
coverage html
coverage-badge > "$TEST_COV_DIR/test_coverage.svg"
cd $DIRNAME_PROJECT_ROOT
cp "$TEST_COV_DIR/test_coverage.svg" "$BADGE_DIR/test_coverage.svg"
echo "Done"
