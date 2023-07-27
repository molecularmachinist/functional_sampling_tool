#!/usr/bin/bash

set -e

# cmd params
noinstall=false
novenv=false
help=false
while test $# -gt 0
do
    case "$1" in
        --no-install) noinstall=true
            ;;
        --no-venv) novenv=true
            ;;
        -h) help=true
            ;;
        --help) help=true
            ;;
        --*) echo "Unrecognised option $1" ; exit 1
            ;;
    esac
    shift
done

if ${help}; then
    echo "Usage: run_tests.sh [-h|--help] [--no-install] [--no-venv]"
    echo "Params:"
    echo "-h|--help:    Print this help message and exit."
    echo "--no-install: Do not install fst or any dependencies before running tests."
    echo "              The old venv will not be cleaned (if the venv is even used)."
    echo "--no-venv:    Do not use the virtual env. No new venv will be made and the"
    echo "              current environment is used to run the tests and possibly install"
    echo "              the tool."
    exit
fi


if ${novenv}; then
    echo "Not runnning in venv" 2>&1 | tee fst_test_run/test_output.log
else
    if ${noinstall}; then
        source fst_test_run/venv/bin/activate
    else
        python -m venv --clear fst_test_run/venv 2>&1 | tee fst_test_run/test_output.log
        source fst_test_run/venv/bin/activate
    fi
fi

if ${noinstall}; then
    echo "Skipping reinstall of fst" 2>&1 | tee fst_test_run/test_output.log
else
    echo "----------------- Running pip install ----------------" 2>&1 | tee -a fst_test_run/test_output.log

    pip install coverage 2>&1 | tee -a fst_test_run/test_output.log
    pip install . 2>&1 | tee -a fst_test_run/test_output.log
fi


cd fst_test_run

echo fst is $(which fst) 2>&1 | tee -a test_output.log
python -c "import functional_sampling_tool as fst; print(fst.__file__)" 2>&1 | tee -a test_output.log

echo "----------------- Running tests -----------------------" 2>&1 | tee -a test_output.log
coverage run --source functional_sampling_tool ../tests/main.py -v 2>&1 | tee -a test_output.log

echo "----------------- Making coverage reports -------------" 2>&1 | tee -a test_output.log
coverage html 2>&1 | tee -a test_output.log
coverage report 2>&1 | tee -a test_output.log
echo "----------------- Done --------------------------------" 2>&1 | tee -a test_output.log
