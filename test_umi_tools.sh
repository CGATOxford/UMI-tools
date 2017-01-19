#!/usr/bin/env bash

# set the dictionary hash 'noise'
export PYTHONHASHSEED=0

# run nosetests
if [ $TEST_STYLE ] ; then
    nosetests -v tests/test_style.py ;
elif [ $TEST_FUNCTIONALITY ] ; then
    nosetests -v tests/test_umi_tools.py ;
fi


