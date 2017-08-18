#!/usr/bin/env bash

# set the dictionary hash 'noise'
export PYTHONHASHSEED=0

# run nosetests
if [ $TEST_STYLE ] ; then
    nosetests -v tests/test_style.py ;
elif [ $TEST_FUNCTIONALITY ] ; then
    nosetests -v tests/test_umi_tools.py ;
elif [ $TEST_HELP ] ; then
    umi_tools whitelist --help;
    umi_tools extract --help;
    umi_tools group --help;
    umi_tools dedup --help;
    umi_tools count --help;
    umi_tools count_tab --help;
fi


