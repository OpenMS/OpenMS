#!/bin/sh
valgrind \
    --tool=memcheck\
    --leak-check=yes\
    --error-limit=no \
    --suppressions=valgrind-python.supp\
    --num-callers=10\
    -v\
    nosetests -w tests
