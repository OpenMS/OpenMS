#!/bin/bash

# pyopenms.so on mac contains absolute pathes to referenced libraries libOpenMS et al.
#
# you can see this using
#    $ otool -L pyopenms.so
#
# here we look at all of them and fix them using install_name_tool as seen below
#

for PYOPENMS in $(find . -name pyopenms.so); do
    echo
    echo "fix $PYOPENMS now:"
    echo
    for LIB in libOpenMS libOpenSwathAlgo libSuperHirn; do
        # find absolute path
        REF=$(otool -L $PYOPENMS | grep $LIB.dylib | cut -d " " -f 1)
        echo "    found link $REF"
        install_name_tool -change $REF @loader_path/$LIB.dylib $PYOPENMS
    done;
    echo
    echo "    now otool -L says:"
    echo
    otool -L $PYOPENMS | while read -r line ; do
        echo "    $line"
        done

done;
