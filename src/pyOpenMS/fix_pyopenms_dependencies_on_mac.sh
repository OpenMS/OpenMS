#!/bin/bash

# pyopenms.so on mac contains absolute pathes to referenced libraries libOpenMS et al.
#
# you can see this using
#    $ otool -L pyopenms.so
#
# here we look at all of them and fix them using install_name_tool as seen below
#

for TOFIX in $(find . -name pyopenms.so); do
    echo
    echo "fix $TOFIX now:"
    echo
    for LIB in libOpenMS libOpenSwathAlgo libSuperHirn; do
        # find absolute path
        REF=$(otool -L $TOFIX | grep $LIB.dylib | cut -d " " -f 1)
        echo "    found link $REF"
        install_name_tool -change $REF @loader_path/$LIB.dylib $TOFIX
    done;
    echo
    echo "    now otool -L says:"
    echo
    otool -L $TOFIX | while read -r line ; do
        echo "    $line"
    done

done

for TOFIX in $(find . -name libOpenMS.dylib); do
    echo
    echo "fix $TOFIX now:"
    echo
    for LIB in libOpenSwathAlgo.dylib QtNetwork QtCore libz.1.dylib; do
        # find absolute path
        REF=$(otool -L $TOFIX | grep $LIB | cut -d " " -f 1)
        echo "    found link $REF"
        install_name_tool -change $REF @loader_path/$LIB $TOFIX
    done;
    echo
    echo "    now otool -L says:"
    echo
    otool -L $TOFIX | while read -r line ; do
        echo "    $line"
        done
done


for TOFIX in $(find . -name libSuperHirn.dylib); do
    echo
    echo "fix $TOFIX now:"
    echo
    for LIB in libOpenMS.dylib libOpenSwathAlgo.dylib QtCore QtNetwork libz.1.dylib; do
        # find absolute path
        REF=$(otool -L $TOFIX | grep $LIB | cut -d " " -f 1)
        echo "    found link $REF"
        install_name_tool -change $REF @loader_path/$LIB $TOFIX
    done;
    echo
    echo "    now otool -L says:"
    echo
    otool -L $TOFIX | while read -r line ; do
        echo "    $line"
        done
done;

for TOFIX in $(find . -name QtNetwork); do
    echo
    echo "fix $TOFIX now:"
    echo
    for LIB in QtCore; do
        # find absolute path
        REF=$(otool -L $TOFIX | grep $LIB | cut -d " " -f 1)
        echo "    found link $REF"
        install_name_tool -change $REF @loader_path/$LIB $TOFIX
    done;
    echo
    echo "    now otool -L says:"
    echo
    otool -L $TOFIX | while read -r line ; do
        echo "    $line"
        done
done;
