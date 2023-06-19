#!/bin/bash

# find "dist" directory
DIST=$(find . -name "dist")
cd $DIST

#  whl inside
NAME=*.whl

# backup .original
cp $NAME $(basename $NAME .whl).original

unzip $NAME  
# move to pyopenms folder and relink the libaries
cd pyopenms

# move QtCore and QtNetwork from framework to pyopenms directory
# and remove the framework dirs
QTN=$(find . -name QtNetwork)
mv $QTN .
rm -rf ./QtNetwork.framework

QTC=$(find . -name QtCore)
mv $QTC .
rm -rf ./QtCore.framework

# call fix script (in pyopenms directory)
# relative path 
sh ../../fix_pyopenms_dependencies_on_mac.sh

cd ../ # back to .dist

# zip back to wheel 
zip -r $NAME pyopenms/ pyopenms*dist-info/  

echo "DONE" 



