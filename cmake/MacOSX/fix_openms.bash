#!/bin/bash
## 
## 04/16/2009 - Stephan Aiche
##  Fixes the association of libOpenMS.dylib to the Qt libraries
##  to point to the installed lib/ path in /Applications/OpenMS-xx/lib/
## 
## Commandline arguments are: 
##  $1 = install_name_tool
##  $2 = QT_LIB_PATH
##  $3 = lib-path

# fix Qt references of OpenMS
for lib in QtOpenGL QtCore QtGui QtXml QtSql QtNetwork
do
    # change association to Qt libs
    $1 -change $2/$lib.framework/Versions/4/$lib \
        $lib.framework/Versions/4/$lib \
        $3/libOpenMS.dylib
done

