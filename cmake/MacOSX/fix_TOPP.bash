#!/bin/bash
## 
## 04/16/2009 - Stephan Aiche
##  Fixes association of TOPP tools (excl. TOPPView and INIFileEditor) to the Qt libraries
##  to point to the installed lib/ path in /Applications/OpenMS-xx/lib/
##
## Commandline arguments are: 
##  $1 -> install_name_tool
##  $2 -> QT_LIBRARY_DIR
##  $3 -> path-to-tool/tool

for lib in QtOpenGL QtCore QtGui QtXml QtSql QtNetwork QtTest QtSvg QtWebKit QtXmlPatterns
do
    # change association of Qt libs
    $1 -change $2/$lib.framework/Versions/4/$lib \
        @executable_path/../lib/$lib.framework/Versions/4/$lib \
        $3
done
