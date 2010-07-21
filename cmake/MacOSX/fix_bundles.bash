#!/bin/bash
## 
## 04/16/2009 - Stephan Aiche
##  Fixes association of TOPPView and INIFileEditor to the Qt libraries
##  to point to the installed lib/ path in /Applications/OpenMS-xx/lib/
##
## Commandline arguments are: 
##  $1 -> install_name_tool
##  $2 -> QT_LIB_PATH
##  $3 -> fullpath to executables
##  $4 -> lib-path

for lib in QtOpenGL QtCore QtGui QtXml QtSql QtNetwork QtTest
do
    # change association to Qt libs
    $1 -change $2/$lib.framework/Versions/4/$lib \
        @executable_path/../../../lib/$lib.framework/Versions/4/$lib \
        $3/TOPPView.app/Contents/MacOS/TOPPView.exe

    $1 -change $2/$lib.framework/Versions/4/$lib \
        @executable_path/../../../lib/$lib.framework/Versions/4/$lib \
        $3/INIFileEditor.app/Contents/MacOS/INIFileEditor.exe

		$1 -change $2/$lib.framework/Versions/4/$lib \
				@executable_path/../../../lib/$lib.framework/Versions/4/$lib \
				$3/TOPPAS.app/Contents/MacOS/TOPPAS.exe
done

# fix openms 
for app in TOPPView INIFileEditor TOPPAS
do
		$1 -change $4/libOpenMS.dylib \
				@executable_path/../../../lib/libOpenMS.dylib \
				$3/$app.app/Contents/MacOS/$app.exe
done
