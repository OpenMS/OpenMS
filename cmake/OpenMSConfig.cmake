### cmake OpenMS config file for external code

### Find the OpenMS includes and library
set(OPENMS_FOUND 1)
set(OPENMS_INCLUDE_DIRS "/share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/include/" "/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/include" "/raid/sturm/OpenMS/include" "/include")
set(OPENMS_LIBRARIES_DIRS "/share/opt/x86_64_sl4/OpenMS_contrib-current-gcc_4.1.1/lib/" "/raid/sturm/OpenMS/lib")
set(OPENMS_LIBRARIES "OpenMS;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtOpenGL.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtOpenGL.so;-lGLU -lGL;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtGui.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtGui.so;/usr/lib64/libpng.a;/usr/X11R6/lib64/libSM.a;/usr/X11R6/lib64/libICE.a;/usr/X11R6/lib64/libXi.a;/usr/X11R6/lib64/libXrender.a;/usr/X11R6/lib64/libXrandr.a;/usr/lib64/libfreetype.a;/usr/lib64/libfontconfig.a;/usr/X11R6/lib64/libXext.a;/usr/X11R6/lib64/libX11.a;/usr/lib64/libm.a;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtXml.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtXml.so;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtSql.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtSql.so;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtNetwork.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtNetwork.so;/usr/lib64/libssl.a;optimized;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtCore.so;debug;/share/opt/x86_64_sl4/qt-4.4.3-gcc_4.1.1/lib/libQtCore.so;/usr/lib64/libz.a;/usr/lib64/librt.a;-lpthread;-ldl")
set(OpenMS_BUILD_SETTINGS_FILE "/raid/sturm/OpenMS/OpenMSBuildSettings.cmake")
set(OpenMS_USE_FILE "/raid/sturm/OpenMS/OpenMSUse.cmake")

