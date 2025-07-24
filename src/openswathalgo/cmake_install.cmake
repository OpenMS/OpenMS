# Install script for directory: /workspaces/OpenMS/src/openswathalgo

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/opt/rh/gcc-toolset-14/root/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "library" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so"
         RPATH "\$ORIGIN/../lib/")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/workspaces/OpenMS/lib/libOpenSwathAlgo.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so"
         OLD_RPATH ":::::::::::::::"
         NEW_RPATH "\$ORIGIN/../lib/")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/rh/gcc-toolset-14/root/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenSwathAlgo.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "library" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/ALGO" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/ALGO" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionHelper.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO/DATAACCESS" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/Transitions.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "OpenSwathAlgo_headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenMS/OPENSWATHALGO" TYPE FILE FILES "/workspaces/OpenMS/src/openswathalgo/include/OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h")
endif()

