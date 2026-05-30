# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
## export a single option indicating if boost static libs should be preferred
option(BOOST_USE_STATIC "Use Boost static libraries." ON)

#------------------------------------------------------------------------------
## export a single option indicating if Arrow/Parquet static libs should be preferred
option(ARROW_USE_STATIC "Use Arrow/Parquet static libraries." ON)

#------------------------------------------------------------------------------
## Wraps the common find boost code into a single call
## @param .. simply add all required components to the call
## @note This macro will define BOOST_MOC_ARGS that should be added to all moc
##       calls (see https://bugreports.qt-project.org/browse/QTBUG-22829)
macro(find_boost)
  set(Boost_USE_STATIC_LIBS ${BOOST_USE_STATIC})
  set(Boost_USE_MULTITHREADED  ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
  add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (boost tends to guess wrong lib names)
  set(Boost_COMPILER "")
  
  ## since boost 1.66 they add an architecture tag if you build with layout=versioned and since 1.69 even when you
  ## build with layout=tagged (which we do in the contrib)
  if(NOT Boost_ARCHITECTURE)
    set(Boost_ARCHITECTURE "-x64")
  endif()

  # help boost finding it's packages
  set(Boost_ADDITIONAL_VERSIONS
    "1.87.1" "1.87.0" "1.87"
    "1.86.1" "1.86.0" "1.86"
    "1.85.1" "1.85.0" "1.85"
    "1.84.1" "1.84.0" "1.84"
    "1.83.1" "1.83.0" "1.83"
    "1.82.1" "1.82.0" "1.82"
    "1.81.1" "1.81.0" "1.81")

  find_package(Boost 1.81.0 COMPONENTS ${ARGN} REQUIRED CONFIG)

endmacro(find_boost)

#------------------------------------------------------------------------------
## Checks if the user supplied package type is valid and aborts if not
## @param package_type The given package type
macro(is_valid_package package_type)
  list(FIND VALID_PACKAGE_TYPES ${package_type} list_pos)
  if( ${list_pos} EQUAL -1 )
  	message(STATUS "The PACKAGE_TYPE ${package_type} is invalid")
  	message(STATUS "Valid PACKAGE_TYPEs are:")
  	foreach( _vpt ${VALID_PACKAGE_TYPES} )
  		message(STATUS " * ${_vpt}")
  	endforeach()
  	message(FATAL_ERROR "Aborting ...")
  endif()
endmacro()
