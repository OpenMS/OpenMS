# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

include(CppcheckTargets)

if(CPPCHECK_FOUND)
  set(SOURCE_FILE_REGEX "\\.C$")

  # certain files need to be excluded, since they cause internal errors
  set(EXCLUDE_FILE_REGEX "ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.C$|FILTERING/DATAREDUCTION/DataFilters.C$")
  set(CPPCHECK_INCLUDEPATH_ARG ${OPENMS_INCLUDE_DIRS})

  # library checks
  set(OpenMS_cppcheck_sources)
  foreach(i ${OpenMS_sources})
    string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${i} )
    string( REGEX MATCH ${EXCLUDE_FILE_REGEX} is_excluded_file ${i} )
    if(is_source_file AND NOT is_excluded_file)
      #add_cppcheck_sources(${i} ${i} STYLE FAIL_ON_WARNINGS)
      list(APPEND OpenMS_cppcheck_sources ${i})
    elseif(is_source_file AND is_excluded_file)
      message(STATUS "Excluded ${i} from cppcheck tests")
    endif()
  endforeach()

  add_cppcheck_sources("OpenMS" ${OpenMS_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # GUI library checks
  set(OpenMSVisual_cppcheck_sources)
  foreach(i ${OpenMSVisual_sources})
    string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${i} )
    if(is_source_file)
      list(APPEND OpenMSVisual_cppcheck_sources ${i})
      # add_cppcheck_sources(${i} ${i} STYLE FAIL_ON_WARNINGS)
    endif()
  endforeach()

  add_cppcheck_sources("OpenMS_GUI" ${OpenMSVisual_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # TOPP checks
  set(TOPP_cppcheck_sources)
  foreach(i ${TOPP_executables})
    list(APPEND TOPP_cppcheck_sources ${TOPP_DIR}/${i}.C)
  endforeach()

  add_cppcheck_sources("TOPP" ${TOPP_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # UTILS checks
  set(UTILS_cppcheck_sources)
  foreach(i ${UTILS_executables})
    list(APPEND UTILS_cppcheck_sources ${UTILS_DIR}/${i}.C)
  endforeach()

  add_cppcheck_sources("UTILS" ${UTILS_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

  # GUI_TOOLS checks
  set(GUI_TOOLS_cppcheck_sources)
  foreach(i ${GUI_executables})
    list(APPEND GUI_TOOLS_cppcheck_sources ${GUI_DIR}/${i}.C)
  endforeach()

  add_cppcheck_sources("GUI_TOOLS" ${GUI_TOOLS_cppcheck_sources} STYLE FAIL_ON_WARNINGS)

else()
  message(STATUS "Missing CPPCHECK executable .. Abort CppCheck")
endif()
