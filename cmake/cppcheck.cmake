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

## we use only source files for cppcheck
set(SOURCE_FILE_REGEX "\\.C$")

MACRO(ADD_CPPCHECK_TEST FILE_TO_TEST)
  string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${FILE_TO_TEST} )
  if(is_source_file)
    add_cppcheck_sources(${FILE_TO_TEST} ${FILE_TO_TEST} STYLE FAIL_ON_WARNINGS)
  endif(is_source_file)
ENDMACRO()

set(CPPCHECK_INCLUDEPATH_ARG ${OPENMS_INCLUDE_DIRS})

# library checks
set(OpenMS_cppcheck_sources)
foreach(i ${OpenMS_sources})
  add_cppcheck_test(${i})
endforeach()

# GUI library checks
foreach(i ${OpenMSVisual_sources})
  add_cppcheck_test(${i})
endforeach()

# TOPP checks
foreach(i ${TOPP_executables})
  add_cppcheck_test(${TOPP_DIR}/${i}.C)
endforeach()

# UTILS checks
foreach(i ${UTILS_executables})
  add_cppcheck_test(${UTILS_DIR}/${i}.C)
endforeach()

foreach(i ${GUI_executables})
  add_cppcheck_test(${GUI_DIR}/${i}.C)
endforeach()
