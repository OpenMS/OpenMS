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
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg, Stephan Aiche $
# --------------------------------------------------------------------------

set(DO_NOT_TEST_THESE_FILES_REGEX "/(moc|ui)_")
set(IGNORE_FILES_IN_BUILD_DIRECTORY "^${PROJECT_BINARY_DIR}")

MACRO(ADD_CPPLINT_TEST FILE_TO_TEST)
  string( REGEX MATCH ${DO_NOT_TEST_THESE_FILES_REGEX} do_not_test ${FILE_TO_TEST} )
  string( REGEX MATCH ${IGNORE_FILES_IN_BUILD_DIRECTORY} is_in_bin_dir ${FILE_TO_TEST})
  if(NOT do_not_test AND NOT is_in_bin_dir)
    add_test(${FILE_TO_TEST}_cpplint_test
      "${PYTHON_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/cmake/cpplint.py"
      "--verbose=5"
      "${PROJECT_SOURCE_DIR}/${FILE_TO_TEST}")

    set_tests_properties(
      ${FILE_TO_TEST}_cpplint_test
      PROPERTIES
      FAIL_REGULAR_EXPRESSION
      "${CPPLINT_FAIL_REGULAR_EXPRESSION}")
  endif()
ENDMACRO()

## create tests for all files in the individual file groups
foreach(i ${OpenMS_sources})
  add_cpplint_test(${i})
endforeach()

foreach(i ${OpenMSVisual_sources})
  add_cpplint_test(${i})
endforeach()

foreach(i ${TOPP_executables})
  add_cpplint_test(${TOPP_DIR}/${i}.C)
endforeach()

foreach(i ${UTILS_executables})
  add_cpplint_test(${UTILS_DIR}/${i}.C)
endforeach()

foreach(i ${GUI_executables})
  add_cpplint_test(${GUI_DIR}/${i}.C)
endforeach()

