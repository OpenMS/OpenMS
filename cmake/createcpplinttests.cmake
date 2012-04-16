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
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# build file list for cpplint
set(OpenMS_cpplint_sources)
foreach(i ${OpenMS_sources})
if("1.${CMAKE_VERSION}" VERSION_LESS "1.2.8.0")
# Older than CMake 2.8.0
add_test(${i}_cpplint_test
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
else()
# CMake 2.8.0 and newer
add_test(NAME
	${i}_cpplint_test
	COMMAND
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
endif()
set_tests_properties(${i}_cpplint_test
PROPERTIES
FAIL_REGULAR_EXPRESSION
"${CPPLINT_FAIL_REGULAR_EXPRESSION}")
endforeach()

foreach(i ${OpenMSVisual_sources})
if("1.${CMAKE_VERSION}" VERSION_LESS "1.2.8.0")
# Older than CMake 2.8.0
add_test(${i}_visual_cpplint_test
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
else()
# CMake 2.8.0 and newer
add_test(NAME
	${i}_visual_cpplint_test
	COMMAND
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
endif()
set_tests_properties(${i}_visual_cpplint_test
PROPERTIES
FAIL_REGULAR_EXPRESSION
"${CPPLINT_FAIL_REGULAR_EXPRESSION}")
endforeach()

foreach(i ${TOPP_executables})
if("1.${CMAKE_VERSION}" VERSION_LESS "1.2.8.0")
# Older than CMake 2.8.0
add_test(${i}_TOPP_cpplint_test
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
else()
# CMake 2.8.0 and newer
add_test(NAME
	${i}_TOPP_cpplint_test
	COMMAND
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
endif()
set_tests_properties(${i}_TOPP_cpplint_test
PROPERTIES
FAIL_REGULAR_EXPRESSION
"${CPPLINT_FAIL_REGULAR_EXPRESSION}")
endforeach()

foreach(i ${UTILS_executables})
if("1.${CMAKE_VERSION}" VERSION_LESS "1.2.8.0")
# Older than CMake 2.8.0
add_test(${i}_UTILS_cpplint_test
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
else()
# CMake 2.8.0 and newer
add_test(NAME
	${i}_UTILS_cpplint_test
	COMMAND
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
endif()
set_tests_properties(${i}_UTILS_cpplint_test
PROPERTIES
FAIL_REGULAR_EXPRESSION
"${CPPLINT_FAIL_REGULAR_EXPRESSION}")
endforeach()
  
foreach(i ${GUI_executables})
if("1.${CMAKE_VERSION}" VERSION_LESS "1.2.8.0")
# Older than CMake 2.8.0
add_test(${i}_GUI_cpplint_test
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
else()
# CMake 2.8.0 and newer
add_test(NAME
	${i}_GUI_cpplint_test
	COMMAND
	"python"
	"cmake/cpplint.py"
	"--verbose=5"
	${i})
endif()
set_tests_properties(${i}_GUI_cpplint_test
PROPERTIES
FAIL_REGULAR_EXPRESSION
"${CPPLINT_FAIL_REGULAR_EXPRESSION}")
endforeach()

