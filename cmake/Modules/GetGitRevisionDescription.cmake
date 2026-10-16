# - Returns a version string from Git
#
# These functions force a re-configure on each git commit so that you can
# trust the values of the variables in your build system.
#
#  get_git_head_revision(<refspecvar> <hashvar> [<additional arguments to git describe> ...])
#
# Returns the refspec and sha hash of the current head revision
#
#  git_describe(<var> [<additional arguments to git describe> ...])
#
# Returns the results of git describe on the source tree, and adjusting
# the output so that it tests false if an error occurs.
#
#  git_get_exact_tag(<var> [<additional arguments to git describe> ...])
#
# Returns the results of git describe --exact-match on the source tree,
# and adjusting the output so that it tests false if there was no exact
# matching tag.
#
# Requires CMake 2.6 or newer (uses the 'function' command)
#
# Original Author:
# 2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
# http://academic.cleardefinition.com
# Iowa State University HCI Graduate Program/VRAC
#
# Copyright Iowa State University 2009-2010.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

if(__get_git_revision_description)
	return()
endif()
set(__get_git_revision_description YES)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(get_git_head_revision _refspecvar _hashvar)
	set(GIT_PARENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
	set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
		set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
		get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
		if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
			# We have reached the root directory, we are not in git
			set(${_refspecvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
			set(${_hashvar} "GITDIR-NOTFOUND" PARENT_SCOPE)
			return()
		endif()
		set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	endwhile()
	# check if this is a submodule
	if(NOT IS_DIRECTORY ${GIT_DIR})
		file(READ ${GIT_DIR} submodule)
		string(REGEX REPLACE "gitdir: (.*)\n$" "\\1" GIT_DIR_RELATIVE ${submodule})
		get_filename_component(SUBMODULE_DIR ${GIT_DIR} PATH)
		get_filename_component(GIT_DIR ${SUBMODULE_DIR}/${GIT_DIR_RELATIVE} ABSOLUTE)
	endif()
	set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
	if(NOT EXISTS "${GIT_DATA}")
		file(MAKE_DIRECTORY "${GIT_DATA}")
	endif()

	if(NOT EXISTS "${GIT_DIR}/HEAD")
		return()
	endif()
	set(HEAD_FILE "${GIT_DATA}/HEAD")
	configure_file("${GIT_DIR}/HEAD" "${HEAD_FILE}" COPYONLY)

	configure_file("${_gitdescmoddir}/GetGitRevisionDescription.cmake.in"
		"${GIT_DATA}/grabRef.cmake"
		@ONLY)
	include("${GIT_DATA}/grabRef.cmake")

	set(${_refspecvar} "${HEAD_REF}" PARENT_SCOPE)
	set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
endfunction()

function(git_describe _var)
	if(NOT GIT_FOUND)
		find_package(Git QUIET)
	endif()
	get_git_head_revision(refspec hash)
	if(NOT GIT_FOUND)
		set(${_var} "GIT-NOTFOUND" PARENT_SCOPE)
		return()
	endif()
	if(NOT hash)
		set(${_var} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
		return()
	endif()

	# TODO sanitize
	#if((${ARGN}" MATCHES "&&") OR
	#	(ARGN MATCHES "||") OR
	#	(ARGN MATCHES "\\;"))
	#	message("Please report the following error to the project!")
	#	message(FATAL_ERROR "Looks like someone's doing something nefarious with git_describe! Passed arguments ${ARGN}")
	#endif()

	#message(STATUS "Arguments to execute_process: ${ARGN}")

	execute_process(COMMAND
		"${GIT_EXECUTABLE}"
		describe
		${hash}
		${ARGN}
		WORKING_DIRECTORY
		"${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE
		res
		OUTPUT_VARIABLE
		out
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT res EQUAL 0)
		set(out "${out}-${res}-NOTFOUND")
	endif()

	set(${_var} "${out}" PARENT_SCOPE)
endfunction()

function(git_get_exact_tag _var)
	git_describe(out --exact-match ${ARGN})
	set(${_var} "${out}" PARENT_SCOPE)
endfunction()

## An optional 4th argument reports whether every git call succeeded. The metadata
## strings cannot signal that on their own: a branch may be named like a sentinel.
function(git_short_info _refspecvar _hashvar _lc_datevar)
	if(NOT GIT_FOUND)
		find_package(Git QUIET)
	endif()
	if(ARGC GREATER 3)
		set(${ARGV3} FALSE PARENT_SCOPE)
	endif()
	get_git_head_revision(refspec hash)
	if(NOT GIT_FOUND)
		set(${_refspecvar} "GIT-NOTFOUND" PARENT_SCOPE)
		set(${_hashvar} "GIT-NOTFOUND" PARENT_SCOPE)
		set(${_lc_datevar} "GIT-NOTFOUND" PARENT_SCOPE)
		return()
	endif()
	if(NOT hash)
		set(${_refspecvar} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
		set(${_hashvar} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
		set(${_lc_datevar} "HEAD-HASH-NOTFOUND" PARENT_SCOPE)
		return()
	endif()

	set(_ok TRUE)

	execute_process(COMMAND
		"${GIT_EXECUTABLE}"
		rev-parse
		--short=7
		${hash}
		WORKING_DIRECTORY
		"${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE
		res
		OUTPUT_VARIABLE
		out
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT res EQUAL 0)
		set(out "${out}-${res}-NOTFOUND")
		set(_ok FALSE)
	endif()

	set(${_hashvar} "${out}" PARENT_SCOPE)

	execute_process(COMMAND
		"${GIT_EXECUTABLE}"
		rev-parse
		--abbrev-ref
		HEAD
		WORKING_DIRECTORY
		"${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE
		res
		OUTPUT_VARIABLE
		out
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT res EQUAL 0)
		set(out "${out}-${res}-NOTFOUND")
		set(_ok FALSE)
	endif()

	set(${_refspecvar} "${out}" PARENT_SCOPE)

	execute_process(COMMAND
		"${GIT_EXECUTABLE}"
		log
		-n 1
		--simplify-by-decoration
		--pretty=%ai
		WORKING_DIRECTORY
		"${CMAKE_SOURCE_DIR}"
		RESULT_VARIABLE
		res
		OUTPUT_VARIABLE
		out
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT res EQUAL 0)
		set(out "${out}-${res}-NOTFOUND")
		set(_ok FALSE)
	endif()

	set(${_lc_datevar} "${out}" PARENT_SCOPE)

	if(ARGC GREATER 3)
		set(${ARGV3} ${_ok} PARENT_SCOPE)
	endif()
endfunction()

## Reads what git archive filled into <_file>, the source tree's .git_archival.txt (export-subst
## in .gitattributes): the commit, its date and the nearest release tag (git describe; empty if
## there is none). <_foundvar> is TRUE only in a source tree made by git archive; in a checkout
## the file still holds its $Format placeholders. Hash and date come in the formats of
## git_short_info. src/tests/git_archival/run.cmake tests it.
function(git_archival_info _file _describevar _hashvar _lc_datevar _foundvar)
	set(${_foundvar} FALSE PARENT_SCOPE)
	if(NOT EXISTS "${_file}")
		return()
	endif()
	file(STRINGS "${_file}" lines)
	set(hash "")
	set(lc_date "")
	set(describe "")
	foreach(line IN LISTS lines)
		if(line MATCHES "^node: ([0-9a-f]+)$")
			set(hash "${CMAKE_MATCH_1}")
		elseif(line MATCHES "^node-date: (.+)$")
			set(lc_date "${CMAKE_MATCH_1}")
		elseif(line MATCHES "^describe-name: (.*)$")
			set(describe "${CMAKE_MATCH_1}")
		endif()
	endforeach()
	if(NOT hash)
		return()
	endif()
	## git before 2.35 knows no describe options and leaves the placeholder as it is
	if(describe MATCHES "^%\\(")
		message(WARNING "${_file}: the git that made this archive did not fill in describe-name "
		                "(that needs git 2.35 or later), so the build cannot tell a release.")
		set(describe "")
	endif()
	string(SUBSTRING "${hash}" 0 7 hash)
	## 2026-09-30T08:25:00+02:00 (%cI) -> 2026-09-30 08:25:00 +0200 (%ai); since git 2.45, %cI
	## writes UTC as Z
	string(REGEX REPLACE "Z$" "+00:00" lc_date "${lc_date}")
	string(REGEX REPLACE "^([0-9-]+)T([0-9:]+)([+-][0-9][0-9]):?([0-9][0-9])$" "\\1 \\2 \\3\\4" lc_date "${lc_date}")
	set(${_describevar} "${describe}" PARENT_SCOPE)
	set(${_hashvar} "${hash}" PARENT_SCOPE)
	set(${_lc_datevar} "${lc_date}" PARENT_SCOPE)
	set(${_foundvar} TRUE PARENT_SCOPE)
endfunction()
