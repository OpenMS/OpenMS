# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

### Helper for the installed-consumer tests of src/tests/external. A plain
### consumer project needs none of this: it links the imported targets and CMake
### resolves their include directories itself.
###
### openms_package_include_dirs(<target> <out_var>)
###
### The include directories an installed OpenMS package hands a consumer, as
### plain paths.
###
### Reading INTERFACE_INCLUDE_DIRECTORIES raw is not enough. The public headers
### are installed as a CMake header file set (cmake/install_macros.cmake), so the
### generated OpenMSTargets.cmake re-creates that file set on the imported target,
### and target_sources() appends the BASE_DIRS of a PUBLIC or INTERFACE header
### set to INTERFACE_INCLUDE_DIRECTORIES wrapped in $<BUILD_INTERFACE:...> (see
### the target_sources() documentation on file sets). Reading the property is
### therefore enough to see every header directory the package reports, named and
### default file sets alike, but the entries are generator expressions rather
### than paths.
###
### The wrapper evaluates to the directory in every context an imported target is
### used in -- it is dropped only where a project exports targets of its own --
### so consumers compile against the right headers either way. Only a
### configure-time reader of the property sees it, and a check that treats the raw
### entry as a path then quietly stops checking anything.
function(openms_package_include_dirs target out_var)
  get_target_property(_dirs ${target} INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT _dirs)
    set(_dirs "")
  endif()

  list(TRANSFORM _dirs REPLACE "^\\$<BUILD_INTERFACE:(.*)>$" "\\1")

  ## Anything the wrapper did not account for is an error rather than an entry
  ## silently skipped: being skipped is how the Eigen check in ../CMakeLists.txt
  ## stopped checking.
  if(_dirs MATCHES "\\$<")
    message(FATAL_ERROR
      "an include directory of ${target} is a generator expression that "
      "openms_package_include_dirs() cannot reduce to a path (${_dirs}). Teach it "
      "about that expression rather than leaving the checks built on it passing on "
      "a string that is not a directory.")
  endif()

  ## ".*" above rather than ".+", so that the empty
  ## "$<BUILD_INTERFACE:${EXTERNAL_INCLUDES}>" openms_add_library() writes for a
  ## library without external includes reduces to nothing and drops out here.
  list(FILTER _dirs EXCLUDE REGEX "^$")
  list(REMOVE_DUPLICATES _dirs)
  set(${out_var} "${_dirs}" PARENT_SCOPE)
endfunction()
