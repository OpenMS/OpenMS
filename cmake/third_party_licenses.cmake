# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

## The packages bundle third-party libraries, and most of their licenses require the
## license text to accompany the binaries. openms_install_third_party_licenses() installs
## these texts under share/OpenMS/LICENSES, in the component 'share' that every
## installation type includes. cmake/package_general.cmake calls it; the top-level
## CMakeLists.txt includes this file in every configuration, so that CI parses it.
##
## - vcpkg/<port>.txt: the license text of every vcpkg port of the target triplet, which
##   vcpkg installs as share/<port>/copyright. The packages carry the code of each port:
##   as a library next to libOpenMS (install(RUNTIME_DEPENDENCY_SET) in
##   package_general.cmake), linked into it (the Windows triplets are static), or compiled
##   in (header-only ports).
## - Qt/: Qt is not a vcpkg port. The Windows and macOS packages bundle it; the DEB depends
##   on the distribution's Qt instead. OpenMS uses Qt under the LGPL version 3, which
##   requires its text, the text of the GPL version 3 that it supplements, and a note on
##   where to get the source code of the bundled Qt.
## - homebrew/<formula>/: on macOS, Qt and the libraries it needs (glib, ICU, freetype, ...)
##   come from Homebrew, and the package bundles them. Homebrew installs the license files
##   of a formula into its keg, and these are installed for the given formulae and every
##   formula they depend on.
##
## openms_install_vendored_licenses() installs the licenses of the third-party code that
## OpenMS carries in its own source tree and compiles into its libraries (vendored/<name>/).
## Every installation of libOpenMS contains that code, the Python wheels too, so the
## top-level CMakeLists.txt calls it in every configuration, not only for packages.
##
## openms_install_contrib_licenses() installs the licenses of the libraries a build takes
## from the OpenMS contrib (OPENMS_CONTRIB_LIBS), as the pyOpenMS wheels do (contrib/<name>/).
## No package manager accounts for these; the texts are those of the contrib's versions, in
## cmake/third_party_licenses/contrib. The top-level CMakeLists.txt calls it in every
## configuration too.

## Sets <formulae_var> to the Homebrew formulae whose kegs hold the given paths, and
## <prefix_var> to the prefix of the Homebrew installation they belong to.
function(openms_homebrew_formulae_of formulae_var prefix_var)
  set(_formulae)
  set(_prefixes)
  foreach(_path IN LISTS ARGN)
    if(NOT EXISTS "${_path}")
      continue()
    endif()
    ## <prefix>/opt/<formula> and the files linked into <prefix>/lib resolve to the keg,
    ## <prefix>/Cellar/<formula>/<version>.
    file(REAL_PATH "${_path}" _real)
    if(_real MATCHES "^(.*)/Cellar/([^/]+)/[^/]+(/|$)")
      list(APPEND _prefixes "${CMAKE_MATCH_1}")
      list(APPEND _formulae "${CMAKE_MATCH_2}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _formulae)
  list(REMOVE_DUPLICATES _prefixes)
  list(LENGTH _prefixes _count)
  if(_count GREATER 1)
    ## A Mac can have two Homebrew installations, /opt/homebrew and /usr/local.
    message(FATAL_ERROR "The paths ${ARGN} belong to more than one Homebrew installation: "
                        "${_prefixes}")
  endif()
  set(${formulae_var} "${_formulae}" PARENT_SCOPE)
  set(${prefix_var} "${_prefixes}" PARENT_SCOPE)
endfunction()

## openms_install_third_party_licenses([QT_VERSION <version>]
##                                     [HOMEBREW_PREFIX <prefix> HOMEBREW_FORMULAE <formula>...])
##   QT_VERSION         the version of the Qt the package bundles; empty if it bundles none
##   HOMEBREW_PREFIX    the prefix of the Homebrew installation the formulae come from
##   HOMEBREW_FORMULAE  the Homebrew formulae the package bundles libraries of; their
##                      dependencies are added
function(openms_install_third_party_licenses)
  cmake_parse_arguments(PARSE_ARGV 0 arg "" "QT_VERSION;HOMEBREW_PREFIX" "HOMEBREW_FORMULAE")
  set(_licenses_dir "${INSTALL_SHARE_DIR}/LICENSES")

  if(OPENMS_USE_VCPKG AND VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET)
    file(GLOB _copyrights LIST_DIRECTORIES false
         "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/*/copyright")
    foreach(_copyright IN LISTS _copyrights)
      get_filename_component(_port_dir "${_copyright}" DIRECTORY)
      get_filename_component(_port "${_port_dir}" NAME)
      install(FILES "${_copyright}"
              DESTINATION "${_licenses_dir}/vcpkg"
              RENAME "${_port}.txt"
              COMPONENT share)
    endforeach()
    list(LENGTH _copyrights _count)
    message(STATUS "Packaging the license texts of ${_count} vcpkg ports")
  endif()

  if(arg_QT_VERSION)
    set(_qt_licenses "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/third_party_licenses/Qt")
    ## The variables README.txt.in refers to.
    set(QT_LICENSE_VERSION "${arg_QT_VERSION}")
    string(REGEX MATCH "^[0-9]+\\.[0-9]+" QT_LICENSE_MAJOR_MINOR "${arg_QT_VERSION}")
    set(QT_LICENSE_HOMEBREW_NOTE "")
    if(arg_HOMEBREW_FORMULAE)
      string(CONCAT QT_LICENSE_HOMEBREW_NOTE
             "\nOn macOS, Qt comes from Homebrew. The Homebrew formulae, including any patches\n"
             "Homebrew applies to Qt, are at https://github.com/Homebrew/homebrew-core; ../homebrew\n"
             "holds the license files of the formulae whose libraries this package bundles.\n")
    endif()
    configure_file("${_qt_licenses}/README.txt.in"
                   "${PROJECT_BINARY_DIR}/third_party_licenses/Qt/README.txt" @ONLY)
    install(FILES "${_qt_licenses}/LGPL-3.0-only.txt"
                  "${_qt_licenses}/GPL-3.0-only.txt"
                  "${PROJECT_BINARY_DIR}/third_party_licenses/Qt/README.txt"
            DESTINATION "${_licenses_dir}/Qt"
            COMPONENT share)
    message(STATUS "Packaging the Qt ${arg_QT_VERSION} licenses")
  endif()

  if(arg_HOMEBREW_FORMULAE)
    ## The package bundles libraries of these formulae and their dependencies, so configuring
    ## fails rather than leave their license files out. The brew of the installation the
    ## formulae come from is used, not the one on the PATH: a brew answers for its own
    ## installation only, and a Mac can have two.
    set(_brew "${arg_HOMEBREW_PREFIX}/bin/brew")
    if(NOT arg_HOMEBREW_PREFIX OR NOT EXISTS "${_brew}")
      message(FATAL_ERROR "The package bundles libraries of the Homebrew formulae "
                          "${arg_HOMEBREW_FORMULAE}, but '${_brew}', which is needed to collect "
                          "their license files, does not exist.")
    endif()
    ## --union: the dependencies of any of the formulae, not only those they all share.
    execute_process(COMMAND "${_brew}" deps --installed --union ${arg_HOMEBREW_FORMULAE}
                    OUTPUT_VARIABLE _dependencies
                    RESULT_VARIABLE _result
                    ERROR_VARIABLE _error
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE)
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "'brew deps' failed (${_error}), so the license files of the "
                          "dependencies of the Homebrew formulae ${arg_HOMEBREW_FORMULAE}, whose "
                          "libraries the package bundles, cannot be collected.")
    endif()
    string(REGEX REPLACE "[ \t\r\n]+" ";" _dependencies "${_dependencies}")
    set(_formulae ${arg_HOMEBREW_FORMULAE} ${_dependencies})
    list(REMOVE_DUPLICATES _formulae)

    set(_without_license)
    foreach(_formula IN LISTS _formulae)
      execute_process(COMMAND "${_brew}" --prefix "${_formula}"
                      OUTPUT_VARIABLE _prefix
                      RESULT_VARIABLE _result
                      ERROR_QUIET
                      OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(NOT _result EQUAL 0 OR NOT IS_DIRECTORY "${_prefix}")
        list(APPEND _without_license "${_formula}")
        continue()
      endif()
      file(REAL_PATH "${_prefix}" _keg)
      ## Homebrew copies a formula's top-level license files (COPYING, LICENSE.md, ...)
      ## into the root of its keg.
      file(GLOB _keg_files LIST_DIRECTORIES false "${_keg}/*")
      set(_license_files)
      foreach(_file IN LISTS _keg_files)
        get_filename_component(_name "${_file}" NAME)
        string(TOUPPER "${_name}" _name)
        if(_name MATCHES "^(COPYING|COPYRIGHT|LICENSE|LICENCE|NOTICE)")
          list(APPEND _license_files "${_file}")
        endif()
      endforeach()
      if(_license_files)
        install(FILES ${_license_files}
                DESTINATION "${_licenses_dir}/homebrew/${_formula}"
                COMPONENT share)
      else()
        list(APPEND _without_license "${_formula}")
      endif()
    endforeach()
    list(LENGTH _formulae _count)
    message(STATUS "Packaging the license files of ${_count} Homebrew formulae")
    if(_without_license)
      ## Expected for Qt's own formulae, whose keg holds its licenses only in a LICENSES
      ## directory; the Qt licenses are installed above.
      message(STATUS "No license file in the Homebrew kegs of: ${_without_license}")
    endif()
  endif()
endfunction()

## Installs the licenses of the libraries under src/openms/extern and src/openms/thirdparty,
## and of the third-party code in other OpenMS source files (libkdtree++ in KDTree.h,
## MSNumpress), under share/OpenMS/LICENSES/vendored. A library of src/openms/extern that the
## build takes from vcpkg or the system instead (USE_EXTERNAL_<name>) is left out; for a vcpkg
## port, openms_install_third_party_licenses() installs its license.
function(openms_install_vendored_licenses)
  set(_licenses_dir "${INSTALL_SHARE_DIR}/LICENSES/vendored")
  set(_openms "${OPENMS_HOST_DIRECTORY}/src/openms")
  set(_texts "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/third_party_licenses")

  ## Percolator is under the Apache License 2.0, whose NOTICE file has to accompany it. The
  ## NOTICE file also holds the license of the LIBLINEAR code that Percolator contains.
  install(FILES "${_openms}/thirdparty/percolator/LICENSE-Apache-2.0.txt"
                "${_openms}/thirdparty/percolator/NOTICE-percolator.txt"
          DESTINATION "${_licenses_dir}/percolator"
          COMPONENT share)
  install(FILES "${_texts}/libkdtree/README.txt"
                "${_texts}/libkdtree/Artistic-2.0.txt"
          DESTINATION "${_licenses_dir}/libkdtree"
          COMPONENT share)
  install(FILES "${_texts}/MSNumpress/LICENSE.txt"
          DESTINATION "${_licenses_dir}/MSNumpress"
          COMPONENT share)

  ## <directory in src/openms/extern>/<its license file>
  set(_extern_licenses evergreen/LICENSE GTE/LICENSE Quadtree/LICENSE)
  if(NOT USE_EXTERNAL_JSON)
    list(APPEND _extern_licenses nlohmann_json/LICENSE.MIT)
  endif()
  if(NOT USE_EXTERNAL_SQLITECPP)
    list(APPEND _extern_licenses SQLiteCpp/LICENSE.txt)
  endif()
  if(NOT USE_EXTERNAL_SIMDE)
    list(APPEND _extern_licenses simde/COPYING)
  endif()
  if(NOT USE_EXTERNAL_ISOSPEC)
    list(APPEND _extern_licenses IsoSpec/LICENSE)
  endif()
  if(NOT USE_EXTERNAL_EOLBSPLINE)
    list(APPEND _extern_licenses eol-bspline/LICENSE)
  endif()
  foreach(_license IN LISTS _extern_licenses)
    get_filename_component(_library "${_license}" DIRECTORY)
    install(FILES "${_openms}/extern/${_license}"
            DESTINATION "${_licenses_dir}/${_library}"
            COMPONENT share)
  endforeach()
  if(ENABLE_TDL)
    install(DIRECTORY "${_openms}/extern/tool_description_lib/LICENSES/"
            DESTINATION "${_licenses_dir}/tool_description_lib"
            COMPONENT share)
  endif()
endfunction()

## Installs contrib/<library> for the libraries this build found inside OPENMS_CONTRIB_LIBS,
## judged by where the find modules found them, and README.txt with their versions and
## sources. Libraries from vcpkg, Homebrew or the system are left to
## openms_install_third_party_licenses() or to the system's package manager.
function(openms_install_contrib_licenses)
  if(NOT OPENMS_CONTRIB_LIBS)
    return()
  endif()
  set(_licenses_dir "${INSTALL_SHARE_DIR}/LICENSES/contrib")
  set(_texts "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/third_party_licenses/contrib")
  cmake_path(SET _contrib NORMALIZE "${OPENMS_CONTRIB_LIBS}")

  ## <folder in _texts>=<variables that hold where a find module found the library>. Arrow
  ## builds and links snappy, zstd, Thrift, xsimd, RapidJSON and mimalloc itself, so they
  ## come with it. Its S3 support (the AWS SDK) is built on Unix but not linked by OpenMS.
  set(_libraries
    "boost=Boost_DIR"
    "bzip2=BZIP2_INCLUDE_DIR"
    "zlib=ZLIB_INCLUDE_DIR"
    "curl=CURL_INCLUDE_DIR,CURL_DIR"
    "eigen=Eigen3_DIR,EIGEN3_INCLUDE_DIR"
    "libsvm=LIBSVM_INCLUDE_DIR"
    "libzip=LIBZIP_INCLUDE_DIR"
    "xerces-c=XercesC_INCLUDE_DIR"
    "coinmp=COIN_INCLUDE_DIR"
    "arrow,snappy,zstd,thrift,xsimd,rapidjson,mimalloc=Arrow_DIR")
  set(_installed)
  foreach(_entry IN LISTS _libraries)
    string(REPLACE "=" ";" _entry "${_entry}")
    list(GET _entry 0 _folders)
    list(GET _entry 1 _variables)
    string(REPLACE "," ";" _folders "${_folders}")
    string(REPLACE "," ";" _variables "${_variables}")
    set(_from_contrib FALSE)
    foreach(_variable IN LISTS _variables)
      if(${_variable})
        cmake_path(IS_PREFIX _contrib "${${_variable}}" NORMALIZE _inside)
        if(_inside)
          set(_from_contrib TRUE)
        endif()
      endif()
    endforeach()
    if(_from_contrib)
      foreach(_folder IN LISTS _folders)
        install(DIRECTORY "${_texts}/${_folder}"
                DESTINATION "${_licenses_dir}"
                COMPONENT share)
      endforeach()
      list(APPEND _installed ${_folders})
    endif()
  endforeach()
  if(_installed)
    install(FILES "${_texts}/README.txt"
            DESTINATION "${_licenses_dir}"
            COMPONENT share)
    list(JOIN _installed ", " _installed)
    message(STATUS "Licenses of the contrib libraries to install: ${_installed}")
  endif()
endfunction()
