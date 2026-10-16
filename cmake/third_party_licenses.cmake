# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

## The packages and the pyOpenMS wheels bundle third-party libraries, and most of their
## licenses require the license text to accompany the binaries. The functions below install
## these texts under share/OpenMS/LICENSES, in the component 'share' that every installation
## type includes.
##
## openms_install_vcpkg_licenses() installs vcpkg/<port>.txt: the license text of every vcpkg
## port of the target triplet, which vcpkg installs as share/<port>/copyright, except the
## build helpers such as vcpkg-cmake. An installation carries the code of each port: as a
## library next to libOpenMS (install(RUNTIME_DEPENDENCY_SET) in package_general.cmake),
## linked into it (the Windows triplets and those of the Linux wheels are static), or
## compiled in (header-only ports).
## The top-level CMakeLists.txt calls it in every configuration, so that the wheels built
## with vcpkg carry the texts too.
##
## openms_install_third_party_licenses() installs the texts of what only the packages bundle;
## cmake/package_general.cmake calls it:
## - Qt/: Qt is not a vcpkg port. The Windows and macOS packages bundle it; the DEB depends
##   on the distribution's Qt instead. OpenMS uses Qt under the LGPL version 3, which
##   requires its text, the text of the GPL version 3 that it supplements, and a note on
##   where to get the source code of the bundled Qt.
## - homebrew/<formula>/: on macOS, Qt and the libraries it needs (glib, ICU, freetype, ...)
##   and the OpenMP runtime libomp come from Homebrew, and the package bundles them. Homebrew
##   installs the license files of a formula into its keg, and these are installed for the
##   given formulae and every installed formula they depend on.
##
## openms_install_system_library_licenses() installs system/<package>/: on Linux, the
## libraries the package takes from the build machine's distribution rather than from vcpkg
## (libgfortran, which vcpkg's LAPACK needs), with the copyright file of the Debian package
## each comes from; cmake/package_general.cmake calls it.
##
## openms_install_vendored_licenses() installs the licenses of the third-party code that
## OpenMS carries in its own source tree and compiles into its libraries (vendored/<name>/).
## Every installation of libOpenMS contains that code, the Python wheels too, so the
## top-level CMakeLists.txt calls it in every configuration, not only for packages.
##
## Each of them also records the files it installs, for share/OpenMS/THIRD-PARTY-NOTICES.txt,
## which gathers the texts in one file (cmake/third_party_notices.cmake).

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

## Installs share/OpenMS/LICENSES/vcpkg/<port>.txt for every port of the target triplet, if
## the build takes its dependencies from vcpkg. Left out are the build helpers such as
## vcpkg-cmake: ports that install nothing but their own share/<port> folder (CMake scripts
## for the builds of other ports), so no installation contains anything of them. They are in
## the target triplet's tree when it is also the host triplet, as in the Linux builds. vcpkg
## lists the files of each port in vcpkg/info/<port>_<version>_<triplet>.list, one
## <triplet>/<path> per line and folders with a trailing slash; a port without such a list
## keeps its license.
function(openms_install_vcpkg_licenses)
  if(NOT (OPENMS_USE_VCPKG AND VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET))
    return()
  endif()
  file(GLOB _copyrights LIST_DIRECTORIES false
       "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/*/copyright")
  set(_count 0)
  set(_helpers)
  foreach(_copyright IN LISTS _copyrights)
    get_filename_component(_port_dir "${_copyright}" DIRECTORY)
    get_filename_component(_port "${_port_dir}" NAME)
    ## Port and triplet names consist of lowercase letters, digits and dashes, so they need
    ## no escaping in the patterns.
    file(GLOB _lists "${VCPKG_INSTALLED_DIR}/vcpkg/info/${_port}_*_${VCPKG_TARGET_TRIPLET}.list")
    if(_lists)
      set(_files)
      foreach(_list IN LISTS _lists)
        file(STRINGS "${_list}" _entries REGEX "[^/]$")
        list(APPEND _files ${_entries})
      endforeach()
      list(FILTER _files EXCLUDE REGEX "^${VCPKG_TARGET_TRIPLET}/share/${_port}/")
      if(NOT _files)
        list(APPEND _helpers "${_port}")
        continue()
      endif()
    endif()
    install(FILES "${_copyright}"
            DESTINATION "${INSTALL_SHARE_DIR}/LICENSES/vcpkg"
            RENAME "${_port}.txt"
            COMPONENT share)
    openms_add_third_party_notice("LICENSES/vcpkg/${_port}.txt" "${_copyright}")
    math(EXPR _count "${_count} + 1")
  endforeach()
  message(STATUS "Installing the license texts of ${_count} vcpkg ports")
  if(_helpers)
    list(JOIN _helpers ", " _helpers)
    message(STATUS "Not installing those of the build helpers ${_helpers}, which install "
                   "nothing but their share folder")
  endif()
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
    set(_qt_files "${_qt_licenses}/LGPL-3.0-only.txt"
                  "${_qt_licenses}/GPL-3.0-only.txt"
                  "${PROJECT_BINARY_DIR}/third_party_licenses/Qt/README.txt")
    install(FILES ${_qt_files}
            DESTINATION "${_licenses_dir}/Qt"
            COMPONENT share)
    foreach(_file IN LISTS _qt_files)
      get_filename_component(_name "${_file}" NAME)
      openms_add_third_party_notice("LICENSES/Qt/${_name}" "${_file}")
    endforeach()
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
    ## Both the dependencies recorded in the installed kegs (--installed) and those the
    ## formulae declare are collected. The recorded ones alone miss formulae: in the Release
    ## run of 2026-10-08, 'brew deps --installed' left out exactly the dependencies that
    ## 'brew install qtbase' had upgraded in the same run (glib, harfbuzz, pcre2, libpng,
    ## cairo, ...), although the package bundles their libraries. homebrew-core's CI checks
    ## that a formula declares every formula it links against, so the declared ones cover
    ## what is bundled.
    set(_dependencies)
    foreach(_installed_only "--installed" "")
      execute_process(COMMAND "${_brew}" deps ${_installed_only} --union ${arg_HOMEBREW_FORMULAE}
                      OUTPUT_VARIABLE _deps
                      RESULT_VARIABLE _result
                      ERROR_VARIABLE _error
                      OUTPUT_STRIP_TRAILING_WHITESPACE
                      ERROR_STRIP_TRAILING_WHITESPACE)
      if(NOT _result EQUAL 0)
        message(FATAL_ERROR "'brew deps ${_installed_only} --union' failed (${_error}), so the "
                            "license files of the dependencies of the Homebrew formulae "
                            "${arg_HOMEBREW_FORMULAE}, whose libraries the package bundles, "
                            "cannot be collected.")
      endif()
      string(REGEX REPLACE "[ \t\r\n]+" ";" _deps "${_deps}")
      list(APPEND _dependencies ${_deps})
    endforeach()
    set(_formulae ${arg_HOMEBREW_FORMULAE} ${_dependencies})
    list(REMOVE_DUPLICATES _formulae)
    list(SORT _formulae)

    set(_without_license)
    set(_not_installed)
    foreach(_formula IN LISTS _formulae)
      execute_process(COMMAND "${_brew}" --prefix "${_formula}"
                      OUTPUT_VARIABLE _prefix
                      RESULT_VARIABLE _result
                      ERROR_QUIET
                      OUTPUT_STRIP_TRAILING_WHITESPACE)
      ## 'brew --prefix' names the opt link of a formula whether or not it is installed. A
      ## declared dependency that is not installed (an optional one, or one the formula
      ## gained after it was installed) cannot have a library in the package.
      if(NOT _result EQUAL 0 OR NOT IS_DIRECTORY "${_prefix}")
        list(APPEND _not_installed "${_formula}")
        continue()
      endif()
      file(REAL_PATH "${_prefix}" _keg)
      ## Homebrew copies a formula's top-level license files (COPYING, LICENSE.md, ...)
      ## into the root of its keg; some formulae install them to share/doc/<formula> or
      ## share/licenses/<formula> instead. The first of these folders that has any is used,
      ## as tools/ci/collect_wheel_licenses.py does for the wheels.
      set(_license_files)
      foreach(_folder "${_keg}" "${_keg}/share/doc/${_formula}" "${_keg}/share/licenses/${_formula}")
        file(GLOB _keg_files LIST_DIRECTORIES false "${_folder}/*")
        foreach(_file IN LISTS _keg_files)
          get_filename_component(_name "${_file}" NAME)
          string(TOUPPER "${_name}" _name)
          if(_name MATCHES "^(COPYING|COPYRIGHT|LICENSE|LICENCE|NOTICE)")
            list(APPEND _license_files "${_file}")
          endif()
        endforeach()
        if(_license_files)
          break()
        endif()
      endforeach()
      if(_license_files)
        install(FILES ${_license_files}
                DESTINATION "${_licenses_dir}/homebrew/${_formula}"
                COMPONENT share)
        foreach(_file IN LISTS _license_files)
          get_filename_component(_name "${_file}" NAME)
          openms_add_third_party_notice("LICENSES/homebrew/${_formula}/${_name}" "${_file}")
        endforeach()
      else()
        list(APPEND _without_license "${_formula}")
      endif()
    endforeach()
    list(LENGTH _formulae _count)
    list(LENGTH _not_installed _count_not_installed)
    math(EXPR _count "${_count} - ${_count_not_installed}")
    message(STATUS "Packaging the license files of ${_count} Homebrew formulae")
    if(_not_installed)
      message(STATUS "Declared Homebrew dependencies that are not installed: ${_not_installed}")
    endif()
    if(_without_license)
      ## Expected for Qt's own formulae, whose keg holds its licenses only in a LICENSES
      ## directory; the Qt licenses are installed above.
      message(STATUS "No license file in the Homebrew kegs of: ${_without_license}")
    endif()
  endif()
endfunction()

## Sets <var> to the libraries <library> needs (the NEEDED entries of its dynamic section).
function(_openms_needed_libraries var objdump library)
  execute_process(COMMAND "${objdump}" -p "${library}"
                  OUTPUT_VARIABLE _output
                  RESULT_VARIABLE _result
                  ERROR_QUIET)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "'${objdump} -p ${library}' failed, so the libraries it needs, whose "
                        "licenses the package has to ship, are unknown.")
  endif()
  string(REGEX MATCHALL "NEEDED +[^ \t\r\n]+" _entries "${_output}")
  list(TRANSFORM _entries REPLACE "^NEEDED +" "")
  set(${var} "${_entries}" PARENT_SCOPE)
endfunction()

## openms_install_system_library_licenses(ROOTS <library>... [EXCLUDE_REGEXES <regex>...])
##   ROOTS            shared libraries that the package bundles and that come with their
##                    licenses, i.e. those of the vcpkg tree; what they need, directly or
##                    through each other, is followed
##   EXCLUDE_REGEXES  paths of libraries the package does not bundle, the
##                    POST_EXCLUDE_REGEXES of its runtime dependency set (the C and C++
##                    runtime, Qt, ...)
## The libraries the roots need that neither another root provides nor an exclusion drops
## come from the build machine's distribution, and the package bundles them: on Ubuntu, the
## GCC Fortran runtime libgfortran, which vcpkg's LAPACK needs. They are found as the
## compiler finds them (-print-file-name), and the libraries they need are followed in turn.
## For each Debian package that installed one, its copyright file goes to
## LICENSES/system/<package>/copyright. Debian's copyright files refer to the full texts of
## common licenses in /usr/share/common-licenses instead of including them, so the texts
## they name go next to it. LICENSES/system/INDEX.txt lists the
## libraries with the package, version and source package each comes from. Configuring
## fails rather than leave a license out: without dpkg, for a library the compiler cannot
## find, or for one no package owns.
function(openms_install_system_library_licenses)
  cmake_parse_arguments(PARSE_ARGV 0 arg "" "" "ROOTS;EXCLUDE_REGEXES")
  set(_provided)
  set(_queue)
  foreach(_root IN LISTS arg_ROOTS)
    get_filename_component(_name "${_root}" NAME)
    list(APPEND _provided "${_name}")
    if(NOT IS_SYMLINK "${_root}")
      list(APPEND _queue "${_root}")
    endif()
  endforeach()
  if(NOT _queue)
    return()
  endif()
  set(_objdump "${CMAKE_OBJDUMP}")
  if(NOT _objdump)
    find_program(OPENMS_OBJDUMP_EXECUTABLE objdump)
    set(_objdump "${OPENMS_OBJDUMP_EXECUTABLE}")
  endif()
  if(NOT _objdump)
    message(FATAL_ERROR "objdump is needed to find the libraries of the build machine that "
                        "the package bundles, whose licenses it has to ship.")
  endif()

  ## Walk the dependencies: <name>|<path> for every library of the build machine.
  set(_seen ${_provided})
  set(_system)
  set(_unresolved)
  while(_queue)
    list(POP_FRONT _queue _library)
    _openms_needed_libraries(_needed "${_objdump}" "${_library}")
    foreach(_name IN LISTS _needed)
      if(_name IN_LIST _seen)
        continue()
      endif()
      list(APPEND _seen "${_name}")
      execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "-print-file-name=${_name}"
                      OUTPUT_VARIABLE _path
                      OUTPUT_STRIP_TRAILING_WHITESPACE
                      ERROR_QUIET)
      if(IS_ABSOLUTE "${_path}" AND EXISTS "${_path}")
        file(REAL_PATH "${_path}" _real)
        set(_candidates "${_path}" "${_real}")
      else()
        ## Not where the compiler looks, but the install step may still find and bundle it
        ## (through a RUNPATH, ld.so.conf or its DIRECTORIES), so only an exclusion that
        ## matches the name lets it pass.
        set(_real)
        set(_candidates "/${_name}")
      endif()
      set(_excluded FALSE)
      foreach(_regex IN LISTS arg_EXCLUDE_REGEXES)
        foreach(_candidate IN LISTS _candidates)
          if(_candidate MATCHES "${_regex}")
            set(_excluded TRUE)
          endif()
        endforeach()
      endforeach()
      if(_excluded)
        continue()
      elseif(NOT _real)
        list(APPEND _unresolved "${_name} (needed by ${_library})")
        continue()
      endif()
      list(APPEND _system "${_name}|${_real}")
      list(APPEND _queue "${_real}")
    endforeach()
  endwhile()
  if(_unresolved)
    list(JOIN _unresolved ", " _list)
    message(FATAL_ERROR "The compiler (-print-file-name) cannot find ${_list}, so whether the "
                        "package bundles them, and the licenses it then has to ship, cannot be "
                        "checked.")
  endif()
  if(NOT _system)
    return()
  endif()

  find_program(OPENMS_DPKG_QUERY_EXECUTABLE dpkg-query)
  set(_dpkg_query "${OPENMS_DPKG_QUERY_EXECUTABLE}")
  if(NOT _dpkg_query)
    string(REPLACE "|" " from " _list "${_system}")
    message(FATAL_ERROR "The package bundles libraries of the build machine (${_list}), but "
                        "without dpkg-query the packages they come from, and so their "
                        "licenses, cannot be found.")
  endif()
  set(_system_dir "${INSTALL_SHARE_DIR}/LICENSES/system")
  set(_packages)
  set(_index)
  foreach(_entry IN LISTS _system)
    string(REPLACE "|" ";" _entry "${_entry}")
    list(GET _entry 0 _name)
    list(GET _entry 1 _real)
    ## With a merged /usr, dpkg may know the file under /lib instead of /usr/lib.
    set(_package)
    set(_architecture)
    foreach(_candidate "${_real}" "")
      if(NOT _candidate)
        string(REGEX REPLACE "^/usr/" "/" _candidate "${_real}")
      endif()
      execute_process(COMMAND "${_dpkg_query}" -S "${_candidate}"
                      OUTPUT_VARIABLE _owner
                      RESULT_VARIABLE _result
                      ERROR_QUIET)
      ## "libgfortran5:amd64: /usr/lib/x86_64-linux-gnu/libgfortran.so.5.0.0"
      if(_result EQUAL 0 AND _owner MATCHES "^([^:, \n]+)(:[^: \n]+)?: ")
        set(_package "${CMAKE_MATCH_1}")
        set(_architecture "${CMAKE_MATCH_2}")
        break()
      endif()
    endforeach()
    if(NOT _package)
      message(FATAL_ERROR "The package bundles ${_name} (${_real}), which no Debian package "
                          "owns, so its license cannot be found.")
    endif()
    ## With the architecture qualifier, as -W lists every installed architecture of a
    ## multiarch package otherwise.
    execute_process(COMMAND "${_dpkg_query}" -W
                            "-f=\${Version}|\${source:Package}|\${source:Version}"
                            "${_package}${_architecture}"
                    OUTPUT_VARIABLE _versions
                    RESULT_VARIABLE _result
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(NOT _result EQUAL 0 OR NOT _versions MATCHES "^[^|\n]+\\|[^|\n]+\\|[^|\n]+$")
      message(FATAL_ERROR "'dpkg-query -W ${_package}${_architecture}' failed or returned "
                          "'${_versions}', so the version of ${_name}'s package cannot be "
                          "recorded.")
    endif()
    string(REPLACE "|" ";" _versions "${_versions}")
    list(GET _versions 0 _version)
    list(GET _versions 1 _source)
    list(GET _versions 2 _source_version)
    list(APPEND _index "${_name}\t${_package} ${_version} (source package ${_source} ${_source_version})")
    if(_package IN_LIST _packages)
      continue()
    endif()
    list(APPEND _packages "${_package}")
    ## The doc folder of a library package is often a link to that of its source package's
    ## base package (libgfortran5 -> gcc-14-base).
    set(_copyright "/usr/share/doc/${_package}/copyright")
    if(NOT EXISTS "${_copyright}")
      message(FATAL_ERROR "The package bundles ${_name} from the Debian package ${_package}, "
                          "which has no ${_copyright}.")
    endif()
    file(REAL_PATH "${_copyright}" _copyright)
    install(FILES "${_copyright}"
            DESTINATION "${_system_dir}/${_package}"
            COMPONENT share)
    openms_add_third_party_notice("LICENSES/system/${_package}/copyright" "${_copyright}")
    file(READ "${_copyright}" _text)
    string(REGEX MATCHALL "/usr/share/common-licenses/[A-Za-z0-9.+-]*[A-Za-z0-9+]" _references "${_text}")
    list(REMOVE_DUPLICATES _references)
    foreach(_reference IN LISTS _references)
      if(NOT EXISTS "${_reference}")
        message(FATAL_ERROR "The copyright file of the Debian package ${_package} refers to "
                            "${_reference}, which does not exist.")
      endif()
      get_filename_component(_text_name "${_reference}" NAME)
      ## GPL is a link to GPL-3; it is installed as a file under the name the text uses.
      file(REAL_PATH "${_reference}" _text_file)
      install(FILES "${_text_file}"
              DESTINATION "${_system_dir}/${_package}"
              RENAME "${_text_name}"
              COMPONENT share)
      openms_add_third_party_notice("LICENSES/system/${_package}/${_text_name}" "${_text_file}")
    endforeach()
  endforeach()

  list(SORT _index)
  list(JOIN _index "\n" _index)
  set(_index_file "${PROJECT_BINARY_DIR}/third_party_licenses/system/INDEX.txt")
  file(WRITE "${_index_file}"
       "Libraries of the build machine's distribution that this package bundles, and the Debian\n"
       "package each comes from. LICENSES/system/<package>/ holds the package's copyright file\n"
       "and the texts in /usr/share/common-licenses that it refers to. The source code of a\n"
       "package is in the archive of the distribution: apt-get source <source package>=<version>.\n\n"
       "${_index}\n")
  install(FILES "${_index_file}"
          DESTINATION "${_system_dir}"
          COMPONENT share)
  openms_add_third_party_notice("LICENSES/system/INDEX.txt" "${_index_file}")
  list(JOIN _packages ", " _packages)
  message(STATUS "Packaging the licenses of the libraries the package takes from the build "
                 "machine, from the Debian packages ${_packages}")
endfunction()

## Installs the licenses of the libraries under src/openms/extern and src/openms/thirdparty,
## and of the third-party code in other OpenMS source files (libkdtree++ in KDTree.h,
## MSNumpress), under share/OpenMS/LICENSES/vendored. A library of src/openms/extern that the
## build takes from vcpkg or the system instead (USE_EXTERNAL_<name>) is left out; for a vcpkg
## port, openms_install_vcpkg_licenses() installs its license.
## Installs <files> to share/OpenMS/LICENSES/vendored/<library> and records them for
## THIRD-PARTY-NOTICES.txt.
function(_openms_install_vendored_license library)
  install(FILES ${ARGN}
          DESTINATION "${INSTALL_SHARE_DIR}/LICENSES/vendored/${library}"
          COMPONENT share)
  foreach(_file IN LISTS ARGN)
    get_filename_component(_name "${_file}" NAME)
    openms_add_third_party_notice("LICENSES/vendored/${library}/${_name}" "${_file}")
  endforeach()
endfunction()

function(openms_install_vendored_licenses)
  set(_licenses_dir "${INSTALL_SHARE_DIR}/LICENSES/vendored")
  set(_openms "${OPENMS_HOST_DIRECTORY}/src/openms")
  set(_texts "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/third_party_licenses")

  ## Percolator is under the Apache License 2.0, whose NOTICE file has to accompany it. The
  ## NOTICE file also holds the license of the LIBLINEAR code that Percolator contains.
  _openms_install_vendored_license(percolator
                                   "${_openms}/thirdparty/percolator/LICENSE-Apache-2.0.txt"
                                   "${_openms}/thirdparty/percolator/NOTICE-percolator.txt")
  _openms_install_vendored_license(libkdtree
                                   "${_texts}/libkdtree/README.txt"
                                   "${_texts}/libkdtree/Artistic-2.0.txt")
  _openms_install_vendored_license(MSNumpress "${_texts}/MSNumpress/LICENSE.txt")

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
    _openms_install_vendored_license("${_library}" "${_openms}/extern/${_license}")
  endforeach()
  if(ENABLE_TDL)
    install(DIRECTORY "${_openms}/extern/tool_description_lib/LICENSES/"
            DESTINATION "${_licenses_dir}/tool_description_lib"
            COMPONENT share)
    openms_add_third_party_notice_directory("LICENSES/vendored/tool_description_lib"
                                            "${_openms}/extern/tool_description_lib/LICENSES")
  endif()
endfunction()
