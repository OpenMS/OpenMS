# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

# a collection of wrapper for install functions that allows easier usage
# throughout the OpenMS build system

#------------------------------------------------------------------------------
# The installed package is layered. Each layer is one export set (the file
# OpenMSConfig.cmake includes) and one pair of install components, so an
# installation may stop at any layer and the package of what is installed stays
# consistent (CMake refuses an export whose library files are missing):
#
#   layer  export set       libraries      cmake files   targets
#   core   OpenMSTargets    library        cmake         OpenMS::OpenMS, OpenMS::OpenSwathAlgo
#                                                        and the bundled third-party libraries
#   CLI    OpenMSCLITargets library_cli    cmake_cli     OpenMS::OpenMS_CLI (links the core layer)
#   GUI    OpenMSGUITargets library_gui    cmake_gui     OpenMS::OpenMS_GUI (links the CLI layer)
#
# Headers keep their own components (OpenMS_headers, OpenMS_CLI_headers,
# OpenMS_GUI_headers, OpenSwathAlgo_headers, thirdparty_headers). A core-only
# installation (e.g. the one the pyOpenMS wheels are built against) is
# library + the core header components + share + cmake.
set(OPENMS_EXPORT_SET "OpenMSTargets")
set(OPENMS_CLI_EXPORT_SET "OpenMSCLITargets")
set(OPENMS_GUI_EXPORT_SET "OpenMSGUITargets")
set(OPENMS_EXPORT_SETS ${OPENMS_EXPORT_SET} ${OPENMS_CLI_EXPORT_SET} ${OPENMS_GUI_EXPORT_SET})

# the install components of each export set: <set>_LIBRARY_COMPONENT holds the
# libraries, <set>_CMAKE_COMPONENT the exported target files
set(${OPENMS_EXPORT_SET}_LIBRARY_COMPONENT library)
set(${OPENMS_EXPORT_SET}_CMAKE_COMPONENT cmake)
set(${OPENMS_CLI_EXPORT_SET}_LIBRARY_COMPONENT library_cli)
set(${OPENMS_CLI_EXPORT_SET}_CMAKE_COMPONENT cmake_cli)
set(${OPENMS_GUI_EXPORT_SET}_LIBRARY_COMPONENT library_gui)
set(${OPENMS_GUI_EXPORT_SET}_CMAKE_COMPONENT cmake_gui)

#------------------------------------------------------------------------------
# Installs the library lib_target_name into the library component of its export
# set and adds it to that export set.
#
# install_library(<target> [HEADERS] [EXPORT_SET <set>])
#
# @param lib_target_name The target name of the library that should be installed
# @param EXPORT_SET      One of ${OPENMS_EXPORT_SETS}; the core set OpenMSTargets
#                        (component 'library') when omitted
# @param HEADERS         Install the public HEADERS file set in <target>_headers
function(install_library lib_target_name)
    cmake_parse_arguments(_install_library "HEADERS" "EXPORT_SET" "" ${ARGN})
    if(_install_library_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "install_library(${lib_target_name}): unexpected arguments ${_install_library_UNPARSED_ARGUMENTS}")
    endif()
    if(NOT _install_library_EXPORT_SET)
      set(_install_library_EXPORT_SET ${OPENMS_EXPORT_SET})
    endif()
    if(NOT _install_library_EXPORT_SET IN_LIST OPENMS_EXPORT_SETS)
      message(FATAL_ERROR "install_library(${lib_target_name}): unknown export set '${_install_library_EXPORT_SET}' (known: ${OPENMS_EXPORT_SETS})")
    endif()
    set(_component ${${_install_library_EXPORT_SET}_LIBRARY_COMPONENT})
    set(_header_install_args)
    if(_install_library_HEADERS)
      set(_header_install_args FILE_SET HEADERS DESTINATION ${INSTALL_INCLUDE_DIR}
        COMPONENT ${lib_target_name}_headers)
    endif()
    install(TARGETS ${lib_target_name}
      RUNTIME_DEPENDENCY_SET OPENMS_DEPS
      EXPORT ${_install_library_EXPORT_SET}
      LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT ${_component}
      ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT ${_component}
      RUNTIME DESTINATION ${INSTALL_LIB_DIR} COMPONENT ${_component}
      ${_header_install_args}
      )
endfunction()

#------------------------------------------------------------------------------
# Reject JSON dependencies in installed headers even when compiler header checks
# could resolve them through a system or shared dependency include directory.
# Keep this configure-time guard independent of the header installation method.
function(openms_validate_public_headers)
  foreach(_header IN LISTS ARGN)
    # nlohmann::json is a PRIVATE dependency of libOpenMS: downstream builds do not have its
    # include directory, so no installed header may include it (configure-time check only).
    # Header lists are relative to the calling source directory.
    # Generated headers in the build tree may not exist yet at this point; they come from OpenMS'
    # own templates, so only what is already on disk is scanned.
    set(_header_to_scan "${_header}")
    if (NOT IS_ABSOLUTE "${_header_to_scan}")
      set(_header_to_scan "${CMAKE_CURRENT_SOURCE_DIR}/${_header_to_scan}")
    endif()
    if (EXISTS "${_header_to_scan}" AND NOT IS_DIRECTORY "${_header_to_scan}")
      file(STRINGS "${_header_to_scan}" _nlohmann_json_hits REGEX "nlohmann/json")
      if (_nlohmann_json_hits)
        message(FATAL_ERROR "Installed header ${_header} includes nlohmann/json. Keep JSON types out of "
                            "public headers: move the header under source/ or expose value types instead.")
      endif()
    endif()
  endforeach()
endfunction()

#------------------------------------------------------------------------------
# Installs the tool tool_target_name
# @param tool_target_name The target name of the tool that should be installed
macro(install_tool tool_target_name)
    install(TARGETS ${tool_target_name} RUNTIME_DEPENDENCY_SET OPENMS_DEPS
      RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT Applications
      BUNDLE DESTINATION ${INSTALL_BIN_DIR} COMPONENT Applications
      )
endmacro()

#------------------------------------------------------------------------------
# Installs a given directory
# @param directory The directory to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the directory belongs
macro(install_directory directory destination component)
    install(DIRECTORY ${directory}
      DESTINATION ${destination}
      COMPONENT ${component}
      FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                            GROUP_READ
                            WORLD_READ
      DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                            GROUP_EXECUTE GROUP_READ
                            WORLD_EXECUTE WORLD_READ
      REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
      REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
        )
endmacro()

#------------------------------------------------------------------------------
# Installs a given file
# @param directory The file to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the file belongs
macro(install_file file destination component)
    install(FILES ${file}
      DESTINATION ${destination}
      COMPONENT ${component})
endmacro()

#------------------------------------------------------------------------------
# Execute the given code while executing the install target
# @param code_snippet The code to execute
# @param component The component to which the code will be associated
macro(install_code code_snippet component)
    install(CODE ${code_snippet}
            COMPONENT ${component})
endmacro()

#------------------------------------------------------------------------------
# Installs the exported target information of every export set that has
# targets: one <set>.cmake file per layer in the cmake component of that layer
# (cmake, cmake_cli, cmake_gui), so an installation without a layer has no file
# claiming its libraries. Consumers see every target as OpenMS::<target>
# (OpenMS::OpenMS, OpenMS::OpenSwathAlgo, OpenMS::OpenMS_CLI, ...); a target of
# one layer refers to the targets of the layer below by these names, and
# OpenMSConfig.cmake includes the files in layer order and adds un-namespaced
# aliases for the OpenMS libraries.
macro(install_export_targets )
    foreach(_export_set IN LISTS OPENMS_EXPORT_SETS)
      if(_OPENMS_EXPORT_TARGETS_${_export_set})
        install(EXPORT ${_export_set}
                NAMESPACE OpenMS::
                DESTINATION ${INSTALL_CMAKE_DIR}
                COMPONENT ${${_export_set}_CMAKE_COMPONENT})
      endif()
    endforeach()
endmacro()

#------------------------------------------------------------------------------
# Returns in ${out_var} the .NET runtime identifiers (RIDs) whose native assets
# this build can load, spelled as NuGet spells them, <os>-<arch>:
# https://learn.microsoft.com/dotnet/core/rid-catalog
function(openms_target_dotnet_rids out_var)
  if(WIN32)
    set(_os "win")
  elseif(APPLE)
    set(_os "osx")
  else()
    set(_os "linux")
  endif()

  ## CMAKE_SYSTEM_PROCESSOR spells the same architecture differently per
  ## platform (x86_64 vs AMD64, aarch64 vs arm64). On macOS
  ## CMAKE_OSX_ARCHITECTURES overrides it and may name more than one.
  set(_processors "${CMAKE_SYSTEM_PROCESSOR}")
  if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    set(_processors ${CMAKE_OSX_ARCHITECTURES})
  endif()

  set(_rids)
  foreach(_processor IN LISTS _processors)
    string(TOLOWER "${_processor}" _processor)
    if(_processor MATCHES "^(x86_64|amd64|x64)$")
      set(_arch "x64")
    elseif(_processor MATCHES "^(aarch64|arm64)$")
      set(_arch "arm64")
    elseif(_processor MATCHES "^(i[3-6]86|x86)$")
      set(_arch "x86")
    elseif(_processor MATCHES "^arm")
      set(_arch "arm")
    else()
      ## Unknown architecture: pass it through rather than guess. A RID nothing
      ## matches only means nothing is kept, which install_thirdparty_folder
      ## reports below.
      set(_arch "${_processor}")
    endif()
    list(APPEND _rids "${_os}-${_arch}")
  endforeach()

  ## An arm64 Mac runs x86_64 processes under Rosetta, so an x86_64 mono on one
  ## still asks for the osx-x64 helper. Mach-O files never reach dpkg-shlibdeps,
  ## so keeping it costs nothing where the exclusion below actually matters.
  if(APPLE AND "osx-arm64" IN_LIST _rids)
    list(APPEND _rids "osx-x64")
  endif()

  list(REMOVE_DUPLICATES _rids)
  set(${out_var} "${_rids}" PARENT_SCOPE)
endfunction()

#------------------------------------------------------------------------------
# Installs Thirdparty folders with executables
macro(install_thirdparty_folder foldername)
  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/${foldername})
    ## The .NET tools here (ThermoRawFileParser) carry NuGet's 'runtimes/<rid>/'
    ## layout: one folder of native helper libraries per runtime identifier the
    ## NuGet package supports. Mono.Unix ships nine of them -- android-arm,
    ## android-arm64, android-x64, android-x86, linux-arm, linux-arm64,
    ## linux-x64, osx-arm64 and osx-x64 -- and at most the one matching the
    ## machine can ever be loaded. ThermoRawFileParser's own Mono.Unix.dll.config
    ## maps the library to runtimes/linux-x64, runtimes/osx-x64 and
    ## runtimes/osx-arm64, and to nothing else.
    ##
    ## Seven of the nine carry a .so, so CPackDeb hands all seven to
    ## dpkg-shlibdeps and five are foreign on any given host; the two osx ones
    ## are Mach-O and never reach it.
    ##
    ## The rest is not merely dead weight. Any foreign-architecture ELF in the
    ## staging tree makes dpkg-shlibdeps fail with "cannot find library
    ## libc.so.6 needed by ... (ELF format: ...)", an error --ignore-missing-info
    ## does not cover, so CPack aborts before writing the .deb. That is what
    ## forced CPACK_DEBIAN_PACKAGE_SHLIBDEPS off once already (on in #10202,
    ## reverted in #10207). Install the RIDs this build targets and leave the
    ## others behind; cmake/package_deb.cmake depends on that.
    openms_target_dotnet_rids(_target_rids)
    set(_foreign_rid_excludes)
    set(_present_rids)
    set(_kept_rids)
    set(_runtimes_dir "${SEARCH_ENGINES_DIRECTORY}/${foldername}/runtimes")
    file(GLOB _rid_entries RELATIVE "${_runtimes_dir}" "${_runtimes_dir}/*")
    foreach(_rid IN LISTS _rid_entries)
      if(NOT IS_DIRECTORY "${_runtimes_dir}/${_rid}")
        continue()
      endif()
      list(APPEND _present_rids "${_rid}")
      if("${_rid}" IN_LIST _target_rids)
        list(APPEND _kept_rids "${_rid}")
      else()
        ## Matches the RID directory itself, so install() never descends into
        ## it. The (/|$) tail is what keeps 'android-arm' from also matching
        ## 'android-arm64'.
        list(APPEND _foreign_rid_excludes REGEX "/runtimes/${_rid}(/|$)" EXCLUDE)
      endif()
    endforeach()
    if(_present_rids AND NOT _kept_rids)
      ## Nothing here this build could load, so drop the tree instead of leaving
      ## an empty runtimes/ behind. Windows is that case: Mono.Unix is a Unix
      ## helper and the vendored tool carries no win-* runtime at all.
      set(_foreign_rid_excludes REGEX "/runtimes(/|$)" EXCLUDE)
    endif()
    if(_foreign_rid_excludes)
      message(STATUS "${foldername}: installing .NET runtimes ${_kept_rids} of ${_present_rids}")
    endif()

    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/${foldername}
            DESTINATION           ${INSTALL_SHARE_DIR}/THIRDPARTY
            COMPONENT             ${foldername}
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            ${_foreign_rid_excludes}
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
    list(APPEND THIRDPARTY_COMPONENT_GROUP ${foldername})
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs Qt plugins. Prefix can be the usual CMAKE_INSTALL_PREFIX or
# if you install to a bundle, the app folder
# Fills _qt_plugins_var with paths to be used e.g. with fixup_bundle at INSTALL time.
macro(install_qt6_plugin_installdir _qt_plugin_name _qt_plugins_var _targetpath _component)
  get_target_property(_qt_plugin_path "${_qt_plugin_name}" LOCATION)
  if(EXISTS "${_qt_plugin_path}")
    get_filename_component(_qt_plugin_file "${_qt_plugin_path}" NAME)
    get_filename_component(_qt_plugin_type "${_qt_plugin_path}" PATH)
    get_filename_component(_qt_plugin_type "${_qt_plugin_type}" NAME)
    set(_qt_plugin_dest "${_targetpath}/${_qt_plugin_type}")
    install(FILES "${_qt_plugin_path}"
      DESTINATION "${_qt_plugin_dest}"
      COMPONENT ${_component})
    set(${_qt_plugins_var}
      "${${_qt_plugins_var}};\${CMAKE_INSTALL_PREFIX}/${_qt_plugin_dest}/${_qt_plugin_file}")
  else()
    message(FATAL_ERROR "Qt plugin ${_qt_plugin_name} not found")
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs Qt plugins. Prefix can be the usual CMAKE_INSTALL_PREFIX or
# if you install to a bundle, the app folder
# Fills _qt_plugins_var with paths to be used e.g. with fixup_bundle at BUILD time.
macro(install_qt6_plugin_builddir _qt_plugin_name _qt_plugins_var _targetpath _component)
  get_target_property(_qt_plugin_path "${_qt_plugin_name}" LOCATION)
  if(EXISTS "${_qt_plugin_path}")
    # Resolve symlinks to get the actual file
    get_filename_component(_qt_plugin_path "${_qt_plugin_path}" REALPATH)
    get_filename_component(_qt_plugin_file "${_qt_plugin_path}" NAME)
    get_filename_component(_qt_plugin_type "${_qt_plugin_path}" PATH)
    get_filename_component(_qt_plugin_type "${_qt_plugin_type}" NAME)
    set(_qt_plugin_dest "${_targetpath}/${_qt_plugin_type}")
    install(FILES "${_qt_plugin_path}"
            DESTINATION "${_qt_plugin_dest}"
            COMPONENT ${_component})
    set(${_qt_plugins_var}
        "${${_qt_plugins_var}};${_qt_plugin_dest}/${_qt_plugin_file}")
  else()
    message(FATAL_ERROR "Qt plugin ${_qt_plugin_name} not found")
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs QT6 libraries to CMAKE_INSTALL_PREFIX based on given components
macro(install_qt6_libs _qt_components _targetpath _install_component)
  foreach (_qt_component ${_qt_components})
    get_target_property(_qt_lib_path "Qt6::${_qt_component}" LOCATION)
    if(_qt_lib_path MATCHES "^.*\\/.*${_qt_component}\\.framework\\/.*$")
    ## we could use if Mac but this is more general
      get_filename_component(_qt_lib_path "${_qt_lib_path}" PATH)
      if(EXISTS "${_qt_lib_path}")
      install(DIRECTORY "${_qt_lib_path}"
        DESTINATION "${_targetpath}"
        COMPONENT ${_install_component})
      else()
        message(FATAL_ERROR "Qt lib ${_qt_component} not found at imported location ${_qt_lib_path} for install/package")
      endif()
    else()
      if(EXISTS "${_qt_lib_path}")
        if (UNIX AND "${_qt_lib_path}" MATCHES ".*\\.[0-9]+\\.[0-9]+\\.[0-9]+$")
          string(REGEX REPLACE "\\.[0-9]+\\.[0-9]+$" "" _qt_lib_path_tgt ${_qt_lib_path})
        else()
          set(_qt_lib_path_tgt ${_qt_lib_path})
        endif()
        get_filename_component(_qt_lib_path_tgt "${_qt_lib_path_tgt}" NAME)
        install(FILES "${_qt_lib_path}"
          DESTINATION "${_targetpath}"
          RENAME "${_qt_lib_path_tgt}"
          COMPONENT ${_install_component})
      else()
        message(FATAL_ERROR "Qt lib ${_qt_component} not found at imported location ${_qt_lib_path} for install/package")
      endif()
    endif()
  endforeach(_qt_component)
endmacro()
