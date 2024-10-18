# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

# a collection of wrapper for install functions that allows easier usage
# throughout the OpenMS build system

set(OPENMS_EXPORT_SET "OpenMSTargets")

#------------------------------------------------------------------------------
# Installs the library lib_target_name and all its headers set via
# set_target_properties(lib_target_name PROPERTIES PUBLIC_HEADER ${headers})
#
# @param lib_target_name The target name of the library that should be installed
macro(install_library lib_target_name)
    install(TARGETS ${lib_target_name}
      RUNTIME_DEPENDENCY_SET OPENMS_DEPS
      EXPORT ${OPENMS_EXPORT_SET}
      LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      RUNTIME DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      )
endmacro()

#------------------------------------------------------------------------------
# Installs the given headers.
#
# @param header_list List of headers to install
macro(install_headers header_list component)
  foreach(_header ${header_list})
    set(_relative_header_path)

    get_filename_component(_target_path ${_header} PATH)
    if ("${_target_path}" MATCHES "^${PROJECT_BINARY_DIR}.*")
      # is generated bin header
      string(REPLACE "${PROJECT_BINARY_DIR}/include/OpenMS" "" _relative_header_path "${_target_path}")
    else()
      # is source header -> strip include/OpenMS
      string(REPLACE "include/OpenMS" "" _relative_header_path "${_target_path}")
    endif()

    # install the header
    install(FILES ${_header}
            # note the missing slash, we need this for file directly located in
            # include/OpenMS (e.g., config.h)
            DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS${_relative_header_path}
            COMPONENT ${component}_headers)
  endforeach()
endmacro()

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
# Installs the exported target information
macro(install_export_targets )
    install(EXPORT ${OPENMS_EXPORT_SET}
            DESTINATION ${INSTALL_LIB_DIR}/cmake/OpenMS
            COMPONENT cmake)
endmacro()

#------------------------------------------------------------------------------
# Installs Thirdparty folders with executables
macro(install_thirdparty_folder foldername)
  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/${foldername})
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/${foldername}
            DESTINATION           ${INSTALL_SHARE_DIR}/THIRDPARTY
            COMPONENT             ${foldername}
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
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

