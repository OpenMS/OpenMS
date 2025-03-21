# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_doc")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

#------------------------------------------------------------------------------
# helper macro
#------------------------------------------------------------------------------
macro(make_doxygen_path doxygen_style_list_var doc_path_list_var)
  set(${doxygen_style_list_var} "")
  foreach(_doc_path ${${doc_path_list_var}})
    set(${doxygen_style_list_var} "\"${_doc_path}\" ${${doxygen_style_list_var}}")
  endforeach()
endmacro()

if(ENABLE_DOCS)
  message(STATUS "Documentation is built (-D ENABLE_DOCS=On).")

  #------------------------------------------------------------------------------
  # decide which doc target is part of all
  #------------------------------------------------------------------------------
  if (HAS_XSERVER)
    set(_DOC_ALL "ALL")
    set(_DOC_MINIMAL_ALL "")
    set(_DOC_MINIMAL_EXTRA_MESSAGE )
  else ()
    set(_DOC_ALL "")
    set(_DOC_MINIMAL_ALL "ALL")
    set(_DOC_MINIMAL_EXTRA_MESSAGE
        COMMAND ${CMAKE_COMMAND} -E echo "We build only a minimal doc since OpenMS was configured without HAS_XSERVER=Off."
        COMMAND ${CMAKE_COMMAND} -E echo "If you want the full doc, either execute cmake --build . --target doc manually"
        COMMAND ${CMAKE_COMMAND} -E echo "or rerun cmake with -DHAS_XSERVER=On."
    )
  endif ()

  #------------------------------------------------------------------------------
  # doc programs
  #------------------------------------------------------------------------------
  # TODO switch to target based linking and avoid re-searching Qt.
  # At least reuse the COMPONENT variables
  set(OpenMS_doc_executables)
  include(doxygen/parameters/executables.cmake)

  set(OPENMS_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "doxygen/parameters")

  add_definitions(/DBOOST_ALL_NO_LIB)

  # build doc executables
  foreach(i ${OpenMS_doc_executables})
    add_executable(${i} EXCLUDE_FROM_ALL doxygen/parameters/${i}.cpp)
    
    # let the doc progs decide how to handle a missing OpenMS_GUI library
    if(WITH_GUI)
      set_target_properties(${i} PROPERTIES COMPILE_DEFINITIONS "WITH_GUI=1")
      target_link_libraries(${i} OpenMS OpenMS_GUI)
    else()
      target_link_libraries(${i} OpenMS)
    endif()
  endforeach(i)

  # collection target
  add_custom_target(doc_progs)
  add_dependencies(doc_progs ${OpenMS_doc_executables} TOPP)

  if(WITH_GUI)
    add_dependencies(doc_progs GUI)
  endif()

  #------------------------------------------------------------------------------
  # Initialize variables needed for a proper doxygen configuration
  #------------------------------------------------------------------------------
  set(CF_OPENMS_BIN_PATH ${PROJECT_BINARY_DIR})
  set(CF_OPENMS_SRC_PATH ${OPENMS_HOST_DIRECTORY})

  # make it a doxygen conform list
  make_doxygen_path(CF_OPENMS_DOCUMENTATION_DIRECTORIES OPENMS_DOCUMENTATION_DIRECTORIES)
  make_doxygen_path(CF_OPENMS_DOCUMENTATION_STRIP_INCLUDES OPENMS_DOCUMENTATION_DIRECTORIES)
  make_doxygen_path(CF_STRIP_PATH OPENMS_DOCUMENTATION_DIRECTORIES)

  #------------------------------------------------------------------------------
  # the doc targets
  #------------------------------------------------------------------------------
  # As soon as there is a TeX formula in the html documentation we need LaTeX, dvips and ghostscript.
  # Probably only affects class documentation.

  if (DOXYGEN_FOUND)
  
    # check doxygen for bug with generated latex files
    set(DOXYGEN_MIN_REQUIRED "1.8.13")

    if (DOXYGEN_VERSION VERSION_LESS DOXYGEN_MIN_REQUIRED )
      message(STATUS "Warning: Doxygen ( vers. installed = ${DOXYGEN_VERSION} < ${DOXYGEN_MIN_REQUIRED}) does not correctly process \\include commands in developer's installation instructions. This will make them incomplete. Upgrade Doxygen to fix.")
    endif()
    
    ## if not found, use empty string to avoid `warning: the dot tool could not be found as 'DOXYGEN_DOT_EXECUTABLE-NOTFOUND/dot.exe'` in doxygen-error.log
    if (DOXYGEN_DOT_EXECUTABLE) 
      set(CF_DOXYGEN_DOT_EXECUTABLE "${DOXYGEN_DOT_EXECUTABLE}")
    else()
      set(CF_DOXYGEN_DOT_EXECUTABLE "")
    endif()
      
    
    #------------------------------------------------------------------------------
    # configure doxygen configuration files
    configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_dot.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_dot)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_noclass.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_noclass)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_xml.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_xml)

    #------------------------------------------------------------------------------
    # doc paths (bin path, topp documenter, defaultparamhandler documenter)
    #------------------------------------------------------------------------------
    # Checks for Multiconfiguration generators like MSVC or XCode
    if(NOT CMAKE_CONFIGURATION_TYPES)
      set(_TOPPDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/TOPPDocumenter")
      set(_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/DefaultParamHandlerDocumenter")
      set(_BINARY_PATH "${OPENMS_RUNTIME_OUTPUT_DIRECTORY}")
    else()
      set(_TOPPDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/${CMAKE_CFG_INTDIR}/TOPPDocumenter")
      set(_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/${CMAKE_CFG_INTDIR}/DefaultParamHandlerDocumenter")
      set(_BINARY_PATH "${OPENMS_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}")
    endif()

    file(TO_NATIVE_PATH "${_TOPPDOCUMENTER_EXECUTABLE}" TOPPDOCUMENTER_EXECUTABLE)
    file(TO_NATIVE_PATH "${_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE}" DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE)
    file(TO_NATIVE_PATH "${_BINARY_PATH}" BINARY_PATH)

    # add bin/ and lib/ to all path that potentially contain dependencies
    list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/bin" OUTPUT_VARIABLE DEP_BIN_DIRS)
    list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/lib" OUTPUT_VARIABLE DEP_LIB_DIRS)

    # 
    if (CMAKE_GENERATOR MATCHES "Visual Studio" AND NOT CMAKE_MSVCIDE_RUN_PATH)
      list(APPEND CMAKE_MSVCIDE_RUN_PATH ${DEP_BIN_DIRS} ${DEP_LIB_DIRS})
    else()
      if(WIN32)
        set(pathsep ";")
      else()
        set(pathsep ":")
      endif()
      string(REPLACE ";" "${pathsep}" DEP_BIN_DIRS "${DEP_BIN_DIRS}")
      string(REPLACE ";" "${pathsep}" DEP_LIB_DIRS "${DEP_LIB_DIRS}")   
    endif()


    #------------------------------------------------------------------------------
    # doc_param_internal targets
    if (CMAKE_GENERATOR MATCHES "Visual Studio")
      add_custom_target(doc_param_internal
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating the algorithm parameter and TOPP parameter documentation"
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Note: A functioning OpenMS/TOPP installation and a running X-server (Unix) is required for this step!"
                      COMMAND ${CMAKE_COMMAND} -E echo "      Windows only: OpenMS and Open_GUI.dll's need to be accessible by the Documenter executables in"
                      COMMAND ${CMAKE_COMMAND} -E echo "      doxygen/parameters/. If the automatic CMake post build copy step did not work, try to add them"
                      COMMAND ${CMAKE_COMMAND} -E echo "      to your Path variable."
                      COMMAND ${CMAKE_COMMAND} -E echo "      If this step fails, use the target 'doc_minimal'."
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Building OpenMS parameter docu:"
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E make_directory doxygen/parameters/output/
                      COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Building TOPP docu:"
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${TOPPDOCUMENTER_EXECUTABLE} ${BINARY_PATH}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the parameters documentation"
                      VERBATIM)      
    else()
      add_custom_target(doc_param_internal
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating the algorithm parameter and TOPP parameter documentation"
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Note: A functioning OpenMS/TOPP installation and a running X-server (Unix) is required for this step!"
                      COMMAND ${CMAKE_COMMAND} -E echo "      Windows only: OpenMS and Open_GUI.dll's need to be accessible by the Documenter executables in"
                      COMMAND ${CMAKE_COMMAND} -E echo "      doxygen/parameters/. If the automatic CMake post build copy step did not work, try to add them"
                      COMMAND ${CMAKE_COMMAND} -E echo "      to your Path variable."
                      COMMAND ${CMAKE_COMMAND} -E echo "      If this step fails, use the target 'doc_minimal'."
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Building OpenMS parameter docu:"
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E make_directory doxygen/parameters/output/
                      COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${CMAKE_COMMAND} -E env "PATH=${DEP_BIN_DIRS}${pathsep}${DEP_LIB_DIRS}${pathsep}$ENV{PATH}" ${DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "Building TOPP docu:"
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${CMAKE_COMMAND} -E env "PATH=${DEP_BIN_DIRS}${pathsep}${DEP_LIB_DIRS}${pathsep}$ENV{PATH}" ${TOPPDOCUMENTER_EXECUTABLE} ${BINARY_PATH}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the parameters documentation"
                      VERBATIM)
    endif()

    add_dependencies(doc_param_internal doc_progs)

    #------------------------------------------------------------------------------
    # doc(_html) generation code reused in two independent targets
    set(_DOC_HTML_CODE
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation"
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E remove_directory html
      COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile
      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created."
      COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser."
      COMMAND ${CMAKE_COMMAND} -E echo ""
      COMMENT "Build the doxygen documentation"
    )

    # regular doc target executed to generate full class and tool documentation
    add_custom_target(doc ${_DOC_ALL}
                      ${_DOC_HTML_CODE}
                      DEPENDS doc_param_internal
                      VERBATIM)
                      
    ## A test which inspects the doxygen-error.log and fails if Doxygen reports errors; also forwards the logfile content to CDash for inspection
    add_executable(Doxygen_Warning_Checker doxygen/test/Doxygen_Warning_Checker.cpp)
    ## this is a compile time dependency, but better than none. Currently, there is no way to have a test depend on a target ('doc') before its run.
    #add_dependencies(Doxygen_Warning_Checker doc) ## disable, since CI seems to re-configure CMake after 'doc' and chokes on it
    add_test(Doxygen_Warning_test ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Doxygen_Warning_Checker "${PROJECT_BINARY_DIR}/doxygen/doxygen-error.log" ${DOXYGEN_VERSION})
    set_tests_properties(Doxygen_Warning_test PROPERTIES ATTACHED_FILES_ON_FAIL "${PROJECT_BINARY_DIR}/doxygen/doxygen-error.log")
    
    # doc_html_only to generate only the doxygen documentation (e.g., when writing new documentation)
    add_custom_target(doc_class_only
                      ${_DOC_HTML_CODE}
                      COMMAND ${CMAKE_COMMAND} -E echo "NOTE: The algorithm/TOPP parameter documentation was not generated/updated."
                      COMMAND ${CMAKE_COMMAND} -E echo "      You will only see updates to the class documentation/doxygen files."
                      COMMAND ${CMAKE_COMMAND} -E echo "      To build the full documentation execute the doc target."
                      VERBATIM)

    #------------------------------------------------------------------------------
    # doc_internal target
    add_custom_target(doc_xml
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating XML documentation"
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E remove_directory xml_output
                      COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_xml
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "The XML documentation has been successfully created."
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the doxygen documentation"
                      VERBATIM)

    #------------------------------------------------------------------------------
    # doc_noclass target
    add_custom_target(doc_noclass
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class documentation"
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E remove_directory html
                      COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_noclass
                      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created."
                      COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser."
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the doxygen documentation"
                      VERBATIM)
    add_dependencies(doc_noclass doc_param_internal)

    #------------------------------------------------------------------------------
    # doc_minimal target
    add_custom_target(doc_minimal ${_DOC_MINIMAL_ALL}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class/TOPP documentation"
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E remove_directory html
                      COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_noclass
                      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created."
                      COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser."
                      ${_DOC_MINIMAL_EXTRA_MESSAGE}
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the doxygen documentation"
                      VERBATIM)
    #------------------------------------------------------------------------------
    # add virtual dependency of doc_minimal on TOPP
    # this is not necessary but defers the generation of doc_minimal to a later
    # stage to avoid confusion if doc_minimal is built first
    add_dependencies(doc_minimal TOPP)

    if (DOXYGEN_DOT_FOUND OR DOXYGEN_DOT_EXECUTABLE)
      #------------------------------------------------------------------------------
      # doc_dot target
      add_custom_target(doc_dot
                        COMMAND ${CMAKE_COMMAND} -E echo ""
                        COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        COMMAND ${CMAKE_COMMAND} -E echo "Creating DOT html documentation"
                        COMMAND ${CMAKE_COMMAND} -E echo ""
                        COMMAND ${CMAKE_COMMAND} -E remove_directory html-dot
                        COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_dot
                        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                        COMMAND ${CMAKE_COMMAND} -E echo ""
                        COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created."
                        COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser."
                        COMMAND ${CMAKE_COMMAND} -E echo ""
                        COMMENT "Build the doxygen documentation"
                        VERBATIM)
      add_dependencies(doc_dot doc_param_internal)
    else()
      message(STATUS "DOT not found. Disabling target 'doc_dot'!")
    endif()

    #------------------------------------------------------------------------------
    # Install documentation
    install(DIRECTORY ${PROJECT_BINARY_DIR}/html
      DESTINATION ${INSTALL_DOC_DIR}
      COMPONENT doc
      OPTIONAL
      FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                            GROUP_READ
                            WORLD_READ
      DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                            GROUP_EXECUTE GROUP_READ
                            WORLD_EXECUTE WORLD_READ
      REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
      REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
      REGEX ".*\\.map\$" EXCLUDE ## Exclude temporary map files
      REGEX ".*\\.md5\$" EXCLUDE ## Exclude temporary md5 files
    )
    # TODO would be cool to somehow have documentation optional
    #install_file(${PROJECT_BINARY_DIR}/index.html ${INSTALL_DOC_DIR} doc)

  else()
    message(STATUS "Doxygen not found. Disabling all documentation targets!")
    message(STATUS "Note that no documentation will be installed alongside OpenMS.")
  endif()

  #------------------------------------------------------------------------------
  # Examples / Tutorials
  #------------------------------------------------------------------------------
  option(BUILD_EXAMPLES "Compile OpenMS code examples (C++ tutorials)" ON)
  message(STATUS "Build C++ tutorials: ${BUILD_EXAMPLES} (-DBUILD_EXAMPLES=ON/OFF);")
  if(BUILD_EXAMPLES AND "${PACKAGE_TYPE}" STREQUAL "none")
    add_subdirectory(code_examples)
    message(STATUS " --> Tutorials in '${PROJECT_SOURCE_DIR}/code_examples' will be built!")
  else()
    message(STATUS " --> Tutorials in '${PROJECT_SOURCE_DIR}/code_examples' will NOT be built! (since PACKAGE_TYPE=${PACKAGE_TYPE}, but must be 'none')")
  endif()


else()
  message(STATUS "Documentation is not built (-D ENABLE_DOCS=Off).")
endif()
