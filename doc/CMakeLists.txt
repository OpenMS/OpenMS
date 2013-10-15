# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#------------------------------------------------------------------------------
# helper macro
#------------------------------------------------------------------------------
macro(make_doxygen_path doxygen_style_list_var doc_path_list_var)
  set(${doxygen_style_list_var} "")
  foreach(_doc_path ${${doc_path_list_var}})
    set(${doxygen_style_list_var} "\"${_doc_path}\" ${${doxygen_style_list_var}}")
  endforeach()
endmacro()

#------------------------------------------------------------------------------
# doc programs
#------------------------------------------------------------------------------
set(OpenMS_doc_executables)
include(doxygen/parameters/executables.cmake)

set(OPENMS_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "doxygen/parameters")

# add OpenMS includes 
include_directories(${OPENMS_INCLUDE_DIRECTORIES})

# build doc executables
foreach(i ${OpenMS_doc_executables})
	add_executable(${i} EXCLUDE_FROM_ALL doxygen/parameters/${i}.C)
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
endforeach(i)

# collection target
add_custom_target(doc_progs)
add_dependencies(doc_progs ${OpenMS_doc_executables} TOPP UTILS GUI)

#------------------------------------------------------------------------------
# Intialize variables needed for a proper doxygen configuration
#------------------------------------------------------------------------------
set(CF_OPENMS_BIN_PATH ${PROJECT_BINARY_DIR})
set(CF_OPENMS_SRC_PATH ${OPENMS_HOST_DIRECTORY})

set(OPENMS_DOCUMENTATION_DIRECTORIES)
list(APPEND OPENMS_DOCUMENTATION_DIRECTORIES "${OPENMS_HOST_DIRECTORY}/openms/include/OpenMS")
list(APPEND OPENMS_DOCUMENTATION_DIRECTORIES "${OPENMS_HOST_DIRECTORY}/openms/source/APPLICATIONS/TOPP/")
list(APPEND OPENMS_DOCUMENTATION_DIRECTORIES "${OPENMS_HOST_DIRECTORY}/openms/source/APPLICATIONS/UTILS/")
list(APPEND OPENMS_DOCUMENTATION_DIRECTORIES "${OPENMS_HOST_DIRECTORY}/openms/source/VISUAL/APPLICATIONS/")

# make it a doxygen conform list
make_doxygen_path(CF_OPENMS_DOCUMENTATION_DIRECTORIES OPENMS_DOCUMENTATION_DIRECTORIES)

# strip include path
set(DOC_OPENMS_INCLUDE_DIRECTORIES)
list(APPEND DOC_OPENMS_INCLUDE_DIRECTORIES "${OPENMS_HOST_DIRECTORY}/openms/include/")

# make it a doxygen conform list
make_doxygen_path(CF_OPENMS_DOCUMENTATION_STRIP_INCLUDES DOC_OPENMS_INCLUDE_DIRECTORIES)

#------------------------------------------------------------------------------
# find the necessary packages
find_package(Doxygen)
find_package(LATEX)

#------------------------------------------------------------------------------
# the doc targets
#------------------------------------------------------------------------------
if (DOXYGEN_FOUND)
  #------------------------------------------------------------------------------
  # configure doxygen configuration files
  configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile)
  configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_dot.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_dot)
  configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_noclass.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_noclass)
  configure_file(${PROJECT_SOURCE_DIR}/doxygen/Doxyfile_xml.in ${PROJECT_BINARY_DIR}/doxygen/Doxyfile_xml)
  configure_file(${PROJECT_SOURCE_DIR}/OpenMS_tutorial/Doxyfile.in ${PROJECT_BINARY_DIR}/OpenMS_tutorial/Doxyfile)
  configure_file(${PROJECT_SOURCE_DIR}/TOPP_tutorial/Doxyfile.in ${PROJECT_BINARY_DIR}/TOPP_tutorial/Doxyfile)

  #------------------------------------------------------------------------------
  # create refman files for PDF tutorials
  configure_file(${PROJECT_SOURCE_DIR}/OpenMS_tutorial/refman_overwrite.tex.in ${PROJECT_BINARY_DIR}/OpenMS_tutorial/refman_overwrite.tex)
  configure_file(${PROJECT_SOURCE_DIR}/TOPP_tutorial/refman_overwrite.tex.in ${PROJECT_BINARY_DIR}/TOPP_tutorial/refman_overwrite.tex)

  #------------------------------------------------------------------------------
  # doc paths (bin path, topp documenter, defaultparamhandler documenter)
  #------------------------------------------------------------------------------  
	if(NOT MSVC)
		set(_TOPPDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/TOPPDocumenter")
		set(_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/DefaultParamHandlerDocumenter")
		set(_BINARY_PATH "${OPENMS_RUNTIME_OUTPUT_DIRECTORY}")
	else()
		set(_TOPPDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/$(ConfigurationName)/TOPPDocumenter")
		set(_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE "${PROJECT_BINARY_DIR}/doxygen/parameters/$(ConfigurationName)/DefaultParamHandlerDocumenter")		
		set(_BINARY_PATH "${OPENMS_RUNTIME_OUTPUT_DIRECTORY}/$(ConfigurationName)")
	endif()
	
	file(TO_NATIVE_PATH "${_TOPPDOCUMENTER_EXECUTABLE}" TOPPDOCUMENTER_EXECUTABLE)
	file(TO_NATIVE_PATH "${_DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE}" DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE)
	file(TO_NATIVE_PATH "${_BINARY_PATH}" BINARY_PATH)
	
  #------------------------------------------------------------------------------
  # doc_param_internal targets
  add_custom_target(doc_param_internal
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating the algorithm parameter and TOPP parameter documentation";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "Note: A functioning OpenMS/TOPP installation and a running X-server (Linux) is required for this step!";
                    COMMAND ${CMAKE_COMMAND} -E echo "      Windows only: OpenMS and Open_GUI.dll's need to be in your PATH!";
                    COMMAND ${CMAKE_COMMAND} -E echo "      If this step fails, use the target 'doc_minimal'.";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "Building OpenMS parameter docu:"
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    COMMAND ${CMAKE_COMMAND} -E make_directory doxygen/parameters/output/
                    COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${DEFAULTPARAMHANDLERDOCUMENTER_EXECUTABLE}
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Building TOPP/UTILS docu:"
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~"
                    COMMAND ${CMAKE_COMMAND} -E chdir doxygen/parameters/ ${TOPPDOCUMENTER_EXECUTABLE} ${BINARY_PATH}
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "Build the parameters documentation"
                    VERBATIM)
	add_dependencies(doc_param_internal doc_progs)

  #------------------------------------------------------------------------------
  # doc target
  add_custom_target(doc
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory html
                    COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doxygen/common/style_ini.css html/style_ini.css
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
                    COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "Build the doxygen documentation"
                    VERBATIM)
  add_dependencies(doc doc_param_internal)

  #------------------------------------------------------------------------------
  # doc_internal target
  add_custom_target(doc_xml
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating XML documentation";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory xml_output
                    COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_xml
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "The XML documentation has been successfully created.";
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "Build the doxygen documentation"
                    VERBATIM)

  #------------------------------------------------------------------------------
  # doc_noclass target
  add_custom_target(doc_noclass
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class documentation";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory html
                    COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_noclass
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doxygen/common/style_ini.css html/style_ini.css
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
                    COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "Build the doxygen documentation"
                    VERBATIM)
  add_dependencies(doc_noclass doc_param_internal)

  #------------------------------------------------------------------------------
  # doc_minimal target
  add_custom_target(doc_minimal
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class/TOPP/UTILS documentation";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory html
                    COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_noclass
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doxygen/common/style_ini.css html/style_ini.css
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
                    COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "Build the doxygen documentation"
                    VERBATIM)

  if (DOXYGEN_DOT_FOUND OR DOXYGEN_DOT_EXECUTABLE)
    #------------------------------------------------------------------------------
    # doc_dot target
    add_custom_target(doc_dot
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                      COMMAND ${CMAKE_COMMAND} -E echo "Creating DOT html documentation";
                      COMMAND ${CMAKE_COMMAND} -E echo "";
                      COMMAND ${CMAKE_COMMAND} -E remove_directory html-dot
                      COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} "${DOXYGEN_EXECUTABLE}" doxygen/Doxyfile_dot
                      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/index.html index.html
                      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doxygen/common/style_ini.css html/style_ini.css
                      COMMAND ${CMAKE_COMMAND} -E echo "";
                      COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                      COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
                      COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
                      COMMAND ${CMAKE_COMMAND} -E echo ""
                      COMMENT "Build the doxygen documentation"
                      VERBATIM)
    add_dependencies(doc_dot doc_param_internal)
  else()
    message(STATUS "DOT not found. Disabling target 'doc_dot'!")
  endif()

else()
  message(STATUS "Doxygen not found. Disabling all documentation targets!")
endif()

#------------------------------------------------------------------------------
# PDF Tutorials
#------------------------------------------------------------------------------
if (DOXYGEN_FOUND AND LATEX_COMPILER AND DVIPS_CONVERTER)
  set(DOC_TUTORIALS_ACTIVE TRUE)
  
  set(DOXYGEN_START_BUGGY "1.6.3")
  set(DOXYGEN_END_BUGGY "1.7.2")
  exec_program(${DOXYGEN_EXECUTABLE}
    ARGS "--version"
    OUTPUT_VARIABLE DOXYGEN_VERSION)

  if (DOXYGEN_VERSION STRGREATER DOXYGEN_START_BUGGY AND DOXYGEN_VERSION STRLESS DOXYGEN_END_BUGGY )
    message(ERROR "Warning, DoxygenBug ( 1.6.? < vers. installed < 1.7.3 ) disguises generated tex inputfiles and files will not be recognized")
  endif ()
  
  #------------------------------------------------------------------------------
  # doc_tutorials target
  add_custom_target(doc_tutorials
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating OpenMS pdf tutorial";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory OpenMS_tutorial/latex_output
                    COMMAND ${CMAKE_COMMAND} -E chdir OpenMS_tutorial/ "${DOXYGEN_EXECUTABLE}" Doxyfile
                    COMMAND ${CMAKE_COMMAND} -E copy OpenMS_tutorial/refman_overwrite.tex OpenMS_tutorial/latex_output/refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir OpenMS_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir OpenMS_tutorial/latex_output/  makeindex refman.idx
                    COMMAND ${CMAKE_COMMAND} -E chdir OpenMS_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir OpenMS_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E copy OpenMS_tutorial/latex_output/refman.pdf OpenMS_tutorial.pdf
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS tutorial in PDF format has been successfully created:";
                    COMMAND ${CMAKE_COMMAND} -E echo "doc/OpenMS_tutorial.pdf";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
                    COMMAND ${CMAKE_COMMAND} -E echo "Creating TOPP/TOPPView pdf tutorial";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E remove_directory TOPP_tutorial/latex_output
                    COMMAND ${CMAKE_COMMAND} -E chdir TOPP_tutorial/ "${DOXYGEN_EXECUTABLE}" Doxyfile
                    COMMAND ${CMAKE_COMMAND} -E copy TOPP_tutorial/refman_overwrite.tex TOPP_tutorial/latex_output/refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir TOPP_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir TOPP_tutorial/latex_output/  makeindex refman.idx
                    COMMAND ${CMAKE_COMMAND} -E chdir TOPP_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E chdir TOPP_tutorial/latex_output/  pdflatex refman.tex
                    COMMAND ${CMAKE_COMMAND} -E copy TOPP_tutorial/latex_output/refman.pdf TOPP_tutorial.pdf
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS tutorial in PDF format has been successfully created:";
                    COMMAND ${CMAKE_COMMAND} -E echo "doc/TOPP_tutorial.pdf";
                    COMMAND ${CMAKE_COMMAND} -E echo "";
                    COMMENT "Build the OpenMS/TOPP pdf tutorial"
                    VERBATIM)
else()
  set(DOC_TUTORIALS_ACTIVE FALSE)
  message(STATUS "Doxygen or Latex missing. Disabling 'doc_tutorials' target!")
endif()

#------------------------------------------------------------------------------
# Examples / Tutorials
#------------------------------------------------------------------------------
option(BUILD_EXAMPLES "Compile OpenMS code examples" ON)
if(BUILD_EXAMPLES AND "${PACKAGE_TYPE}" STREQUAL "none")
  add_subdirectory(code_examples)
endif(BUILD_EXAMPLES AND "${PACKAGE_TYPE}" STREQUAL "none")
