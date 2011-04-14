#######################################################################
#create Doxyfiles for html documentation
set(CF_OPENMS_BIN_PATH ${PROJECT_BINARY_DIR})
set(CF_OPENMS_SRC_PATH ${PROJECT_SOURCE_DIR})

find_package(Doxygen)
find_package(LATEX)

if (DOXYGEN_FOUND)
	configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile)
	configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile_internal.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile_internal)
	configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile_dot.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile_dot)
	configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile_noclass.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile_noclass)
	configure_file(${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile_xml.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile_xml)
	configure_file(${PROJECT_SOURCE_DIR}/doc/OpenMS_tutorial/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial/Doxyfile)
	configure_file(${PROJECT_SOURCE_DIR}/doc/TOPP_tutorial/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/TOPP_tutorial/Doxyfile)

	#######################################################################
	#create refman files for PDF tutorials
	configure_file(${PROJECT_SOURCE_DIR}/doc/OpenMS_tutorial/refman_overwrite.tex.in ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial/refman_overwrite.tex)
	configure_file(${PROJECT_SOURCE_DIR}/doc/TOPP_tutorial/refman_overwrite.tex.in ${PROJECT_BINARY_DIR}/doc/TOPP_tutorial/refman_overwrite.tex)


	#######################################################################
	##doc_param_internal target
	add_custom_target(doc_param_internal
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating the algorithm parameter and TOPP parameter documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "Note: A functioning OpenMS/TOPP installation and a running X-server (Linux) is required for this step!";
										COMMAND ${CMAKE_COMMAND} -E echo "      If this step fails, use the target 'doc_minimal'.";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "Building OpenMS parameter docu:"
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
										COMMAND ${CMAKE_COMMAND} -E make_directory doc/doxygen/parameters/output/
										COMMAND ${CMAKE_COMMAND} -E chdir doc/doxygen/parameters/ ./DefaultParamHandlerDocumenter
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Building TOPP/UTILS docu:"
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~"
										COMMAND ${CMAKE_COMMAND} -E chdir doc/doxygen/parameters/ ./TOPPDocumenter
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the parameters documentation"
										VERBATIM)
	if (MSVC_IDE)
		##copy required executables:
		add_custom_target(doc_prepare
											COMMAND ${CMAKE_COMMAND} -E copy  ${PROJECT_BINARY_DIR}/doc/doxygen/parameters/$(OutDir)/DefaultParamHandlerDocumenter.exe ${PROJECT_BINARY_DIR}/doc/doxygen/parameters/DefaultParamHandlerDocumenter.exe
											COMMAND ${CMAKE_COMMAND} -E copy  ${PROJECT_BINARY_DIR}/doc/doxygen/parameters/$(OutDir)/TOPPDocumenter.exe ${PROJECT_BINARY_DIR}/doc/doxygen/parameters/TOPPDocumenter.exe
											VERBATIM)
		add_dependencies(doc_prepare doc_progs)
		add_dependencies(doc_param_internal doc_prepare)
	else()
		add_dependencies(doc_param_internal doc_progs)
	endif()

	#######################################################################
	## doc target
	add_custom_target(doc
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html
										COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile
										COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doc/index.html doc/index.html
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
										COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the doxygen documentation"
										VERBATIM)
	add_dependencies(doc doc_param_internal)

	#######################################################################
	## doc_internal target
	add_custom_target(doc_internal
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating internal html documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html
										COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile_internal
										COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doc/index.html doc/index.html
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
										COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the doxygen documentation"
										VERBATIM)
	add_dependencies(doc_internal doc_param_internal)

	#######################################################################
	## doc_internal target
	add_custom_target(doc_xml
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating XML documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/xml_output
										COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile_xml
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "The XML documentation has been successfully created.";
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the doxygen documentation"
										VERBATIM)

	#######################################################################
	## doc_noclass target
	add_custom_target(doc_noclass
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html
										COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile_noclass
										COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doc/index.html doc/index.html
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
										COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the doxygen documentation"
										VERBATIM)
	add_dependencies(doc_noclass doc_param_internal)

	#######################################################################
	## doc_minimal target
	add_custom_target(doc_minimal
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation without class/TOPP/UTILS documentation";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html
										COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile_noclass
										COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doc/index.html doc/index.html
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
										COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "Build the doxygen documentation"
										VERBATIM)

	if (DOXYGEN_DOT_FOUND OR DOXYGEN_DOT_EXECUTABLE)
		#######################################################################
		## doc_dot target
		add_custom_target(doc_dot
											COMMAND ${CMAKE_COMMAND} -E echo ""
											COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
											COMMAND ${CMAKE_COMMAND} -E echo "Creating DOT html documentation";
											COMMAND ${CMAKE_COMMAND} -E echo "";
											COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html-dot
											COMMAND ${CMAKE_COMMAND} -E chdir doc doxygen doxygen/Doxyfile_dot
											COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/doc/index.html doc/index.html
											COMMAND ${CMAKE_COMMAND} -E echo "";
											COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
											COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created.";
											COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/index.html' in a web browser.";
											COMMAND ${CMAKE_COMMAND} -E echo ""
											COMMENT "Build the doxygen documentation"
											VERBATIM)
		add_dependencies(doc_dot doc_param_internal)
	else()
		Message(STATUS "DOT not found. Disabling target 'doc_dot'!")
	endif()

else()
	Message(STATUS "Doxygen not found. Disabling all documentation targets!")
endif()

if (DOXYGEN_FOUND AND LATEX_COMPILER AND DVIPS_CONVERTER)
	#######################################################################
	# doc_tutorials target
	set(DOC_TUTORIALS_ACTIVE TRUE)

  set(DOXYGEN_START_BUGGY "1.6.3")
  set(DOXYGEN_END_BUGGY "1.7.2")
  EXEC_PROGRAM(${DOXYGEN_EXECUTABLE}
    ARGS "--version"
    OUTPUT_VARIABLE DOXYGEN_VERSION)
if (DOXYGEN_VERSION STRGREATER DOXYGEN_START_BUGGY AND DOXYGEN_VERSION STRLESS DOXYGEN_END_BUGGY )
  MESSAGE(ERROR "Warning, DoxygenBug ( 1.6.? < vers. installed < 1.7.3 ) disguises generated tex inputfiles and files will not be recognized")
endif ()

	add_custom_target(doc_tutorials
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating OpenMS pdf tutorial";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/OpenMS_tutorial/latex_output
										COMMAND ${CMAKE_COMMAND} -E chdir doc/OpenMS_tutorial/ doxygen Doxyfile
										COMMAND ${CMAKE_COMMAND} -E copy doc/OpenMS_tutorial/refman_overwrite.tex doc/OpenMS_tutorial/latex_output/refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/OpenMS_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/OpenMS_tutorial/latex_output/	makeindex refman.idx
										COMMAND ${CMAKE_COMMAND} -E chdir doc/OpenMS_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/OpenMS_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E copy doc/OpenMS_tutorial/latex_output/refman.pdf doc/OpenMS_tutorial.pdf
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS tutorial in PDF format has been successfully created:";
										COMMAND ${CMAKE_COMMAND} -E echo "doc/OpenMS_tutorial.pdf";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
										COMMAND ${CMAKE_COMMAND} -E echo "Creating TOPP/TOPPView pdf tutorial";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E remove_directory doc/TOPP_tutorial/latex_output
										COMMAND ${CMAKE_COMMAND} -E chdir doc/TOPP_tutorial/ doxygen Doxyfile
										COMMAND ${CMAKE_COMMAND} -E copy doc/TOPP_tutorial/refman_overwrite.tex doc/TOPP_tutorial/latex_output/refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/TOPP_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/TOPP_tutorial/latex_output/	makeindex refman.idx
										COMMAND ${CMAKE_COMMAND} -E chdir doc/TOPP_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E chdir doc/TOPP_tutorial/latex_output/	pdflatex refman.tex
										COMMAND ${CMAKE_COMMAND} -E copy doc/TOPP_tutorial/latex_output/refman.pdf doc/TOPP_tutorial.pdf
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS tutorial in PDF format has been successfully created:";
										COMMAND ${CMAKE_COMMAND} -E echo "doc/TOPP_tutorial.pdf";
										COMMAND ${CMAKE_COMMAND} -E echo "";
										COMMENT "Build the OpenMS/TOPP pdf tutorial"
										VERBATIM)
else()
	set(DOC_TUTORIALS_ACTIVE FALSE)
	Message(STATUS "Doxygen or Latex missing. Disabling 'doc_tutorials' target!")
endif()

