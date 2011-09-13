
##### targets list #####
if (MSVC)
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    ALL_BUILD       [Visual Studio only] builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     [NMake only]         builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    UTILS           builds the UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    RUN_TESTS       [Visual Studio only] executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            [NMake only]         executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "                    *) make sure they are built using the 'test_build' target (see below)"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_exec  executes the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen documentation and tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the pdf tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools and UTILS have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo "The class tests have their own project in ./source/TEST (project test_build)."
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM
		)
else()
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    UTILS           builds the UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test_build      builds the OpenMS tests"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            executes OpenMS and TOPP tests"
										COMMAND ${CMAKE_COMMAND} -E echo "                    make sure they are built using the 'test_build' target"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_exec  executes the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen documentation and tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the pdf tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    help            list all available targets (very long)"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools and UTILS have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo "The class tests have their own project in ./source/TEST."
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM
		)
endif()

##### Message after OpenMS has been built #####
if (MSVC)
	## copy OpenMS.dll to test executables dir
	get_target_property(WIN32_DLLLOCATION OpenMS LOCATION)
	get_filename_component(WIN32_DLLPATH ${WIN32_DLLLOCATION} PATH)

	## copy OpenMS.dll to test executables dir "$(TargetFileName)" is a placeholder filled by VS at runtime
	file(TO_NATIVE_PATH "${WIN32_DLLPATH}/$(TargetFileName)" DLL_SOURCE)
	file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}/source/TEST/bin/$(ConfigurationName)/$(TargetFileName)" DLL_TARGET)
	# create target path if not exists
	add_custom_command(TARGET OpenMS
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS library has been built."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "You should now build the TOPP tools and tests."
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Execute the 'targets' project to see prominent targets!"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/source/TEST/bin/$(ConfigurationName)"
										COMMAND copy ${DLL_SOURCE} ${DLL_TARGET} /Y
                    COMMENT "Library build. Copying OpenMS(d).dll to test binary dir"
                    VERBATIM)

	## copy OpenMS_GUI.dll
	add_custom_command(TARGET OpenMS_GUI
                    POST_BUILD
					COMMAND copy ${DLL_SOURCE} ${DLL_TARGET} /Y
                    COMMENT "GUI Library build. Copying OpenMS_GUI(d).dll to test binary dir"
                    VERBATIM)

					
else()
  add_custom_command(TARGET OpenMS
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS library has been built."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "You should now build the TOPP tools and tests."
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Make sure to add the OpenMS and contrib lib/ path"
                    COMMAND ${CMAKE_COMMAND} -E echo "to your LD_LIBRARY_PATH environment variable."
                    COMMAND ${CMAKE_COMMAND} -E echo "Otherwise the tests and TOPP tools will not work."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "For a full list of targets execute:"
                    COMMAND ${CMAKE_COMMAND} -E echo "make targets"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    VERBATIM)
endif()

##### Message after TOPP has been built #####
if (MSVC)
  add_custom_command(TARGET TOPP
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The TOPP tools have been built and installed to the bin/ folder."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Execute the 'targets' project to see prominent targets!"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    VERBATIM)
else()
  add_custom_command(TARGET TOPP
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The TOPP tools have been built and installed to the bin/ folder."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Make sure to add the OpenMS and contrib lib/ path"
                    COMMAND ${CMAKE_COMMAND} -E echo "to your LD_LIBRARY_PATH environment variable."
                    COMMAND ${CMAKE_COMMAND} -E echo "Otherwise the tests and TOPP tools will not work."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "For a full list of make targets execute:"
                    COMMAND ${CMAKE_COMMAND} -E echo "make targets"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    VERBATIM)
endif()


##### Messages at the end of cmake #####
MESSAGE(STATUS "")
MESSAGE(STATUS "-----------------------------------------------------------------")
MESSAGE(STATUS "")
MESSAGE(STATUS "You have successfully configured OpenMS and TOPP.")
MESSAGE(STATUS "")
if (MSVC)
  MESSAGE(STATUS "Execute the 'targets' project to see prominent targets!")
else()
  MESSAGE(STATUS "For a full list of make targets execute:")
  MESSAGE(STATUS "'make targets'")
endif()
MESSAGE(STATUS "")
MESSAGE(STATUS "-----------------------------------------------------------------")
MESSAGE(STATUS "")
