# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Chris Bielow, Julianus Pfeuffer $
# --------------------------------------------------------------------------



# Optional targets
# --------------------------------------------------------------------------
if (PYOPENMS)
	set(pyopenms_targets
										COMMAND ${CMAKE_COMMAND} -E echo "    pyopenms           builds pyOpenMS inplace"
										COMMAND ${CMAKE_COMMAND} -E echo "    pyopenms_bdist_egg builds pyOpenMS bdist_egg"
										COMMAND ${CMAKE_COMMAND} -E echo "    pyopenms_bdist     builds pyOpenMS bdist as zip file"
										COMMAND ${CMAKE_COMMAND} -E echo "    pyopenms_rpm       builds pyOpenMS rpm"

	)
else()
	set(pyopenms_targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "    (Disabled) pyopenms targets are not enabled (to enable use -D PYOPENMS=ON)."
										COMMAND ${CMAKE_COMMAND} -E echo ""
    )
endif()

if (OPENMS_COVERAGE)
	set(coverage_target
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS_coverage    generates a coverage report in html format."
										COMMAND ${CMAKE_COMMAND} -E echo "                       Requires the test target to be executed at least once."
	)
else()
	set(coverage_target
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "    (Disabled) OpenMS_coverage reporting target is not enabled (to enable use -D OPENMS_COVERAGE=ON)."
										COMMAND ${CMAKE_COMMAND} -E echo "               Caution: Building with debug and coverage info uses a lot of disk space (>40GB)"
										COMMAND ${CMAKE_COMMAND} -E echo ""
    )
endif()



# --------------------------------------------------------------------------
# Main targets list
if (MSVC)
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    ALL_BUILD       [Visual Studio only] builds the OpenMS library, TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     [NMake only]         builds the OpenMS library, TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    RUN_TESTS       [Visual Studio only] executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            [NMake only]         executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "                    *) make sure they are built using the ALL_BUILD/all target."
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the code snippets of the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen and class documentation, parameters"
										COMMAND ${CMAKE_COMMAND} -E echo "                    documentation, and tutorial PDFs"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_class_only  builds only the doxygen and class documentation"
										COMMAND ${CMAKE_COMMAND} -E echo "                    (faster then doc and very useful when writing"
										COMMAND ${CMAKE_COMMAND} -E echo "                    documentation)."
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the PDF tutorials"
										${pyopenms_targets}
										${coverage_target}
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM)
else()
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     builds the OpenMS library, TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            executes OpenMS and TOPP tests"
										COMMAND ${CMAKE_COMMAND} -E echo "                    make sure they are built using the 'all' target"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the code snippets of the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen and class documentation, parameters"
										COMMAND ${CMAKE_COMMAND} -E echo "                    documentation, and tutorial PDFs"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_class_only  builds only the doxygen and class documentation"
										COMMAND ${CMAKE_COMMAND} -E echo "                    (faster then doc and very useful when writing"
										COMMAND ${CMAKE_COMMAND} -E echo "                    documentation)."
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the PDF tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    help            list all available targets (very verbose)"
										${pyopenms_targets}
										${coverage_target}
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM)
endif()
