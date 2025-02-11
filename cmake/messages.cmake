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
