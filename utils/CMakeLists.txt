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

project("OpenMS_UTILS")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

# --------------------------------------------------------------------------
# OpenMS' UTILS tools

# add OpenMS includes 
include_directories(${OPENMS_INCLUDE_DIRECTORIES})

# add all the tools
set(UTILS_executables)
include(executables.cmake)
foreach(i ${UTILS_executables})
	add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)

# collection target 
add_custom_target(UTILS)
add_dependencies(UTILS ${UTILS_executables})

## some regular UTILS need the GUI lib
foreach(i ${UTILS_executables_with_GUIlib})
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
  target_include_directories(${i} PUBLIC ${OPENMS_GUI_INCLUDE_DIRECTORIES})  
endforeach(i)

## export the list of UTILS tools into CACHE
set(UTILS_TOOLS ${UTILS_executables}
    CACHE INTERNAL "OpenMS' UTILS tools" FORCE)