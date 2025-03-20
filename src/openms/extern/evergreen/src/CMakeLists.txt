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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS-evergreen")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# export include directory
add_library(Evergreen INTERFACE)
target_include_directories(Evergreen SYSTEM INTERFACE 
                            "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
                            "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen>")

# Export it under OpenMSTargets.cmake instead of a separate EvergreenTargets file
#install(EXPORT EvergreenTargets DESTINATION lib/cmake/evergreen)
#install(TARGETS Evergreen EXPORT EvergreenTargets PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen COMPONENT library)

install_library(Evergreen)
openms_register_export_target(Evergreen)

get_filename_component(EVERGREEN_HEADERS "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(EVERGREEN_FOLDER "${EVERGREEN_HEADERS}" PATH)
file(APPEND ${TPLICENSEFILE} "\nevergreen:\n")
cat(${EVERGREEN_FOLDER}/LICENSE ${TPLICENSEFILE})

# Install headers
install(DIRECTORY ${EVERGREEN_HEADERS}/
  DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen
  COMPONENT thirdparty_headers
  FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                        GROUP_READ
                        WORLD_READ
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ
                        WORLD_EXECUTE WORLD_READ
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
  PATTERN "demos/*" EXCLUDE
  )
