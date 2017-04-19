# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
## export a single option indicating if boost static libs should be preferred
option(BOOST_USE_STATIC "Use Boost static libraries." ON)

#------------------------------------------------------------------------------
## Wraps the common find boost code into a single call
## @param .. simply add all required components to the call
## @note This macro will define BOOST_MOC_ARGS that should be added to all moc
##       calls (see https://bugreports.qt-project.org/browse/QTBUG-22829)
macro(find_boost)
  set(Boost_USE_STATIC_LIBS ${BOOST_USE_STATIC})
  set(Boost_USE_MULTITHREADED  ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
  add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (boost tends to guess wrong lib names)
  set(Boost_COMPILER "")

  # help boost finding it's packages
  set(Boost_ADDITIONAL_VERSIONS "1.48.0" "1.49.0" "1.50.0" "1.51.0" "1.52.0" "1.53.0" "1.54.0")

  # 1st attempt does not explicitly requires boost to enable second check (see below)
  find_package(Boost 1.48.0 COMPONENTS ${ARGN})

  set(BOOST_MOC_ARGS "")

  if(Boost_FOUND)
    # see: https://bugreports.qt-project.org/browse/QTBUG-22829
    # Confirmed only on mac os x and leads to problems on win32 and lnx
    # so we handle it for now only on mac os x and boost versions > 1.52
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR ${Boost_MINOR_VERSION} GREATER "52")
      set(BOOST_MOC_ARGS "-DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION")
    endif()
  endif()
endmacro(find_boost)

#------------------------------------------------------------------------------
## Checks if the user supplied package type is valid and aborts if not
## @param package_type The given package type
macro(is_valid_package package_type)
  list(FIND VALID_PACKAGE_TYPES ${package_type} list_pos)
  if( ${list_pos} EQUAL -1 )
  	message(STATUS "The PACKAGE_TYPE ${package_type} is invalid")
  	message(STATUS "Valid PACKAGE_TYPEs are:")
  	foreach( _vpt ${VALID_PACKAGE_TYPES} )
  		message(STATUS " * ${_vpt}")
  	endforeach()
  	message(FATAL_ERROR "Aborting ...")
  endif()
endmacro()
