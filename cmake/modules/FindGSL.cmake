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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

# Based on the FindGLPK module

SET(GSL_ROOT_DIR "" CACHE PATH "GSL root directory")

SET(GSL_REGKEY "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Gsl;InstallPath]")
GET_FILENAME_COMPONENT(GSL_ROOT_PATH ${GSL_REGKEY} ABSOLUTE)

FIND_PATH(GSL_INCLUDE_DIR
  gsl_version.h
  PATHS ${GSL_REGKEY}/include
  HINTS ${GSL_ROOT_DIR}/include
  PATH_SUFFIXES gsl
)
FIND_LIBRARY(GSL_LIBRARY_RELEASE
  gsl
  PATHS ${GSL_REGKEY}/lib
  HINTS ${GSL_ROOT_DIR}/lib
)

FIND_LIBRARY(GSL_LIBRARY_DEBUG
  gsl_d
  PATHS ${GSL_REGKEY}/lib
  HINTS ${GSL_ROOT_DIR}/lib
)

FIND_LIBRARY(GSL_CBLAS_LIBRARY_RELEASE
  gslcblas
	cblas
  PATHS ${GSL_REGKEY}/lib
  HINTS ${GSL_ROOT_DIR}/lib
)

FIND_LIBRARY(GSL_CBLAS_LIBRARY_DEBUG
  gslcblas_d
	cblas_d
  PATHS ${GSL_REGKEY}/lib
  HINTS ${GSL_ROOT_DIR}/lib
)

include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
SELECT_LIBRARY_CONFIGURATIONS(GSL)
SELECT_LIBRARY_CONFIGURATIONS(GSL_CBLAS)

IF(GSL_INCLUDE_DIR AND GSL_LIBRARY AND GSL_CBLAS_LIBRARY)
  FILE(READ ${GSL_INCLUDE_DIR}/gsl_version.h GSL_GSL_VERSION_H)

  string(REGEX MATCH "define[ ]+GSL_VERSION[ ]+\"[0-9]+\\.[0-9]+\"" GSL_VERSION_LINE "${GSL_GSL_VERSION_H}")
  string(REGEX REPLACE "define[ ]+GSL_VERSION[ ]+\"([0-9.]+)\"" "\\1" GSL_VERSION_FULL "${GSL_VERSION_LINE}")
  string(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\1" GSL_VERSION_MAJOR "${GSL_VERSION_FULL}")
  string(REGEX REPLACE "([0-9]+)\\.([0-9]+)" "\\2" GSL_VERSION_MINOR "${GSL_VERSION_FULL}")  

  SET(GSL_VERSION_STRING "${GSL_VERSION_MAJOR}.${GSL_VERSION_MINOR}")

  IF(GSL_FIND_VERSION)
    IF(GSL_FIND_VERSION_COUNT GREATER 2)
      MESSAGE(SEND_ERROR "unexpected version string")
    ENDIF(GSL_FIND_VERSION_COUNT GREATER 2)

    MATH(EXPR GSL_REQUESTED_VERSION "${GSL_FIND_VERSION_MAJOR}*100 + ${GSL_FIND_VERSION_MINOR}")
    MATH(EXPR GSL_FOUND_VERSION "${GSL_VERSION_MAJOR}*100 + ${GSL_VERSION_MINOR}")

    IF(GSL_FOUND_VERSION LESS GSL_REQUESTED_VERSION)
      SET(GSL_PROPER_VERSION_FOUND FALSE)
    ELSE(GSL_FOUND_VERSION LESS GSL_REQUESTED_VERSION)
      SET(GSL_PROPER_VERSION_FOUND TRUE)
    ENDIF(GSL_FOUND_VERSION LESS GSL_REQUESTED_VERSION)
  ELSE(GSL_FIND_VERSION)
    SET(GSL_PROPER_VERSION_FOUND TRUE)
  ENDIF(GSL_FIND_VERSION)
ENDIF(GSL_INCLUDE_DIR AND GSL_LIBRARY AND GSL_CBLAS_LIBRARY)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL 
                                  DEFAULT_MSG 
                                  GSL_INCLUDE_DIR
                                  GSL_LIBRARY 
                                  GSL_CBLAS_LIBRARY 
                                  GSL_PROPER_VERSION_FOUND)

IF(GSL_FOUND)
  SET(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
  SET(GSL_LIBRARIES "${GSL_LIBRARY};${GSL_CBLAS_LIBRARY}")
  SET(GSL_BIN_DIR ${GSL_ROOT_PATH}/bin)
ENDIF(GSL_FOUND)

MARK_AS_ADVANCED(GSL_LIBRARY 
                 GSL_CBLAS_LIBRARY 
                 GSL_INCLUDE_DIR 
                 GSL_BIN_DIR)
