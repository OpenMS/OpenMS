# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

# --------------------------------------------------------------------------
# Custom wrapper of Qt's UI tool
# --------------------------------------------------------------------------
macro (qt5_extract_options _qt5_files _qt5_options)
  set(${_qt5_files})
  set(${_qt5_options})
  set(_QT5_DOING_OPTIONS FALSE)
  foreach(_currentArg ${ARGN})
    if ("${_currentArg}" STREQUAL "OPTIONS")
      set(_QT5_DOING_OPTIONS TRUE)
    else ()
      if(_QT5_DOING_OPTIONS)
        list(APPEND ${_qt5_options} "${_currentArg}")
      else()
        list(APPEND ${_qt5_files} "${_currentArg}")
      endif()
    endif ()
  endforeach()
endmacro ()

macro (QT5_WRAP_UI_OWN outfiles )
  qt5_extract_options(ui_files ui_options ${ARGN})

  # create output directory (will not exist for out-of-source builds)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${directory})

  # wrap all files and put them into the
  #  -> ${PROJECT_BINARY_DIR}/${directory}/ui_${outfile}.h
  foreach (_it ${ui_files})
    get_filename_component(outfile ${_it} NAME_WE)
    get_filename_component(infile ${_it} ABSOLUTE)
    set(outfile ${PROJECT_BINARY_DIR}/${directory}/ui_${outfile}.h)
    add_custom_command(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    set(${outfiles} ${${outfiles}} ${outfile})
  endforeach ()
endmacro (QT5_WRAP_UI_OWN)
