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
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow $
# --------------------------------------------------------------------------


## Check if link libraries have right architecture (.lib can be either a static or link library for a dll)
    ## Motivation: if 32bit libs are accidentally linked against a 64bit OpenMS, weird linker errors (symbol not found, etc) will appear.
    ##             These are really hard to track down, so we check in advance.
function(check_lib_architecture link_libs)

    if (MSVC)
      message(STATUS "  Checking architecture of link libraries ...")
      ## check for 'dumpbin.exe' which comes with every VS to extract architecture
      find_program(PROG_DUMPBIN "dumpbin" DOC "dumpbin allows to extract information from (link) libraries")
      if (NOT PROG_DUMPBIN)
        message(FATAL_ERROR "\n\nCould not find 'dumpbin.exe' (${PROG_DUMPBIN})! Please make sure it's in %PATH%.\n\n")
      endif()
      foreach(lib_file ${${link_libs}})
        if (${lib_file} MATCHES "\\.lib$") ## ends with '.lib'
          ## call 'dumpbin.exe' on each lib
          execute_process(COMMAND "${PROG_DUMPBIN}" "/HEADERS" "${lib_file}"
                          RESULT_VARIABLE DUMP_RESULT
                          OUTPUT_VARIABLE PIPE_OUT_LIB
                          ERROR_VARIABLE PIPE_OUT_LIB
                          TIMEOUT 3) ## three seconds should be enough to get the first lines, which suffice to know what we need; big libs may need much longer!
          ## do not check '${DUMP_RESULT} EQUAL 0' here since a timeout could set it to 'Process terminated due to timeout'
          ## get a shorter version for (faster) regex matching below
          string(SUBSTRING "${PIPE_OUT_LIB}" 0 600 PIPE_OUT_LIB_S)
          ## we expect consistent architectures, output can be
          ## 'Machine      : 8664 (x64)'
          ## or
          ## '8664 machine (x64)'
          ## (similar for x86)
          set(M_64bit "Machine      : 8664 \\(x64\\)|8664 machine \\(x64\\)")
          set(M_32bit "Machine      :  14C \\(x86\\)|14C machine \\(x86\\)")
          if ((OPENMS_64BIT_ARCHITECTURE EQUAL 1) AND NOT (${PIPE_OUT_LIB_S} MATCHES ${M_64bit}))
            message(FATAL_ERROR "\n####\n${PIPE_OUT_LIB_S}\n####\n\nWe detected an external 32-bit link library ('${lib_file}') which cannot be linked to this 64-bit OpenMS! Find or compile a 64-bit version and run CMake again! Evidence see above. Expected substring regex was '${M_64bit}'.")
          elseif((OPENMS_64BIT_ARCHITECTURE EQUAL 0) AND NOT(${PIPE_OUT_LIB_S} MATCHES ${M_32bit}))
            message(FATAL_ERROR "\n####\n${PIPE_OUT_LIB_S}\n####\n\nWe detected an external 64-bit link library ('${lib_file}') which cannot be linked to this 32-bit OpenMS! Find or compile a 32-bit version and run CMake again! Evidence see above. Expected substring was '${M_32bit}'.")
          endif()
        endif()
      endforeach(lib_file)
      message(STATUS "  Checking architecture of link libraries ... done")
    else()
      ## not implemented yet for Linux/MacOSX
    endif()
    
 endfunction()