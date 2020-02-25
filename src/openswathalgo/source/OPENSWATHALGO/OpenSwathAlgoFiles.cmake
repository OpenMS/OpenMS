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
# $Maintainer: Hannes Röst $
# $Authors: Hannes Röst $
# --------------------------------------------------------------------------

### the directory name
set(directory source/OPENSWATHALGO)

### list all files of the directory here
set(sources_algo_list
  ALGO/MRMScoring.cpp
  ALGO/Scoring.cpp
  ALGO/StatsHelpers.cpp
)

set(sources_dataaccess_list
  DATAACCESS/DataFrameWriter.cpp
  DATAACCESS/ISpectrumAccess.cpp
  DATAACCESS/MockObjects.cpp
  DATAACCESS/SpectrumHelpers.cpp
  DATAACCESS/TransitionHelper.cpp
  DATAACCESS/Transitions.cpp
)

### add path to the filenames
set(sources_algo)
foreach(i ${sources_algo_list})
	list(APPEND sources_algo ${directory}/${i})
endforeach(i)

set(sources_dataaccess)
foreach(i ${sources_dataaccess_list})
	list(APPEND sources_dataaccess ${directory}/${i})
endforeach(i)

set(OpenSwathAlgoFiles)
list(APPEND OpenSwathAlgoFiles ${sources_algo})
list(APPEND OpenSwathAlgoFiles ${sources_dataaccess})

source_group("Source Files\\ANALYSIS\\OPENSWATH\\OPENSWATHALGO\\DATAACESS" FILES ${sources_dataaccess})
source_group("Source Files\\ANALYSIS\\OPENSWATH\\OPENSWATHALGO\\ALGO" FILES ${sources_algo})

## the header directory
set(header_directory include/OpenMS/OPENSWATHALGO)

## add groups for headers
set(header_algo_list
  ALGO/MRMScoring.h
  ALGO/Scoring.h
  ALGO/StatsHelpers.h
)
set(header_dataaccess_list
  DATAACCESS/DataFrameWriter.h
  DATAACCESS/DataStructures.h
  DATAACCESS/ISpectrumAccess.h
  DATAACCESS/ITransition.h
  DATAACCESS/MockObjects.h
  DATAACCESS/SpectrumHelpers.h
  DATAACCESS/TransitionExperiment.h
  DATAACCESS/TransitionHelper.h
  DATAACCESS/Transitions.h
)

### add path to the filenames
set(header_algo)
foreach(i ${header_algo_list})
	list(APPEND header_algo ${header_directory}/${i})
endforeach(i)

set(header_dataaccess)
foreach(i ${header_dataaccess_list})
	list(APPEND header_dataaccess ${header_directory}/${i})
endforeach(i)

list(APPEND OpenSwathAlgoFiles ${header_algo})
list(APPEND OpenSwathAlgoFiles ${header_dataaccess})

# define list of headers related to openswathalgo needed for
# installation and export
set(OpenSwathAlgoHeaders ${header_algo} ${header_dataaccess})

source_group("Header Files\\ANALYSIS\\OPENSWATH\\OPENSWATHALGO\\DATAACESS" FILES ${header_dataaccess})
source_group("Header Files\\ANALYSIS\\OPENSWATH\\OPENSWATHALGO\\ALGO" FILES ${header_algo})
