# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Hannes Röst $
# $Authors: Hannes Röst $
# --------------------------------------------------------------------------

### the directory name
set(directory source/OPENSWATHALGO)

### list all files of the directory here
set(sources_algo_list

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
  DATAACCESS/SwathMap.h
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
