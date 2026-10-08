### the directory name
set(directory source/ANALYSIS/MAPMATCHING/PIPECHO)

### list all filenames of the directory here
set(sources_list_h
  FeatureTypes.h
  GridWithStorage.h
  Impl.h
  MzDiff.h
  TransferFDRModel.h
  Run.h
  RunStatistics.h
  Score.h
  Util.h
  Window.h
)

set(sources_list
  FeatureTypes.cpp
  Impl.cpp
  TransferFDRModel.cpp
  Run.cpp
  RunStatistics.cpp
  Util.cpp
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
  list(APPEND sources_h ${directory}/${i})
endforeach(i)

set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})
# Private implementation headers that live next to the .cpp under source/.  They
# are registered on OpenMS_private_headers (declared in includes.cmake), which is
# merged into OpenMS_sources for IDE listing only -- never installed, exported,
# or wrapped.  Public headers must not include them.
set(OpenMS_private_headers ${OpenMS_private_headers} ${sources_h})

### source group definition
source_group("Source Files\\ANALYSIS\\MAPMATCHING\\PIPECHO" FILES ${sources})
source_group("Header Files\\OpenMS\\ANALYSIS\\MAPMATCHING\\PIPECHO" FILES ${sources_h})
