### the directory name
set(directory source/ANALYSIS/MAPMATCHING/PIPECHO)

### list all filenames of the directory here
set(sources_list_h
  FeatureTypes.h
  GridWithStorage.h
  Impl.h
  MzDiff.h
  Pep.h
  Run.h
  RunStatistics.h
  Score.h
  Util.h
  Window.h
)

set(sources_list
  FeatureTypes.cpp
  Impl.cpp
  Pep.cpp
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
# These are private implementation headers that live next to the .cpp under
# source/.  includes.cmake resets OpenMS_sources_h after this (source) phase, so
# register them on the master source list (IDE listing only; not compiled).
set(OpenMS_sources ${OpenMS_sources} ${sources_h})

### source group definition
#source_group("Source Files\\ANALYSIS\\MAPMATCHING\\PIPECHO" FILES ${sources})
