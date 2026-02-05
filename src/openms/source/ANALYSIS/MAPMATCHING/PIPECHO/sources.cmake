### the directory name
set(directory source/ANALYSIS/MAPMATCHING/PIPECHO)

### list all filenames of the directory here
set(sources_list_h
  GridWithStorage.h
  MzDiff.h
  PeakTypes.h
  Run.h
  RunStatistics.h
  Score.h
  Util.h
  Window.h
)

set(sources_list
  PeakTypes.cpp
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
set(OpenMS_sources_h ${OpenMS_private_sources_h} ${sources_h})

### source group definition
#source_group("Source Files\\ANALYSIS\\MAPMATCHING\\PIPECHO" FILES ${sources})
