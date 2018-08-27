### the directory name
set(directory source/FILTERING/CALIBRATION)

### list all filenames of the directory here
set(sources_list
InternalCalibration.cpp
MZTrafoModel.cpp
PrecursorCorrection.cpp
TOFCalibration.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FILTERING\\CALIBRATION" FILES ${sources})

