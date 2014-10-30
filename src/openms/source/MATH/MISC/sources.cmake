### the directory name
set(directory source/MATH/MISC)

### list all filenames of the directory here
set(sources_list
BilinearInterpolation.cpp
BSpline2d.cpp
CubicSpline2d.cpp
LinearInterpolation.cpp
MathFunctions.cpp
NonNegativeLeastSquaresSolver.cpp
MSNumpress.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\MATH\\MISC" FILES ${sources})
