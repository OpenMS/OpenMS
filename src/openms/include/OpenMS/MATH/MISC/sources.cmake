### the directory name
set(directory include/OpenMS/MATH/MISC)

### list all header files of the directory here
set(sources_list_h
BilinearInterpolation.h
BSpline2d.h
CubicSpline2d.h
LinearInterpolation.h
MathFunctions.h
NonNegativeLeastSquaresSolver.h
MSNumpress.h
RANSAC.h
Spline2d.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\MATH\\MISC" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
