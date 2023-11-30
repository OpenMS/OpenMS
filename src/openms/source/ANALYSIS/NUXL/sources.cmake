### the directory name
set(directory source/ANALYSIS/NUXL)

### list all filenames of the directory here
set(sources_list
NuXLAnnotateAndLocate.cpp
NuXLFDR.cpp
NuXLFeatureAugmentation.cpp
NuXLFragmentAdductDefinition.cpp
NuXLFragmentAnnotationHelper.cpp
NuXLFragmentIonGenerator.cpp
NuXLMarkerIonExtractor.cpp
NuXLModificationsGenerator.cpp
NuXLParameterParsing.cpp
NuXLReport.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\NUXL" FILES ${sources})

