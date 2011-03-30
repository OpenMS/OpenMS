### the directory name
set(directory source/VISUAL/ANNOTATION)

### list all filenames of the directory here
set(sources_list
Annotation1DDistanceItem.C
Annotation1DItem.C
Annotation1DPeakItem.C
Annotation1DTextItem.C
Annotations1DContainer.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\ANNOTATION" FILES ${sources})

