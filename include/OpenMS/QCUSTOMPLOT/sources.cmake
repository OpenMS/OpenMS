### the directory name
set(directory include/OpenMS/QCUSTOMPLOT)

### list all MOC filenames of the directory here
set(sources_list
qcustomplot.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\QCUSTOMPLOT" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
qcustomplot.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\QCUSTOMPLOT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

