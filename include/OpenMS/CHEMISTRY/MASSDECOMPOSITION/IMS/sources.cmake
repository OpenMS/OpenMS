### Files from IMSLIB
### the directory name
set(directory include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS)

### list all header files of the directory here
set(sources_list_h
Alphabet.h
IntegerMassDecomposer.h
MassDecomposer.h
RealMassDecomposer.h
Weights.h
Element.h
IMSIsotopeDistribution.h
AlphabetParser.h
compose_f_gx_hy_t.h
compose_f_gx_t.h
AlphabetTextParser.h
)


### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CHEMISTRY\\MASSDECOMPOSITION\\IMS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

