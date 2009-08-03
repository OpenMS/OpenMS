### the directory name
set(directory include/OpenMS/ANALYSIS/DENOVO)

### list all header files of the directory here
set(sources_list_h
DeNovoIonScoring.h
DeNovoAlgorithm.h
DeNovoPostScoring.h
DeNovoIdentification.h
MassDecomposition.h
MassDecompositionAlgorithm.h
CompNovoIdentificationBase.h
CompNovoIdentification.h
CompNovoIonScoringCID.h
CompNovoIdentificationCID.h
CompNovoIonScoringBase.h
CompNovoIonScoring.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\DENOVO" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

