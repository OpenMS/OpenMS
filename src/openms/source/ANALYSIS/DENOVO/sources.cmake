### the directory name
set(directory source/ANALYSIS/DENOVO)

### list all filenames of the directory here
set(sources_list
DeNovoIonScoring.cpp
DeNovoAlgorithm.cpp
DeNovoPostScoring.cpp
DeNovoIdentification.cpp
CompNovoIdentificationBase.cpp
CompNovoIdentificationCID.cpp
CompNovoIonScoring.cpp
CompNovoIdentification.cpp
CompNovoIonScoringBase.cpp
CompNovoIonScoringCID.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\DENOVO" FILES ${sources})

