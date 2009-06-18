### the directory name
set(directory source/SIMULATION)

### list all header files of the directory here
set(sources_list
DetectabilitySimulation.C
DigestSimulation.C
ElutionModel.C
IonizationSimulation.C
IsotopeModelGeneral.C
MSSim.C
MixtureModel.C
PTMSimulation.C
RTSimulation.C
RawMSSignalSimulation.C
RawTandemMSSignalSimulation.C
sources.cmake
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\SIMULATION" FILES ${sources})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)
