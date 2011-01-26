### the directory name
set(directory include/OpenMS/CHEMISTRY)

### list all header files of the directory here
set(sources_list_h
AASequence.h
EdwardsLippertIterator.h
EdwardsLippertIteratorTryptic.h
Element.h
ElementDB.h
EmpiricalFormula.h
EnzymaticDigestion.h
IsotopeDistribution.h
ModificationDefinition.h
ModificationDefinitionsSet.h
ModificationsDB.h
ModifierRep.h
PepIterator.h
Residue.h
ResidueDB.h
ResidueModification.h
TheoreticalSpectrumGenerator.h
SvmTheoreticalSpectrumGenerator.h
TrypticIterator.h
AAIndex.h
WeightWrapper.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CHEMISTRY" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

