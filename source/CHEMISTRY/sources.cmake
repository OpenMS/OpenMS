### the directory name
set(directory source/CHEMISTRY)

### list all filenames of the directory here
set(sources_list
AASequence.C
EdwardsLippertIterator.C
EdwardsLippertIteratorTryptic.C
Element.C
ElementDB.C
EmpiricalFormula.C
EnzymaticDigestion.C
IsotopeDistribution.C
ModificationDefinition.C
ModificationDefinitionsSet.C
ModificationsDB.C
ModifierRep.C
PepIterator.C
Residue.C
ResidueDB.C
ResidueModification.C
TheoreticalSpectrumGenerator.C
AdvancedTheoreticalSpectrumGenerator.C
SvmTheoreticalSpectrumGenerator.C
TrypticIterator.C
WeightWrapper.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\CHEMISTRY" FILES ${sources})

