### the directory name
set(directory include/OpenMS/CHEMISTRY)

### list all header files of the directory here
set(sources_list_h
AAIndex.h
AASequence.h
AdductInfo.h
CrossLinksDB.h
DecoyGenerator.h
Element.h
ElementDB.h
EmpiricalFormula.h
EnzymaticDigestion.h
DigestionEnzyme.h
DigestionEnzymeProtein.h
DigestionEnzymeRNA.h
DigestionEnzymeDB.h
ModificationDefinition.h
ModificationDefinitionsSet.h
ModifiedNASequenceGenerator.h
ModifiedPeptideGenerator.h
ModificationsDB.h
NASequence.h
NucleicAcidSpectrumGenerator.h
ProteaseDB.h
ProteaseDigestion.h
Residue.h
ResidueDB.h
ResidueModification.h
RNaseDB.h
RNaseDigestion.h
Ribonucleotide.h
RibonucleotideDB.h
SimpleTSGXLMS.h
SpectrumAnnotator.h
Tagger.h
TheoreticalSpectrumGenerator.h
TheoreticalSpectrumGeneratorXLMS.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CHEMISTRY" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
