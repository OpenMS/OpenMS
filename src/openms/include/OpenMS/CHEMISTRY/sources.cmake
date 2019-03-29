### the directory name
set(directory include/OpenMS/CHEMISTRY)

### list all header files of the directory here
set(sources_list_h
AAIndex.h
AASequence.h
CrossLinksDB.h
Element.h
ElementDB.h
EmpiricalFormula.h
EnzymaticDigestionLogModel.h
EnzymaticDigestion.h
DigestionEnzyme.h
DigestionEnzymeProtein.h
DigestionEnzymeRNA.h
DigestionEnzymeDB.h
ModificationDefinition.h
ModificationDefinitionsSet.h
ModifiedNASequenceGenerator.h
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
SvmTheoreticalSpectrumGenerator.h
SvmTheoreticalSpectrumGeneratorSet.h
SvmTheoreticalSpectrumGeneratorTrainer.h
TheoreticalSpectrumGenerator.h
TheoreticalSpectrumGeneratorXLMS.h
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
