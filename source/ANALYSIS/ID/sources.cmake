### the directory name
set(directory source/ANALYSIS/ID)

### list all filenames of the directory here
set(sources_list
AScore.C
ConsensusID.C
FalseDiscoveryRate.C
HiddenMarkovModel.C
IDMapper.C
PILISIdentification.C
PILISModel.C
PILISModelGenerator.C
PILISNeutralLossModel.C
PILISScoring.C
PILISCrossValidation.C
ProtonDistributionModel.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\ID" FILES ${sources})

