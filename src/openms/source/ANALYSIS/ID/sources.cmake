### the directory name
set(directory source/ANALYSIS/ID)

### list all filenames of the directory here
set(sources_list
AccurateMassSearchEngine.cpp
AhoCorasickAmbiguous.cpp
AScore.cpp
BasicProteinInferenceAlgorithm.cpp
BayesianProteinInferenceAlgorithm.cpp
ConsensusIDAlgorithm.cpp
ConsensusIDAlgorithmAverage.cpp
ConsensusIDAlgorithmBest.cpp
ConsensusIDAlgorithmIdentity.cpp
ConsensusIDAlgorithmPEPIons.cpp
ConsensusIDAlgorithmPEPMatrix.cpp
ConsensusIDAlgorithmRanks.cpp
ConsensusIDAlgorithmSimilarity.cpp
ConsensusIDAlgorithmWorst.cpp
ConsensusMapMergerAlgorithm.cpp
FalseDiscoveryRate.cpp
HiddenMarkovModel.cpp
IDBoostGraph.cpp
IDConflictResolverAlgorithm.cpp
IDMapper.cpp
IDMergerAlgorithm.cpp
IDRipper.cpp
IDDecoyProbability.cpp
IDScoreGetterSetter.cpp
MessagePasserFactory.cpp
MetaboliteSpectralMatching.cpp
PeptideProteinResolution.cpp
PrecursorPurity.cpp
ProtonDistributionModel.cpp
PeptideIndexing.cpp
PercolatorFeatureSetHelper.cpp
SimpleSearchEngineAlgorithm.cpp
SiriusAdapterAlgorithm.cpp
SiriusMSConverter.cpp
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
