### the directory name
set(directory include/OpenMS/ANALYSIS/ID)

### list all header files of the directory here
set(sources_list_h
AccurateMassSearchEngine.h
AhoCorasickAmbiguous.h
AScore.h
BasicProteinInferenceAlgorithm.h
BayesianProteinInferenceAlgorithm.h
ConsensusIDAlgorithm.h
ConsensusIDAlgorithmAverage.h
ConsensusIDAlgorithmBest.h
ConsensusIDAlgorithmIdentity.h
CometModification.h
ConsensusIDAlgorithmPEPIons.h
ConsensusIDAlgorithmPEPMatrix.h
ConsensusIDAlgorithmRanks.h
ConsensusIDAlgorithmSimilarity.h
ConsensusIDAlgorithmWorst.h
ConsensusMapMergerAlgorithm.h
FalseDiscoveryRate.h
FIAMSDataProcessor.h
FIAMSScheduler.h
FragmentIndex.h
HyperScore.h
IDBoostGraph.h
IDDecoyProbability.h
IDConflictResolverAlgorithm.h
IDMapper.h
IDMergerAlgorithm.h
IDRipper.h
IDScoreGetterSetter.h
IDScoreSwitcherAlgorithm.h
IonIdentityMolecularNetworking.h
MetaboliteSpectralMatching.h
ProSEAlgorithm.h
MorpheusScore.h
NeighborSeq.h
OpenSearchModificationAnalysis.h
PeptideIndexing.h
PeptideProteinResolution.h
Percolator.h
PercolatorFeatureSetHelper.h
PercolatorTypes.h
PrecursorPurity.h
PScore.h
Scores.h
SimpleSearchEngineAlgorithm.h
SiriusExportAlgorithm.h
SiriusMSConverter.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\ID" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
