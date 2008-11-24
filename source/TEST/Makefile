# -*- mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------

include ../config_defs.mak

DIRECTORY=TEST

ADD_CXXFLAGS=-D_FILE_OFFSET_BITS=64

# ClassTest_test should be the first one, because it tests the testing ...

#CONCEPT

EXECUTABLES_CONCEPT = \
	ClassTest_test\
	Exception_Base_test\
	FactoryBase_test\
	Factory_test\
	FactoryProduct_test\
	PreciseTime_test\
	SingletonRegistry_test\
	TimeStamp_test\
	VersionInfo_test\
	FuzzyStringComparator_test

#DATASTRUCTURES

EXECUTABLES_DATASTRUCTURES = \
	String_test\
	StringList_test\
	IntList_test\
	DoubleList_test\
	Date_test\
	Map_test\
	DPosition_test\
	DIntervalBase_test\
	DRange_test\
	DBoundingBox_test\
	DataValue_test\
	DateTime_test\
	RangeManager_test\
	Matrix_test\
	ConvexHull2D_test\
	DefaultParamHandler_test\
	SuffixArray_test\
	SuffixArraySeqan_test\
	SuffixArrayTrypticSeqan_test\
	SuffixArrayTrypticCompressed_test\
	SuffixArrayPeptideFinder_test\
	SparseVector_test\
	DistanceMatrix_test

#METADATA

EXECUTABLES_METADATA = \
	MetaInfoRegistry_test\
	MetaInfo_test\
	MetaInfoInterface_test\
	IonSource_test\
	IonDetector_test\
	MassAnalyzer_test\
	Instrument_test\
	ProteinHit_test\
	PeptideHit_test\
	SampleTreatment_test\
	Digestion_test\
	Modification_test\
	Tagging_test\
	Sample_test\
	Acquisition_test\
	ProteinIdentification_test\
	PeptideIdentification_test\
	AcquisitionInfo_test\
	Precursor_test\
	ContactPerson_test\
	SourceFile_test\
	Software_test\
	MetaInfoDescription_test\
	DataProcessing_test\
	InstrumentSettings_test\
	SpectrumSettings_test\
	Gradient_test\
	HPLC_test\
	ExperimentalSettings_test\
	DocumentIdentifier_test

#SYSTEM

EXECUTABLES_SYSTEM = \
	StopWatch_test\
	File_test\
	ProcessResource_test\
	ExternalAllocator_test\
	ExternalAllocatorUnique_test\
	MemoryMap_test

#KERNEL

EXECUTABLES_KERNEL =\
	PeakIndex_test\
	Peak1D_test\
	Peak2D_test\
	DPeak_test\
	RichPeak1D_test\
	RichPeak2D_test\
	DRichPeak_test\
	ConstRefVector_test\
	DSpectrum_test\
	ComparatorUtils_test\
	StandardTypes_test\
	Feature_test\
	FeatureMap_test\
	RangeUtils_test\
	MSSpectrum_test\
	MSExperiment_test\
	ConsensusMap_test\
	ConsensusFeature_test\
	AreaIterator_test

#VISUAL


EXECUTABLES_VISUAL =\
	MultiGradient_test\
	AxisTickCalculator_test

#FORMAT

ifeq ($(ANDIMS_SUPPORT),true)
	ANDI_EXECUTABLES_FORMAT=\
		ANDIFile_test
endif

EXECUTABLES_FORMAT = \
	XMLValidator_test\
	TextFile_test\
	XMLFile_test\
	$(ANDI_EXECUTABLES_FORMAT)\
	LibSVMEncoder_test\
	DTAFile_test\
	DTA2DFile_test\
	MascotInfile_test\
	MascotOutfile_test\
	Base64_test\
	UniqueIdGenerator_test\
	PersistentObject_test\
	FASTAFile_test\
	FileHandler_test\
	MzXMLFile_test\
	MzMLFile_test\
	MzDataFile_test\
	Param_test\
	FeatureXMLFile_test\
	MascotXMLFile_test\
	PeakFileOptions_test\
	PeakTypeEstimator_test\
	InspectInfile_test\
	InspectOutfile_test\
	SequestInfile_test\
	SequestOutfile_test\
	PepNovoInfile_test\
	PepNovoOutfile_test\
	PTMXMLFile_test\
	ConsensusXMLFile_test\
	DBConnection_test\
	DBAdapter_test\
	OMSSAXMLFile_test\
	IdXMLFile_test\
	BigString_test\
	XTandemXMLFile_test\
	MSPFile_test\
	OMSSACSVFile_test\
	XTandemInfile_test\
	UnimodXMLFile_test\
	PepXMLFile_test\
	ControlledVocabulary_test\
	TransformationXMLFile_test\
	SemanticValidator_test\
	MS2File_test

#MATH

EXECUTABLES_MATH = \
	LinearInterpolation_test\
	AveragePosition_test\
	BasicStatistics_test\
	AsymmetricStatistics_test\
	StatisticFunctions_test\
	ROCCurve_test\
	Histogram_test\
	MathFunctions_test\
	LinearRegression_test\
	BilinearInterpolation_test\
	GammaDistributionFitter_test\
	GaussFitter_test\
	NonNegativeLeastSquaresSolver_test


#FILTERING

EXECUTABLES_FILTERING = \
	LinearResampler_test\
	TopHatFilter_test\
	MorphFilter_test\
	SignalToNoiseEstimator_test\
	SignalToNoiseEstimatorMedian_test\
	SignalToNoiseEstimatorMeanIterative_test\
	GaussFilter_test\
	PreprocessingFunctor_test\
	NLargest_test\
	Scaler_test\
	Normalizer_test\
	ParentPeakMower_test\
	SqrtMower_test\
	ThresholdMower_test\
	WindowMower_test\
	FilterFunctor_test\
	ComplementFilter_test\
	GoodDiffFilter_test\
	IntensityBalanceFilter_test\
	IsotopeDiffFilter_test\
	NeutralLossDiffFilter_test\
	PeakDensityFilter_test\
	TICFilter_test\
	PeakMarker_test\
	ComplementMarker_test\
	IsotopeMarker_test\
	NeutralLossMarker_test\
	MarkerMower_test\
	SavitzkyGolayFilter_test\
	InternalCalibration_test\
	IDFilter_test\
	TOFCalibration_test\
	ParentFilter_test\
	BernNorm_test\
	DataFilters_test


#COMPARISON

EXECUTABLES_COMPARISON=\
	SpectrumCheapDPCorr_test\
	HierarchicalClustering_test\
	PeakSpectrumCompareFunctor_test\
	SpectrumPrecursorComparator_test\
	ZhangSimilarityScore_test\
	SpectrumAlignment_test\
	SpectrumAlignmentScore_test\
	SteinScottImproveScore_test\
	CompareFouriertransform_test\
	BinnedSpectrum_test\
	BinnedSpectrumCompareFunctor_test\
	BinnedSumAgreeingIntensities_test\
	BinnedSpectralContrastAngle_test\
	BinnedSharedPeakCount_test\
	ClusterFunctor_test\
	SingleLinkage_test\
	CompleteLinkage_test\
	AverageLinkage_test\
	PeakAlignment_test\
	ClusterHierarchical_test\
	ClusterAnalyzer_test\
	EuclideanSimilarity_test

#CHEMISTRY

EXECUTABLES_CHEMISTRY = \
	IsotopeDistribution_test\
	Element_test\
	ElementDB_test\
	EmpiricalFormula_test\
	ResidueModification_test\
	Residue_test\
	ResidueDB_test\
	AASequence_test\
	EnzymaticDigestion_test\
	TheoreticalSpectrumGenerator_test\
	ModifierRep_test\
	FastaIterator_test\
	FastaIteratorIntern_test\
	PepIterator_test\
	EdwardsLippertIterator_test\
	EdwardsLippertIteratorTryptic_test\
	TrypticIterator_test\
	ModificationsDB_test\
	ModificationDefinition_test\
	ModificationDefinitionsSet_test

#ANALYSIS

EXECUTABLES_ANALYSIS_QUANTITATION = \
	ItraqChannelExtractor_test \
	ItraqQuantifier_test
				
EXECUTABLES_ANALYSIS_SVM = \
	SVMWrapper_test

EXECUTABLES_ANALYSIS_MAPMATCHING=\
	TransformationDescription_test\
	FeatureHandle_test\
	MapAlignmentAlgorithm_test\
	MapAlignmentAlgorithmSpectrumAlignment_test\
	MapAlignmentAlgorithmPoseClustering_test\
	LabeledPairFinder_test\
	BaseGroupFinder_test\
	DelaunayPairFinder_test\
	SimplePairFinder_test\
	BaseSuperimposer_test\
	PoseClusteringAffineSuperimposer_test\
	PoseClusteringShiftSuperimposer_test\
	FeatureGroupingAlgorithm_test\
	FeatureGroupingAlgorithmLabeled_test\
	FeatureGroupingAlgorithmUnlabeled_test

EXECUTABLES_ANALYSIS_DECHARGING = \
	FeatureDecharger_test

EXECUTABLES_ANALYSIS_ID =\
	IDMapper_test\
	ConsensusID_test\
	ProteinInference_test\
	ProtonDistributionModel_test\
	PILISScoring_test\
	HiddenMarkovModel_test\
	IDDecoyProbability_test\
	FalseDiscoveryRate_test

# APPLICATIONS

EXECUTABLES_APPLICATIONS = \
	TOPPBase_test

# TRANSFORMATIONS

EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER = \
	ModelDescription_test\
	BaseModel_test\
	EmgModel_test\
	LmaGaussModel_test\
	ExtendedIsotopeModel_test\
	InterpolationModel_test\
	GaussModel_test\
	BiGaussModel_test\
	IsotopeModel_test\
	LmaIsotopeModel_test\
	ProductModel_test\
	FeatureFinder_test\
	FeaFiModule_test\
	FeatureFinderAlgorithm_test\
	FeatureFinderAlgorithmPicked_test\
	FeatureFinderAlgorithmSimplest_test\
	FeatureFinderAlgorithmSimple_test\
	FeatureFinderAlgorithmWavelet_test\
	FeatureFinderAlgorithmIsotopeWavelet_test\
	FeatureFinderAlgorithmWatershed_test\
	SimpleSeeder_test\
	SimpleExtender_test\
  GaussFitter1D_test\
  BiGaussFitter1D_test\
  IsotopeFitter1D_test\
  EmgFitter1D_test\
  LmaGaussFitter1D_test\
	LmaIsotopeFitter1D_test\
	ExtendedIsotopeFitter1D_test\
	Fitter1D_test\
  LevMarqFitter1D_test\
  MaxLikeliFitter1D_test\
	ModelFitter_test\
	IsotopeWavelet_test\
	IsotopeWaveletTransform_test

EXECUTABLES_TRANSFORMATIONS_PEAKPICKER = \
	ContinuousWaveletTransform_test\
	PeakShape_test\
	ContinuousWaveletTransformNumIntegration_test\
	OptimizePick_test\
	PeakPickerCWT_test\
	OptimizePeakDeconvolution_test\
	TwoDOptimization_test

EXECUTABLES = \
	$(EXECUTABLES_COMPARISON)\
	$(EXECUTABLES_MATH)\
	$(EXECUTABLES_FILTERING)\
	$(EXECUTABLES_CONCEPT)\
	$(EXECUTABLES_METADATA)\
	$(EXECUTABLES_SYSTEM)\
	$(EXECUTABLES_KERNEL)\
	$(EXECUTABLES_FORMAT)\
	$(EXECUTABLES_DATASTRUCTURES)\
	$(EXECUTABLES_VISUAL)\
	$(EXECUTABLES_CHEMISTRY)\
	$(EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER)\
	$(EXECUTABLES_TRANSFORMATIONS_PEAKPICKER)\
	$(EXECUTABLES_ANALYSIS_ID)\
	$(EXECUTABLES_ANALYSIS_DECHARGING)\
	$(EXECUTABLES_ANALYSIS_SVM)\
	$(EXECUTABLES_ANALYSIS_QUANTITATION)\
	$(EXECUTABLES_ANALYSIS_MAPMATCHING)\
	$(EXECUTABLES_APPLICATIONS)


# default target
all:	$(EXECUTABLES)

# required for clean and depend
CPP_SOURCES=$(EXECUTABLES:=.C)
OBJECTS=$(EXECUTABLES:=.o)
VALGRINDS=$(EXECUTABLES:=.vgr)
#VALGRINDS=$(EXECUTABLES_CHEMISTRY:=.vgr)

include ../common.mak
-include .Dependencies

# remove the library specific flags
# and replace them by the application specific
# flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS)

# do not optimize tests (reduces compile time drastically)
CPP_MODE_FLAGS= $(CPP_MODE_FLAGS_NO_OPTIMIZATION) -O0

.o:
	$(CXX) $(CXXFLAGS) $*.o $(LDFLAGS) -o $@ $(LIBS)

%%.C:
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(OPENMS_CXXFLAGS) $(OPENMS_INCLUDES) $*.C $(LDFLAGS) -o $@ $(LIBS)

.o.vgr:
	@$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(OPENMS_CXXFLAGS) $(OPENMS_INCLUDES) $*.o $(LDFLAGS) -o $@ $(LIBS)
	$(VALGRIND) $(VALGRIND_OPTS) $* > $*.vgr 2>&1 || true

dummy:

# concept_test should be the first one, because it tests the testing ...
test: info concept_test system_test math_test datastructures_test metadata_test kernel_test format_test filtering_test visual_test comparison_test chemistry_test transformations_featurefinder_test transformations_peakpicker_test analysis_id_test analysis_decharging_test analysis_svm_test analysis_quantitation_test analysis_mapmatching_test applications_test
	@./runtests -success

info:
	@./runtests -info

concept_test:	$(EXECUTABLES_CONCEPT) dummy
	@./runtests -section "CONCEPT" $(TEST_OPTIONS) $(EXECUTABLES_CONCEPT)

math_test:	$(EXECUTABLES_MATH) dummy
	@./runtests -section "MATH" $(TEST_OPTIONS) $(EXECUTABLES_MATH)

filtering_test:	$(EXECUTABLES_FILTERING) dummy
	@./runtests -section "FILTERING" $(TEST_OPTIONS) $(EXECUTABLES_FILTERING)

comparison_test:	$(EXECUTABLES_COMPARISON) dummy
	@./runtests -section "COMPARISON" $(TEST_OPTIONS) $(EXECUTABLES_COMPARISON)

kernel_test:	$(EXECUTABLES_KERNEL) dummy
	@./runtests -section "KERNEL" $(TEST_OPTIONS) $(EXECUTABLES_KERNEL)

datastructures_test:	$(EXECUTABLES_DATASTRUCTURES) dummy
	@./runtests -section "DATASTRUCTURES" $(TEST_OPTIONS) $(EXECUTABLES_DATASTRUCTURES)

format_test:	$(EXECUTABLES_FORMAT) dummy
	@./runtests -section "FORMAT" $(TEST_OPTIONS) $(EXECUTABLES_FORMAT)

system_test:	$(EXECUTABLES_SYSTEM) dummy
	@./runtests -section "SYSTEM" $(TEST_OPTIONS) $(EXECUTABLES_SYSTEM)

metadata_test:	$(EXECUTABLES_METADATA) dummy
	@./runtests -section "METADATA" $(TEST_OPTIONS) $(EXECUTABLES_METADATA)

visual_test:	$(EXECUTABLES_VISUAL) dummy
	@./runtests -section "VISUAL" $(TEST_OPTIONS) $(EXECUTABLES_VISUAL)

chemistry_test:	$(EXECUTABLES_CHEMISTRY) dummy
	@./runtests -section "CHEMISTRY" $(TEST_OPTIONS) $(EXECUTABLES_CHEMISTRY)

transformations_featurefinder_test:	$(EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER) dummy
	@./runtests -section "TRANSFORMATIONS_FEATUREFINDER" $(TEST_OPTIONS) $(EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER)

transformations_peakpicker_test:	$(EXECUTABLES_TRANSFORMATIONS_PEAKPICKER) dummy
	@./runtests -section "TRANSFORMATIONS_PEAKPICKER" $(TEST_OPTIONS) $(EXECUTABLES_TRANSFORMATIONS_PEAKPICKER)

analysis_id_test:	$(EXECUTABLES_ANALYSIS_ID) dummy
	@./runtests -section "ANALYSIS_ID" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_ID)

analysis_decharging_test:	$(EXECUTABLES_ANALYSIS_DECHARGING) dummy
	@./runtests -section "ANALYSIS_DECHARGING" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_DECHARGING)

analysis_svm_test:	$(EXECUTABLES_ANALYSIS_SVM) dummy
	@./runtests -section "ANALYSIS_SVM" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_SVM)

analysis_quantitation_test:	$(EXECUTABLES_ANALYSIS_QUANTITATION) dummy
	@./runtests -section "ANALYSIS_QUANTITATION" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_QUANTITATION)

analysis_mapmatching_test:	$(EXECUTABLES_ANALYSIS_MAPMATCHING) dummy
	@./runtests -section "ANALYSIS_MAPMATCHING" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_MAPMATCHING)

applications_test:	$(EXECUTABLES_APPLICATIONS) dummy
	@./runtests -section "APPLICATIONS" $(TEST_OPTIONS) $(EXECUTABLES_APPLICATIONS)

clean:	default_clean
	@$(RM) -rf -- *.o *~ *.log ii_files Templates.DB	SunWS_cache *.TMP *.tmp core a.out region* cut* feature* model* *.vgr *.output $(EXECUTABLES)
	@cd TOPP && make clean

#creates a .vgr-file with valgrind check output for each test
valgrind: $(VALGRINDS)
