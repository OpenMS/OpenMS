# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../config_defs.mak

DIRECTORY=TEST

ADD_CXXFLAGS=-D_FILE_OFFSET_BITS=64

# ClassTest_test should be the first one, because it tests the testing ...

#CONCEPT

EXECUTABLES_CONCEPT = \
	ClassTest_test\
	Exception_Base_test\
	PreciseTime_test\
	TimeStamp_test\
	VersionInfo_test\
	FactoryProduct_test\
	Factory_test

#DATASTRUCTURES

EXECUTABLES_DATASTRUCTURES = \
	String_test\
	Date_test\
	HashMap_test\
	DPosition_test\
	DIntervalBase_test\
	DRange_test\
	DBoundingBox_test\
	DataValue_test\
	SparseVector_test\
	DateTime_test\
	RangeManager_test\
	Matrix_test\
	ConvexHull2D_test\
	DefaultParamHandler_test

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
	ProcessingMethod_test\
	InstrumentSettings_test\
	SpectrumSettings_test\
	Gradient_test\
	HPLC_test\
	ExperimentalSettings_test

#SYSTEM

EXECUTABLES_SYSTEM = \
	StopWatch_test\
	File_test

#KERNEL

EXECUTABLES_KERNEL =\
	RawDataPoint1D_test\
	RawDataPoint2D_test\
	DRawDataPoint_test\
	DPeak_test\
	Peak1D_test\
	Peak2D_test\
  PickedPeak1D_test\
	DPeakArray_test\
	DPeakConstReferenceArray_test\
	DSpectrum_test\
	ComparatorUtils_test\
	StandardTypes_test\
	Feature_test\
	FeatureMap_test\
	RangeUtils_test\
	MSSpectrum_test\
	MSExperiment_test\
	MSExperimentExtern_test\
	ConsensusFeature_test\
	ConsensusPeak_test\
	ConsensusMap_test\
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

ifeq ($(LIBSVM_SUPPORT),true)
	LIBSVM_EXECUTABLES_FORMAT=\
		LibSVMEncoder_test
endif

EXECUTABLES_FORMAT = \
	XMLValidator_test\
	TextFile_test\
	$(ANDI_EXECUTABLES_FORMAT)\
	${LIBSVM_EXECUTABLES_FORMAT}\
	DTAFile_test\
	DTA2DFile_test\
	MascotInfile_test\
	MascotOutfile_test\
	Base64_test\
	UniqueIdGenerator_test\
	PersistentObject_test\
	FASTAFile_test\
	FileHandler_test\
	MzDataExpSettHandler_test\
	MzXMLFile_test\
	MzDataFile_test\
	Param_test\
	FeatureXMLFile_test\
	FeaturePairsXMLFile_test\
	MascotXMLFile_test\
	GridFile_test\
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
	OMSSACSVFile_test\
	MSPFile_test\
	IdXMLFile_test\
	XTandemXMLFile_test\
	Schema2CV_test\
	ControlledVocabulary_test\
	XTandemXMLFile_test


#MATH

EXECUTABLES_MATH = \
	LinearInterpolation_test\
	AveragePosition_test\
	BasicStatistics_test\
	AsymmetricStatistics_test\
	ROCCurve_test\
	Histogram_test\
	MathFunctions_test\
	LinearRegression_test\
	BilinearInterpolation_test\

#FILTERING

EXECUTABLES_FILTERING = \
	LinearResampler_test\
	TopHatFilter_test\
	MorphFilter_test\
	SignalToNoiseEstimator_test\
	SignalToNoiseEstimatorMedian_test\
	SignalToNoiseEstimatorMeanIterative_test\
	SmoothFilter_test\
	GaussFilter_test\
	MaxReducer_test\
	SumReducer_test\
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
	SavitzkyGolayQRFilter_test\
	SavitzkyGolaySVDFilter_test\
	InternalCalibration_test\
	IDFilter_test\
	ExternalCalibration_test\
	ParentFilter_test\
	TradSeqQuality_test\
	PeakPosBins_test\
	BernNorm_test\
	PeakDiffBins_test\
	IntensityDistBins_test


#COMPARISON

EXECUTABLES_COMPARISON=\
	ClusterSpectrum_test\
	SpectrumCheapDPCorr_test\
	BinnedRep_test\
	HierarchicalClustering_test\
	BinnedRepCompareFunctor_test\
	BinnedRepMutualInformation_test\
	BinnedRepSpectrumContrastAngle_test\
	BinnedRepSharedPeakCount_test\
	BinnedRepSumAgreeingIntensities_test\
	PeakSpectrumCompareFunctor_test\
	SpectrumPrecursorComparator_test\
	ZhangSimilarityScore_test\
	SpectrumAlignment_test\
	SpectrumAlignmentScore_test

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
	TheoreticalSpectrumGenerator_test

#ANALYSIS

ifeq ($(LIBSVM_SUPPORT),true)
	EXECUTABLES_ANALYSIS_SVM = \
	SVMWrapper_test
endif



EXECUTABLES_ANALYSIS_MAPMATCHING = \
	IndexTuple_test\
	Group_test\
	GridCell_test\
	Grid_test\
	BaseMapping_test\
	LinearMapping_test\
	ElementPair_test\
	PairMatcher_test\
	BaseMapMatcher_test\
	MapMatcherRegression_test\
	MapDewarper_test\
	BasePairFinder_test\
	SimplePairFinder_test\
	BaseSuperimposer_test\
	PoseClusteringAffineSuperimposer_test\
	PoseClusteringShiftSuperimposer_test\
	BasePairwiseMapMatcher_test\
	PoseClusteringPairwiseMapMatcher_test\
	BaseAlignment_test\
	StarAlignment_test\
	DelaunayPairFinder_test

EXECUTABLES_ANALYSIS_DECHARGING = \
	FeatureDecharger_test

EXECUTABLES_ANALYSIS_ID =\
	IDSpectrumMapper_test\
	IDFeatureMapper_test\
	ConsensusID_test\
	HiddenMarkovModelLight_test\
	HiddenMarkovModel_test\
	ProtonDistributionModel_test\
	PILISSequenceDB_test\
	PILISModel_test\
	PILISModelGenerator_test\
	PILISScoring_test\
	PILISIdentification_test

# APPLICATIONS

EXECUTABLES_APPLICATIONS = \
	TOPPBase_test

# TRANSFORMATIONS

EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER = \
	ModelDescription_test\
	BaseModel_test\
	EmgModel_test\
	LmaGaussModel_test\
	LogNormalModel_test\
	InterpolationModel_test\
	GaussModel_test\
	BiGaussModel_test\
	IsotopeModel_test\
	ProductModel_test\
	FeatureFinder_test\
	FeaFiModule_test\
	FeatureFinderAlgorithm_test\
	SimpleSeeder_test\
	SimpleExtender_test\
	SimpleModelFitter_test

#	EuclideanDistance_test\
#	PickedPeakSeeder_test\
#	MarrWaveletSeeder_test\
#	IsotopeWaveletSeeder_test\
#	RankCorrelation_test\
#	BaseSweepSeeder_test\
#	AveragineMatcher_test

EXECUTABLES_TRANSFORMATIONS_PEAKPICKER = \
	ContinuousWaveletTransform_test\
	PeakShape_test\
	PeakPicker_test\
	ContinuousWaveletTransformNumIntegration_test\
	OptimizePick_test\
	PeakPickerCWT_test\
	OptimizePeakDeconvolution_test\
	TwoDOptimization_test

EXECUTABLES_PEPTIDEMASSINDEX = \
	FastaIterator_test \
	FastaIteratorIntern_test \
	EdwardsLippertIterator_test \
	EdwardsLippertIteratorTryptic_test \
	TrypticIterator_test \
	SuffixArrayTrypticCompressed_test \
	SuffixArraySeqan_test \
	SuffixArrayTrypticSeqan_test  \
	BigString_test \
  ModifierRep_test\
	SuffixArrayPeptideFinder_test

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
	$(EXECUTABLES_ANALYSIS_MAPMATCHING)\
	$(EXECUTABLES_APPLICATIONS)\
	$(EXECUTABLES_PEPTIDEMASSINDEX)


# default target
all:	$(EXECUTABLES)

# required for clean and depend
CPP_SOURCES=$(EXECUTABLES:=.C)
OBJECTS=$(EXECUTABLES:=.o)
VALGRINDS=$(EXECUTABLES:=.vgr)

include ../common.mak
-include .Dependencies

# remove the library specific flags
# and replace them by the application specific
# flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS)

# do not optimize tests (reduces compile time drastically)
CPP_MODE_FLAGS= $(CPP_MODE_FLAGS_NO_OPTIMIZATION) -O0


LIBS=$(GUI_LIBS)


.o:
	$(CXX) $(CXXFLAGS) $*.o $(LDFLAGS) -o $@ $(LIBS)

%%.C:
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(OPENMS_CXXFLAGS) $(OPENMS_INCLUDES) $*.C $(LDFLAGS) -o $@ $(LIBS)

.o.vgr:
	@$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(OPENMS_CXXFLAGS) $(OPENMS_INCLUDES) $*.o $(LDFLAGS) -o $@ $(LIBS)
	$(VALGRIND) $(VALGRIND_OPTS) $* >& $*.vgr || true

dummy:

# concept_test should be the first one, because it tests the testing ...
test: info concept_test system_test math_test datastructures_test metadata_test kernel_test format_test filtering_test visual_test comparison_test chemistry_test transformations_test analysis_id_test analysis_decharging_test analysis_svm_test analysis_mapmatching_test applications_test
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

transformations_test:	$(EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER) $(EXECUTABLES_TRANSFORMATIONS_PEAKPICKER) dummy
	@./runtests -section "TRANSFORMATIONS" $(TEST_OPTIONS) $(EXECUTABLES_TRANSFORMATIONS_FEATUREFINDER) $(EXECUTABLES_TRANSFORMATIONS_PEAKPICKER)

analysis_id_test:	$(EXECUTABLES_ANALYSIS_ID) dummy
	@./runtests -section "ANALYSIS_ID" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_ID)

analysis_decharging_test:	$(EXECUTABLES_ANALYSIS_DECHARGING) dummy
	@./runtests -section "ANALYSIS_DECHARGING" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_DECHARGING)

analysis_svm_test:	$(EXECUTABLES_ANALYSIS_SVM) dummy
	@./runtests -section "ANALYSIS_SVM" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_SVM)

analysis_mapmatching_test:	$(EXECUTABLES_ANALYSIS_MAPMATCHING) dummy
	@./runtests -section "ANALYSIS_MAPMATCHING" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_MAPMATCHING)

applications_test:	$(EXECUTABLES_APPLICATIONS) dummy
	@./runtests -section "APPLICATIONS" $(TEST_OPTIONS) $(EXECUTABLES_APPLICATIONS)

peptidemassindex_test:	$(EXECUTABLES_PEPTIDEMASSINDEX) dummy
	@./runtests -section "PEPTIDEMASSINDEX" $(TEST_OPTIONS) $(EXECUTABLES_PEPTIDEMASSINDEX)

clean:	default_clean
	@$(RM) -rf -- *.o *~ *.log ii_files Templates.DB	SunWS_cache *.TMP *.tmp core a.out region* cut* feature* model* *.vgr *.output $(EXECUTABLES)
	@cd TOPP && make clean

# Create data files for MzData, MzXML, ANDIFile with variable number of scans
# Call: format_data [#scans] [filename] (i.e. format_data 100 tmp)
format_data:	format_data.o
	$(CXX)	$(LDFLAGS) -o format_data format_data.o $(LIBS)

#creates a .vgr-file with valgrind check output for each test
valgrind: $(VALGRINDS)
