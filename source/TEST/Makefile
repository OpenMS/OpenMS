
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

# ClassTest_test should be the first one, because it tests the testing ...
EXECUTABLES_CONCEPT = \
	ClassTest_test\
	Exception_Base_test\
	PreciseTime_test\
	TimeStamp_test\
	VersionInfo_test\
	FactoryProduct_test\
	Factory_test

EXECUTABLES_DATASTRUCTURES = \
	${CGAL_EXECUTABLES_DATASTRUCTURES}\
	String_test\
	Date_test\
	HashMap_test\
	IndexSet_test\
	RunningAveragePosition_test\
	ScanIndex_test\
	DIntervalBase_test\
	DRange_test\
	DBoundingBox_test\
	DataValue_test\
	Matrix_test\
	SparseVector_test\
	BinnedSparseVector_test\
	DateTime_test\
	QuadTree_test\
	RangeManager_test

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
	Identification_test\
	ProteinIdentification_test\
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

EXECUTABLES_SYSTEM = \
	StopWatch_test

EXECUTABLES_KERNEL =\
	DPosition_test\
	DRawDataPoint_test\
	DPeak_test\
	DPickedPeak_test\
	DPeakArray_test\
	DPeakList_test\
	DPeakListIterator_test\
	DPeakArrayNonPolymorphic_test\
	DPeakArrayIterator_test\
	DSpectrum_test\
	ComparatorUtils_test\
	StandardTypes_test\
	DFeature_test\
	DFeatureMap_test\
	RangeUtils_test\
	MSSpectrum_test\
	MSExperiment_test\
	MSExperimentExtern_test

#VISUAL

ifeq ($(ENABLE_GUI),true)
	EXECUTABLES_VISUAL =\
		MultiGradient_test\
		MappingInfo_test
endif

#FORMAT

ifeq ($(ANDIMS_SUPPORT),true)
	ifeq ($(NETCDF_SUPPORT),true)
		ANDI_EXECUTABLES_FORMAT=\
			ANDIFile_test
	endif
endif

#conditional tests
ifeq ($(ENABLE_DB),true)
	DB_EXECUTABLES_FORMAT=\
	DBConnection_test\
	DBAdapter_test
endif

EXECUTABLES_FORMAT = \
	RTTI_test\
	TextFile_test\
	$(ANDI_EXECUTABLES_FORMAT)\
	${DB_EXECUTABLES_FORMAT}\
	DTAFile_test\
	DTA2DFile_test\
	MascotInfile_test\
	MascotOutfile_test\
	Param_test\
	Base64_test\
	MzXMLFile_test\
	MzDataFile_test\
	UniqueIdGenerator_test\
	PersistentObject_test\
	PersistenceManager_test\
	DFeatureMapFile_test\
	DFeaturePairsFile_test\
	FASTAFile_test\
	MzDataExpSettHandler_test\
	AnalysisXMLFile_test\
	MascotXMLFile_test\
	DGridFile_test\
	FileHandler_test

#MATH

ifeq ($(GSL_SUPPORT),true)
	GSL_EXECUTABLES_MATH = \
		LinearRegression_test
endif

EXECUTABLES_MATH = \
	LinearInterpolation_test\
	BasicStatistics_test\
	ROCCurve_test\
	Histogram_test\
	MathFunctions_test\
	${GSL_EXECUTABLES_MATH}

#FILTERING

ifeq ($(GSL_SUPPORT),true)
	GSL_EXECUTABLES_FILTERING = \
		DSavitzkyGolaySVDFilter_test\
		DSavitzkyGolayQRFilter_test
endif

EXECUTABLES_FILTERING = \
	IDFilter_test\
	LinearResampler_test\
	TopHatFilter_test\
	MorphFilter_test\
	DSignalToNoiseEstimatorWindowing_test\
	PreprocessingFunctor_test\
	NLargest_test\
	Scaler_test\
	BernNorm_test\
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
	ParentFilter_test\
	TradSeqQuality_test\
	PeakPosBins_test\
	PeakDiffBins_test\
	IntensityDistBins_test\
	DSignalToNoiseEstimatorMedian_test\
	KellerQuality_test
	
#	DSmoothFilter_test\
	#${GSL_EXECUTABLES_FILTERING}\
	#DGaussFilter_test\

#COMPARISON

EXECUTABLES_COMPARISON = \
	BinnedRep_test\
	ClNode_test\
	Cluster_test\
	ClusterSpectrum_test\
	ClusterFunctor_test\
  ClusterExperiment_test\
	AnalysisFunctor_test\
	CompareFunctor_test\
	BinnedRepMutualInformation_test\
	BinnedRepSpectrumContrastAngle_test\
	SpectrumCheapDPCorr_test\
	BinnedRepSharedPeakCount_test\
	BinnedRepSumAgreeingIntensities_test\
	SequestCompareFunctor_test\
	SpectrumPrecursorComparator_test

EXECUTABLES_CHEMISTRY = \
	IsotopeDistribution_test\
	Element_test\
	ElementDB_test\
	EmpiricalFormula_test\
	Residue_test\
	ResidueDB_test\
	PeptideSequence_test


ifeq ($(GSL_SUPPORT),true)
	GSL_EXECUTABLES_ANALYSIS_MAPMATCHING = \
		DMapMatcherRegression_test\
		DMapDewarper_test
endif


EXECUTABLES_ANALYSIS_MAPMATCHING = \
	DFeaturePair_test\
	DFeaturePairVector_test\
	DGridCell_test\
	DGrid_test\
	DLinearMapping_test\
	$(GSL_EXECUTABLES_ANALYSIS_MAPMATCHING)

# TRANSFORMATIONS

ifeq ($(GSL_SUPPORT),true)
	GSL_EXECUTABLES_TRANSFORMATIONS = \
		ContinuousWaveletTransformNumIntegration_test\
		OptimizePick_test

endif

EXECUTABLES_TRANSFORMATIONS = \
	${GSL_EXECUTABLES_TRANSFORMATIONS}\
	ContinuousWaveletTransform_test\
	PeakShape_test\
	BaseModel_test\
	GaussModel_test\
	BiGaussModel_test\
	IsotopeModel_test\
	ProductModel_test\
	SimpleSeeder_test\
	SimpleExtender_test\
	SimpleModelFitter_test\
	FeaFiTraits_test\
	PeakPicker_test\
	PeakPickerCWT_test\
	ExtractSignalRegions_test


EXECUTABLES_ANALYSIS_ID =\
	MSExperimentAnnotator_test

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
	$(EXECUTABLES_TRANSFORMATIONS)\
	$(EXECUTABLES_ANALYSIS_ID)\
	$(EXECUTABLES_ANALYSIS_MAPMATCHING)


# default target
all:	$(EXECUTABLES)

# test for completeness of tests
check_test:	../config/check_test
	@L=`wc -l check_test.lst|awk '{print $$1}'` ; I=1 ; while test $${I} -le $${L} ; do ../config/check_test ../../include/OpenMS/`cat check_test.lst|sed -n $${I}p` ; I=`expr $${I} + 1`; done

# required for clean and depend
CPP_SOURCES=$(EXECUTABLES:=.C)
OBJECTS=$(EXECUTABLES:=.o)

include ../common.mak
-include .Dependencies

# remove the library specific flags
# and replace them by the application specific
# flags
LIB_CXXFLAGS=$(NONLIB_CXXFLAGS)

# do not optimize tests (reduces compile time drastically)
CPP_MODE_FLAGS=$(CPP_MODE_FLAGS_NO_OPTIMIZATION) -O0

#add libraries for GUI to LIBS
ifeq ($(ENABLE_GUI),true)
	LIBS=$(GUI_LIBS)
	OPENMS_INCLUDES=$(GUI_INCLUDES)
endif

.o:
	$(CXX) $(CXXFLAGS) $*.o $(LDFLAGS) -o $@ $(LIBS)
# 
%%.C:
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) $(CPP_MODE_FLAGS) $(OPENMS_CXXFLAGS) $(OPENMS_INCLUDES) $*.C $(LDFLAGS) -o $@ $(LIBS)

dummy:

# concept_test should be the first one, because it tests the testing ...
test: info concept_test system_test math_test datastructures_test metadata_test kernel_test format_test filtering_test visual_test comparison_test chemistry_test transformations_test analysis_id_test analysis_mapmatching_test
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

transformations_test:	$(EXECUTABLES_TRANSFORMATIONS) dummy
	@./runtests -section "TRANSFORMATIONS" $(TEST_OPTIONS) $(EXECUTABLES_TRANSFORMATIONS)

analysis_id_test:	$(EXECUTABLES_ANALYSIS_ID) dummy
	@./runtests -section "ANALYSIS_ID" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_ID)
	
analysis_mapmatching_test:	$(EXECUTABLES_ANALYSIS_MAPMATCHING) dummy
	@./runtests -section "ANALYSIS_MAPMATCHING" $(TEST_OPTIONS) $(EXECUTABLES_ANALYSIS_MAPMATCHING)

clean:	default_clean
	@$(RM) -rf *.o *~ *.log ii_files Templates.DB	SunWS_cache _*.TMP core a.out $(EXECUTABLES)

# Create data files for MzData, MzXML, ANDIFile with variable number of scans
# Call: format_data [#scans] [filename] (i.e. format_data 100 tmp)
format_data:	format_data.o 
	$(CXX)	$(LDFLAGS) -o format_data format_data.o $(LIBS)
