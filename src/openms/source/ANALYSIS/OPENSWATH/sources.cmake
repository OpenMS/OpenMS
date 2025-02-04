### the directory name
set(directory source/ANALYSIS/OPENSWATH)

### list all header files of the directory here
set(sources_list
  ChromatogramExtractor.cpp
  ChromatogramExtractorAlgorithm.cpp
  ConfidenceScoring.cpp
  DIAHelper.cpp
  DIAPrescoring.cpp
  DIAScoring.cpp
  IonMobilityScoring.cpp
  MasstraceCorrelator.cpp
  MRMAssay.cpp
  MRMDecoy.cpp
  MRMFeatureFilter.cpp
  MRMFeatureFinderScoring.cpp
  MRMFeaturePicker.cpp
  MRMFeatureQC.cpp
  MRMBatchFeatureSelector.cpp
  MRMFeatureSelector.cpp
  MRMIonSeries.cpp
  MRMRTNormalizer.cpp
  MRMScoring.cpp
  MRMTransitionGroupPicker.cpp
  OpenSwathHelper.cpp
  OpenSwathScores.cpp
  OpenSwathScoring.cpp
  OpenSwathOSWWriter.cpp
  OpenSwathWorkflow.cpp
  PeakIntegrator.cpp
  PeakPickerChromatogram.cpp
  SwathMapMassCorrection.cpp
  SwathWindowLoader.cpp
  SwathQC.cpp
  SpectrumAddition.cpp
  TargetedSpectraExtractor.cpp
  TransitionTSVFile.cpp
  TransitionPQPFile.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\OPENSWATH" FILES ${sources})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)
