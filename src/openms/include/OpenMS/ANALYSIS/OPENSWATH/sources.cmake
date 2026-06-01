### the directory name
set(directory include/OpenMS/ANALYSIS/OPENSWATH)

### list all header files of the directory here
set(sources_list_h
  CalibrationWorkflow.h
  ChromatogramExtractor.h
  ChromatogramExtractorAlgorithm.h
  ConfidenceScoring.h
  DIAHelper.h
  DIAPrescoring.h
  DIAScoring.h
  IonMobilityScoring.h
  LevelContextInference.h
  MasstraceCorrelator.h
  MRMAssay.h
  MRMDecoy.h
  MRMFeatureFilter.h
  MRMFeatureFinderScoring.h
  MRMFeaturePicker.h
  MRMFeatureQC.h
  MRMBatchFeatureSelector.h
  MRMFeatureSelector.h
  MRMIonSeries.h
  MRMRTNormalizer.h
  MRMScoring.h
  MRMTransitionGroupPicker.h
  OpenSwathGeneInference.h
  OpenSwathHelper.h
  OpenSwathExportConfig.h
  OpenSwathExportData.h
  OpenSwathMatrixExporter.h
  OpenSwathParquetExporter.h
  OpenSwathInferenceConfig.h
  OpenSwathInferenceData.h
  OpenSwathScores.h
  OpenSwathScoring.h
  OpenSwathOSWParquetReader.h
  OpenSwathOSWParquetWriter.h
  OpenSwathOSWWriter.h
  OpenSwathResultsExporter.h
  OpenSwathWorkflow.h
  OpenSwathWorkflowScheduler.h
  PeakIntegrator.h
  PeakPickerChromatogram.h
  PeakPickerMobilogram.h
  OpenSwathPeptideInference.h
  OpenSwathPeptidoformInference.h
  OpenSwathProteinInference.h
  SwathMapMassCorrection.h
  SwathWindowLoader.h
  SwathQC.h
  SpectrumAddition.h
  TargetedSpectraExtractor.h
  TransitionTSVFile.h
  TransitionParquetFile.h
  TransitionPQPFile.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\OPENSWATH" FILES ${sources_h})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
