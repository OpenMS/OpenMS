// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>
#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ML/SVM/SimpleSVM.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/COMPARISON/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/PeakAlignment.h>
#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>
#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResampler.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/PROCESSING/SCALING/SqrtScaler.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/MATH/MISC/EmgGradientDescent.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/QC/DBSuitability.h>
#include <OpenMS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/FEATUREFINDER/MaxLikeliFitter1D.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>

// those are only added if GUI is enabled
#ifdef WITH_GUI
#include <QApplication>

#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>
#include <OpenMS/VISUAL/Plot3DCanvas.h>
#include <OpenMS/VISUAL/SpectraIDViewTab.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#endif


#include <fstream>

using namespace std;
using namespace OpenMS;

// this weird piece of code is required to avoid the following linker errors in VS2019
/*
Error	LNK2001	unresolved external symbol "public: virtual void * __cdecl OpenMS::MSExperiment::`scalar deleting destructor'(unsigned int)" (??_GMSExperiment@OpenMS@@UEAAPEAXI@Z)	DefaultParamHandlerDocumenter	C:\dev\openms_test_build19\doc\DefaultParamHandlerDocumenter.obj	1	
Error	LNK2019	unresolved external symbol "public: virtual void * __cdecl OpenMS::MSExperiment::`vector deleting destructor'(unsigned int)" (??_EMSExperiment@OpenMS@@UEAAPEAXI@Z) referenced in function "[thunk]:public: virtual void * __cdecl OpenMS::MSExperiment::`vector deleting destructor'`adjustor{72}' (unsigned int)" (??_EMSExperiment@OpenMS@@WEI@EAAPEAXI@Z)	DefaultParamHandlerDocumenter	C:\dev\openms_test_build19\doc\DefaultParamHandlerDocumenter.obj	1
see https://stackoverflow.com/a/74235019/1913074
Alternatively, define ~MSExperiment(){}; instead of using ' = default;'
*/
void foo()
{
  auto p = new MSExperiment();
  delete p;
}

//**********************************************************************************
//Helper method - use this method to generate the actual parameter documentation
//**********************************************************************************
void writeParameters(const String& class_name, const Param& param, bool table_only = false)
{
  const String filename = String("output/OpenMS_") + class_name + ".parameters";
  ofstream f(filename.c_str());

  if (!f)
  {
    std::cerr << "Cannot open file '" << filename << "'. Check for invalid characters in filename and permissions.\n";
    exit(1);
  }
    
  if (!table_only)
  {
    f << "<B>Parameters of this class are:</B><BR><BR>\n";
  }
  f << R"(<table class="doxtable" border="1" width="100%" cellpadding="4">)" << endl;
  f << "<tr><th>Name</th><th>Type</th><th>Default</th><th>Restrictions</th><th>Description</th></tr>" << endl;
  String type, description, restrictions;
  for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
  {
    restrictions = "";
    if (it->value.valueType() == ParamValue::INT_VALUE || it->value.valueType() == ParamValue::INT_LIST)
    {
      type = "int";
      if (it->value.valueType() == ParamValue::INT_LIST)
      {
        type += " list";
      }

      //restrictions
      bool first = true;
      if (it->min_int != -(numeric_limits<Int>::max)())
      {
        restrictions += String("min: ") + it->min_int;
        first = false;
      }
      if (it->max_int != (numeric_limits<Int>::max)())
      {
        if (!first)
        {
          restrictions += ' ';
        }          
        restrictions += String("max: ") + it->max_int;
      }
    }
    else if (it->value.valueType() == ParamValue::DOUBLE_VALUE || it->value.valueType() == ParamValue::DOUBLE_LIST)
    {
      type = "float";
      if (it->value.valueType() == ParamValue::DOUBLE_LIST)
        type += " list";

      //restrictions
      bool first = true;
      if (it->min_float != -(numeric_limits<double>::max)())
      {
        restrictions += String("min: ") + it->min_float;
        first = false;
      }
      if (it->max_float != (numeric_limits<double>::max)())
      {
        if (!first)
          restrictions += ' ';
        restrictions += String("max: ") + it->max_float;
      }
    }
    else if (it->value.valueType() == ParamValue::STRING_VALUE || it->value.valueType() == ParamValue::STRING_LIST)
    {
      type = "string";
      if (it->value.valueType() == ParamValue::STRING_LIST)
        type += " list";

      //restrictions
      if (!it->valid_strings.empty())
      {
        String valid_strings;
        valid_strings.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ", ");
        restrictions += valid_strings;
      }
    }
    if (restrictions == "")
    {
      restrictions = "&nbsp;";
    }
    //replace #, @ and newline in description
    description = param.getDescription(it.getName());
    description.substitute("@", "XXnot_containedXX");
    description.substitute("XXnot_containedXX", "@@");
    description.substitute("#", "XXnot_containedXX");
    description.substitute("XXnot_containedXX", "@#");
    description.substitute("\n", "<BR>");

    //create tooltips for sections if they are documented
    String name = it.getName();
    vector<String> parts;
    name.split(':', parts);
    String prefix = "";
    for (Size i = 0; i + 1 < parts.size(); ++i)
    {
      if (i == 0)
      {
        prefix = parts[i];
      }
      else
      {
        prefix = prefix + ":" + parts[i];
      }
      String docu = param.getSectionDescription(prefix);
      if (docu != "")
      {
        parts[i] = String("<span title=\"") + docu + "\">" + parts[i] + "</span>";
      }
    }
    if (parts.size() != 0)
    {
      name.concatenate(parts.begin(), parts.end(), ":");
    }

    //replace # and @ in values
    String value = it->value.toString(true);
    value.substitute("@", "XXnot_containedXX");
    value.substitute("XXnot_containedXX", "@@");
    value.substitute("#", "XXnot_containedXX");
    value.substitute("XXnot_containedXX", "@#");

    //make the advanced parameters cursive, the normal ones bold
    String style = "b";
    if (it->tags.count("advanced") == 1)
      style = "i";

    //final output
    f << "<tr>\n"
      << "  <td style=\"vertical-align:top\"><" << style << ">" << name << "</" << style << "></td>\n"
      << "  <td style=\"vertical-align:top\">" << type << "</td><td style=\"vertical-align:top\">" << value <<  "</td>\n"
      << "  <td style=\"vertical-align:top\">" << restrictions << "</td><td style=\"vertical-align:top\">" << description <<  "</td>\n"
      << "</tr>\n";
  }
  f << "</table>" << "\n";
  if (!table_only)
  {
    f << "<br>" << "\n"
      << "<b>Note:</b>" << "\n"
      << "<UL style=\"margin-top:0px;\">" << "\n"
      << "  <LI> If a section name is documented, the documentation is displayed as tooltip." << "\n"
      << "  <LI> Advanced parameter names are italic." << "\n"
      << "</UL>" << "\n";
  }
  f.close();
}

//**********************************************************************************
//Helper macros that can be used for easy classes
//**********************************************************************************

// For classes that have a default-constructor, simply use this macro with the
// class name
#define DOCME(class_name) \
  writeParameters("" # class_name, class_name().getDefaults());

// For class templates and classes without default constructor use this macro
// with (1.) the class name and (2.) a class instance.
#define DOCME2(class_template_name, instantiation) \
  writeParameters("" # class_template_name, (instantiation).getDefaults());

//**********************************************************************************
//Main method - add your class here
//**********************************************************************************
int main(int argc, char** argv)
{
  //////////////////////////////////
  // Simple cases
  //////////////////////////////////

  DOCME(AScore);
  DOCME(BasicProteinInferenceAlgorithm);
  DOCME(BayesianProteinInferenceAlgorithm);
  DOCME(TransitionPQPFile);
  DOCME(BiGaussFitter1D);
  DOCME(BiGaussModel);
  DOCME(BinnedSharedPeakCount);
  DOCME(BinnedSpectralContrastAngle);
  DOCME(BinnedSumAgreeingIntensities);

  DOCME(ConsensusIDAlgorithmAverage);
  DOCME(ConsensusIDAlgorithmBest);
  DOCME(ConsensusIDAlgorithmPEPIons);
  DOCME(ConsensusIDAlgorithmPEPMatrix);
  DOCME(ConsensusIDAlgorithmRanks);
  DOCME(ConsensusIDAlgorithmWorst);
  DOCME(DBSuitability);
  DOCME(DiaPrescore);
  DOCME(DIAScoring);
  DOCME(ElutionModelFitter);
  DOCME(EmgFitter1D);
  DOCME(EmgGradientDescent)
  DOCME(EmgModel);
  DOCME(ExtendedIsotopeFitter1D);
  DOCME(ExtendedIsotopeModel);
  DOCME(FalseDiscoveryRate);
  DOCME(FeatureDeconvolution);
  DOCME(FeatureDistance);
  DOCME(FeatureFinderAlgorithmMetaboIdent);
  DOCME(ElutionPeakDetection);
  DOCME(FeatureFindingMetabo);
  DOCME(FeatureGroupingAlgorithmLabeled);
  DOCME(FeatureGroupingAlgorithmQT);
  DOCME(FeatureGroupingAlgorithmKD);
  DOCME(FeatureGroupingAlgorithmUnlabeled);
  DOCME(MapAlignmentAlgorithmIdentification);
  DOCME(MapAlignmentAlgorithmTreeGuided);
  DOCME(MassTraceDetection);
  DOCME(FIAMSDataProcessor);
  DOCME(GaussFilter);
  DOCME(GaussFitter1D);
  DOCME(GaussModel);
  DOCME(IDMapper);
  DOCME(IDRipper);
  DOCME(InterpolationModel);
  DOCME(IsotopeFitter1D);
  DOCME(IsotopeModel);
  DOCME(TMTSixPlexQuantitationMethod);
  DOCME(TMTTenPlexQuantitationMethod);
  DOCME(TMTSixteenPlexQuantitationMethod);
  DOCME(TMTEighteenPlexQuantitationMethod);
  DOCME(ItraqEightPlexQuantitationMethod);
  DOCME(ItraqFourPlexQuantitationMethod);
  DOCME(LabeledPairFinder);
  DOCME(LinearResampler);
  DOCME(MSPFile);
  DOCME(MapAlignmentAlgorithmPoseClustering);
  DOCME(SpectrumAnnotator);
  DOCME(TheoreticalSpectrumGeneratorXLMS);
  DOCME(MRMDecoy);
  DOCME(MetaboliteFeatureDeconvolution);
  DOCME(MRMFeatureFilter);
  DOCME(MRMFeatureFinderScoring);
  DOCME(MRMTransitionGroupPicker);
  DOCME(MultiplexDeltaMassesGenerator);
  DOCME(NucleicAcidSpectrumGenerator);
  DOCME(NLargest);
  DOCME(Normalizer);
  DOCME(PeakAlignment);
  DOCME(PeakIntegrator);
  DOCME(PeakPickerHiRes);
  DOCME(PeakPickerIterative);
  DOCME(PeakPickerChromatogram);
  DOCME(PeptideIndexing);
  DOCME(PoseClusteringAffineSuperimposer);
  DOCME(PoseClusteringShiftSuperimposer);
  DOCME(QTClusterFinder);
  DOCME(SavitzkyGolayFilter);
  DOCME(LowessSmoothing);
  DOCME(SimpleSVM);
  DOCME(StablePairFinder);
  DOCME(SpectrumAlignment);
  DOCME(SpectrumAlignmentScore);
  DOCME(SpectrumCheapDPCorr);
  DOCME(SpectrumPrecursorComparator);
  DOCME(SqrtScaler);
  DOCME(SteinScottImproveScore);
  DOCME(SpectraMerger);
  DOCME(TheoreticalSpectrumGenerator);
  DOCME(ThresholdMower);
  DOCME(TransitionTSVFile);
  DOCME(IDDecoyProbability);
  DOCME(WindowMower);
  DOCME(ZhangSimilarityScore);
  DOCME(MorphologicalFilter);
  DOCME(MassDecompositionAlgorithm);
  DOCME(MascotRemoteQuery);
  DOCME(MascotGenericFile);
  DOCME(Fitter1D);  
  DOCME(PeptideAndProteinQuant);
  DOCME(SimpleTSGXLMS);
  // workarounds for documenting model parameters in MapAligners:
  writeParameters("MapAlignerIdentificationModel", MapAlignerBase::getModelDefaults("interpolated"), true);
  writeParameters("MapAlignerPoseClusteringModel", MapAlignerBase::getModelDefaults("linear"), true);
  writeParameters("MapRTTransformerModel", MapAlignerBase::getModelDefaults("none"), true);

  //////////////////////////////////
  // More complicated cases
  //////////////////////////////////

  // ConsensusIDAlgorithm...: abstract base classes, get param. from subclass:
  DOCME2(ConsensusIDAlgorithm, (ConsensusIDAlgorithmBest()));
  DOCME2(ConsensusIDAlgorithmIdentity, (ConsensusIDAlgorithmBest()));
  DOCME2(ConsensusIDAlgorithmSimilarity, (ConsensusIDAlgorithmBest()));
  DOCME2(FeatureFinderAlgorithmPicked, (FeatureFinderAlgorithmPicked()));
  DOCME2(SignalToNoiseEstimatorMeanIterative, SignalToNoiseEstimatorMeanIterative<>());
  DOCME2(SignalToNoiseEstimatorMedian, SignalToNoiseEstimatorMedian<>());
  DOCME2(SignalToNoiseEstimator, SignalToNoiseEstimatorMedian<>()); //SignalToNoiseEstimator is a base class, get parameters from subclass SignalToNoiseEstimatorMedian
  DOCME2(GaussTraceFitter, (GaussTraceFitter()));
  DOCME2(EGHTraceFitter, (EGHTraceFitter()));
  DOCME2(TraceFitter, (GaussTraceFitter())); //TraceFitter is an abstract base class, get parameters from subclass GaussTraceFitter
  DOCME2(BinnedSpectrumCompareFunctor, (BinnedSharedPeakCount())); //BaseModel is a base class, get parameters from subclass BinnedSharedPeakCount
  ItraqFourPlexQuantitationMethod itraq4;
  DOCME2(IsobaricChannelExtractor, (IsobaricChannelExtractor(&itraq4)))
  DOCME2(IsobaricQuantifier, (IsobaricQuantifier(&itraq4)))
  DOCME2(PosteriorErrorProbabilityModel, Math::PosteriorErrorProbabilityModel());
  

  // handle GUI documentation separately
#ifdef WITH_GUI
  // some classes require a QApplication
  QApplication app(argc, argv);

  DOCME(TOPPASBase);

  DOCME2(TOPPViewBase, TOPPViewBase(TOPPViewBase::TOOL_SCAN::SKIP_SCAN));
  DOCME2(PlotCanvas, Plot1DCanvas(Param()));
  DOCME2(Plot1DCanvas, Plot1DCanvas(Param()));
  DOCME2(Plot2DCanvas, Plot2DCanvas(Param()));
  DOCME2(Plot3DCanvas, Plot3DCanvas(Param()));
  DOCME2(SpectraIDViewTab, SpectraIDViewTab(Param()));

#endif

  return 0;
}
