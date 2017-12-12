// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>
//#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorTrainer.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/SIMULATION/DigestSimulation.h>
#include <OpenMS/SIMULATION/IonizationSimulation.h>
#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/SIMULATION/MSSim.h>
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/SIMULATION/EGHFitter1D.h>
#include <OpenMS/SIMULATION/EGHModel.h>
#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>
#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>
#include <OpenMS/SIMULATION/LABELING/ICPLLabeler.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

// those are only added if GUI is enabled
#ifdef WITH_GUI
#include <QtGui/QApplication>

#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#endif

// include this file after the GUI stuff, or there will be a conflict between
// "LayerData.h" (via "Spectrum1DCanvas.h") and "SeqanIncludeWrapper.h"!
// (see https://github.com/OpenMS/OpenMS/issues/1327)
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>


#include <fstream>

using namespace std;
using namespace OpenMS;

//**********************************************************************************
//Helper method - use this method to generate the actual parameter documentation
//**********************************************************************************
void writeParameters(const String& class_name, const Param& param, bool table_only = false)
{
  ofstream f((String("output/OpenMS_") + class_name + ".parameters").c_str());

  if (!table_only)
    f << "<B>Parameters of this class are:</B><BR><BR>\n";
  f << "<table border=\"1\" style=\"border-style:solid; border-collapse:collapse; border-color:#c0c0c0;\" width=\"100%\" cellpadding=\"4\">" << endl;
  f << "<tr style=\"border-bottom:1px solid black; background:#fffff0\"><th>Name</th><th>Type</th><th>Default</th><th>Restrictions</th><th>Description</th></tr>" << endl;
  String type, description, restrictions;
  for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
  {
    restrictions = "";
    if (it->value.valueType() == DataValue::INT_VALUE || it->value.valueType() == DataValue::INT_LIST)
    {
      type = "int";
      if (it->value.valueType() == DataValue::INT_LIST)
        type += " list";

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
          restrictions += ' ';
        restrictions += String("max: ") + it->max_int;
      }
    }
    else if (it->value.valueType() == DataValue::DOUBLE_VALUE || it->value.valueType() == DataValue::DOUBLE_LIST)
    {
      type = "float";
      if (it->value.valueType() == DataValue::DOUBLE_LIST)
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
    else if (it->value.valueType() == DataValue::STRING_VALUE || it->value.valueType() == DataValue::STRING_LIST)
    {
      type = "string";
      if (it->value.valueType() == DataValue::STRING_LIST)
        type += " list";

      //restrictions
      if (it->valid_strings.size() != 0)
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
    String value = it->value;
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
  DOCME(BernNorm);
  DOCME(BiGaussFitter1D);
  DOCME(BiGaussModel);
  DOCME(BinnedSharedPeakCount);
  DOCME(BinnedSpectralContrastAngle);
  DOCME(BinnedSumAgreeingIntensities);
  DOCME(ComplementFilter);
  DOCME(ComplementMarker);
  DOCME(ConsensusIDAlgorithmAverage);
  DOCME(ConsensusIDAlgorithmBest);
  DOCME(ConsensusIDAlgorithmPEPIons);
  DOCME(ConsensusIDAlgorithmPEPMatrix);
  DOCME(ConsensusIDAlgorithmRanks);
  DOCME(ConsensusIDAlgorithmWorst);
  DOCME(DetectabilitySimulation);
  DOCME(DIAScoring);
  DOCME(DigestSimulation);
  DOCME(ElutionModelFitter);
  DOCME(EmgFitter1D);
  DOCME(EmgModel);
  DOCME(ExtendedIsotopeFitter1D);
  DOCME(ExtendedIsotopeModel);
  DOCME(FalseDiscoveryRate);
  DOCME(FeatureDeconvolution);
  DOCME(FeatureDistance);
  DOCME(FeatureGroupingAlgorithmLabeled);
  DOCME(FeatureGroupingAlgorithmQT);
  DOCME(FeatureGroupingAlgorithmUnlabeled);
  DOCME(GaussFilter);
  DOCME(GaussFitter1D);
  DOCME(GaussModel);
  DOCME(GoodDiffFilter);
  DOCME(IDMapper);
  DOCME(InterpolationModel);
  DOCME(IsotopeDiffFilter);
  DOCME(IsotopeFitter1D);
  DOCME(IsotopeMarker);
  DOCME(IsotopeModel);
  DOCME(TMTSixPlexQuantitationMethod);
  DOCME(ItraqEightPlexQuantitationMethod);
  DOCME(ItraqFourPlexQuantitationMethod);
  DOCME(LabeledPairFinder);
  DOCME(LinearResampler);
  DOCME(MSPFile);
  DOCME(MSSim);
  DOCME(MapAlignmentAlgorithmPoseClustering);
  DOCME(MapAlignmentAlgorithmSpectrumAlignment);
  DOCME(MRMFeatureFinderScoring);
  DOCME(MRMTransitionGroupPicker);
  DOCME(NLargest);
  DOCME(NeutralLossDiffFilter);
  DOCME(NeutralLossMarker);
  DOCME(Normalizer);
  DOCME(OptimizePeakDeconvolution);
  DOCME(ParentPeakMower);
  DOCME(PeakAlignment);
  DOCME(PeakPickerCWT);
  DOCME(PeakPickerHiRes);
  DOCME(PeakPickerMRM);
  DOCME(PoseClusteringAffineSuperimposer);
  DOCME(PoseClusteringShiftSuperimposer);
  DOCME(QTClusterFinder);
  DOCME(SavitzkyGolayFilter);
  DOCME(LowessSmoothing);
  DOCME(SimplePairFinder);
  DOCME(SimpleSVM);
  DOCME(StablePairFinder);
  DOCME(SpectrumAlignment);
  DOCME(SpectrumAlignmentScore);
  DOCME(SpectrumCheapDPCorr);
  DOCME(SpectrumPrecursorComparator);
  DOCME(SteinScottImproveScore);
  DOCME(SpectraMerger);
  DOCME(SvmTheoreticalSpectrumGenerator);
  //DOCME(SvmTheoreticalSpectrumGeneratorSet);
  DOCME(SvmTheoreticalSpectrumGeneratorTrainer);
  DOCME(TICFilter);
  DOCME(TheoreticalSpectrumGenerator);
  DOCME(ThresholdMower);
  DOCME(TwoDOptimization);
  DOCME(WindowMower);
  DOCME(ZhangSimilarityScore);
  DOCME(PrecursorIonSelection);
  DOCME(PrecursorIonSelectionPreprocessing);
  DOCME(MorphologicalFilter);
  DOCME(CompNovoIonScoring);
  DOCME(CompNovoIonScoringCID);
  DOCME(CompNovoIdentification);
  DOCME(CompNovoIdentificationCID);
  DOCME(MassDecompositionAlgorithm);
  DOCME(MRMFragmentSelection);
  DOCME(ProtonDistributionModel);
  DOCME(MascotRemoteQuery);
  DOCME(MascotGenericFile);
  DOCME(OfflinePrecursorIonSelection);
  DOCME(Fitter1D);
  DOCME(EGHModel);
  DOCME(EGHFitter1D);
  DOCME(O18Labeler);
  DOCME(ITRAQLabeler);
  DOCME(SILACLabeler);
  DOCME(ICPLLabeler);
  DOCME(PeptideAndProteinQuant);
  DOCME(Math::PosteriorErrorProbabilityModel);
  // workarounds for documenting model parameters in MapAligners:
  writeParameters("MapAlignerIdentificationModel", TOPPMapAlignerBase::getModelDefaults("interpolated"), true);
  writeParameters("MapAlignerPoseClusteringModel", TOPPMapAlignerBase::getModelDefaults("linear"), true);
  writeParameters("MapAlignerSpectrumModel", TOPPMapAlignerBase::getModelDefaults("interpolated"), true);
  writeParameters("MapRTTransformerModel", TOPPMapAlignerBase::getModelDefaults("none"), true);

  //////////////////////////////////
  // More complicated cases
  //////////////////////////////////

  // ConsensusIDAlgorithm...: abstract base classes, get param. from subclass:
  DOCME2(ConsensusIDAlgorithm, (ConsensusIDAlgorithmBest()));
  DOCME2(ConsensusIDAlgorithmIdentity, (ConsensusIDAlgorithmBest()));
  DOCME2(ConsensusIDAlgorithmSimilarity, (ConsensusIDAlgorithmBest()));
  DOCME2(FeatureFinderAlgorithmIsotopeWavelet, (FeatureFinderAlgorithmIsotopeWavelet()));
  DOCME2(FeatureFinderAlgorithmPicked, (FeatureFinderAlgorithmPicked()));
  DOCME2(FeatureFinderAlgorithmMRM, (FeatureFinderAlgorithmMRM()))
  DOCME2(ProductModel, ProductModel<2>());
  DOCME2(SignalToNoiseEstimatorMeanIterative, SignalToNoiseEstimatorMeanIterative<>());
  DOCME2(SignalToNoiseEstimatorMedian, SignalToNoiseEstimatorMedian<>());
  DOCME2(IonizationSimulation, IonizationSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr()));
  DOCME2(RawMSSignalSimulation, RawMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr()));
  DOCME2(RawTandemMSSignalSimulation, RawTandemMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr()))
  DOCME2(RTSimulation, RTSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr()))
  DOCME2(GaussTraceFitter, (GaussTraceFitter()))
  DOCME2(EGHTraceFitter, (EGHTraceFitter()))

  // handle GUI documentation separately
#ifdef WITH_GUI
  // some classes require a QApplication
  QApplication app(argc, argv);

  DOCME(TOPPViewBase);
  DOCME(TOPPASBase);

  DOCME2(Spectrum1DCanvas, Spectrum1DCanvas(Param(), 0));
  DOCME2(Spectrum2DCanvas, Spectrum2DCanvas(Param(), 0));
  DOCME2(Spectrum3DCanvas, Spectrum3DCanvas(Param(), 0));
#endif

  return 0;
}
