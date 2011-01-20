// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <fstream>
#include <QtGui/QApplication>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>
#include <OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>
#include <OpenMS/ANALYSIS/ID/PILISCrossValidation.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
#include <OpenMS/ANALYSIS/DENOVO/MassDecompositionAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmApplyGivenTrafo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimple.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimplest.h>
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
#include <OpenMS/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>

using namespace std;
using namespace OpenMS;

//**********************************************************************************
//Helper method - use this method to generate the actual parameter documentation
//**********************************************************************************
void writeParameters(const String& class_name, const Param& param)
{
	ofstream f((String("output/OpenMS_") + class_name + ".parameters").c_str());

	f << "<B>Parameters of this class are:</B><BR><BR>\n";
	f << "<table border=\"1\" style=\"border-style:solid; border-collapse:collapse; border-color:#c0c0c0;\" width=\"100%\" cellpadding=\"4\">" << endl;
	f <<"<tr style=\"border-bottom:1px solid black; background:#fffff0\"><th>Name</th><th>Type</th><th>Default</th><th>Restrictions</th><th>Description</th></tr>" << endl;
	String type, description, restrictions;
	for(Param::ParamIterator it = param.begin(); it != param.end();++it)
	{
		restrictions = "";
		if (it->value.valueType()==DataValue::INT_VALUE || it->value.valueType()==DataValue::INT_LIST)
		{
			type = "int";
			if (it->value.valueType()==DataValue::INT_LIST) type += " list";

			//restrictions
			bool first = true;
			if (it->min_int!=-numeric_limits<Int>::max())
			{
				restrictions += String("min: ") + it->min_int;
				first = false;
			}
			if (it->max_int!=numeric_limits<Int>::max())
			{
				if (!first) restrictions += ' ';
				restrictions += String("max: ") + it->max_int;
			}
		}
		else if (it->value.valueType()==DataValue::DOUBLE_VALUE || it->value.valueType()==DataValue::DOUBLE_LIST)
		{
			type = "float";
			if (it->value.valueType()==DataValue::DOUBLE_LIST) type += " list";

			//restrictions
			bool first = true;
			if (it->min_float!=-numeric_limits<DoubleReal>::max())
			{
				restrictions += String("min: ") + it->min_float;
				first = false;
			}
			if (it->max_float!=numeric_limits<DoubleReal>::max())
			{
				if (!first) restrictions += ' ';
				restrictions += String("max: ") + it->max_float;
			}
		}
		else if (it->value.valueType()==DataValue::STRING_VALUE || it->value.valueType()==DataValue::STRING_LIST)
		{
			type = "string";
			if (it->value.valueType()==DataValue::STRING_LIST) type += " list";

			//restrictions
			if (it->valid_strings.size()!=0)
			{
				String valid_strings;
				valid_strings.concatenate(it->valid_strings.begin(),it->valid_strings.end(),", ");
				restrictions += valid_strings;
			}
		}
		if (restrictions=="")
		{
			restrictions="&nbsp;";
		}
		//replace #, @ and newline in description
		description = param.getDescription(it.getName());
		description.substitute("@","XXnot_containedXX");
		description.substitute("XXnot_containedXX","@@");
		description.substitute("#","XXnot_containedXX");
		description.substitute("XXnot_containedXX","@#");
		description.substitute("\n","<BR>");

		//create tooltips for sections if they are documented
		String name = it.getName();
		vector<String> parts;
		name.split(':', parts);
		String prefix = "";
		for (Size i=0; i+1< parts.size(); ++i)
		{
			if (i==0)
			{
				prefix = parts[i];
			}
			else
			{
				prefix = prefix + ":" + parts[i];
			}
			String docu = param.getSectionDescription(prefix);
			if (docu!="")
			{
				parts[i] = String("<span title=\"") + docu + "\">" + parts[i] + "</span>";
			}
		}
		if (parts.size()!=0)
		{
			name.concatenate(parts.begin(), parts.end(), ":");
		}

		//replace # and @ in values
		String value = it->value;
		value.substitute("@","XXnot_containedXX");
		value.substitute("XXnot_containedXX","@@");
		value.substitute("#","XXnot_containedXX");
		value.substitute("XXnot_containedXX","@#");

		//make the advanced parameters cursive, the normal ones bold
		String style = "b";
		if (it->tags.count("advanced")==1) style = "i";

		//final output
		f << "<tr>\n"
      << "  <td style=\"vertical-align:top\"><" << style << ">"<< name << "</" << style << "></td>\n"
      << "  <td style=\"vertical-align:top\">" << type << "</td><td style=\"vertical-align:top\">" << value <<  "</td>\n"
      << "  <td style=\"vertical-align:top\">" << restrictions << "</td><td style=\"vertical-align:top\">" << description <<  "</td>\n"
      << "</tr>\n";
	}
	f << "</table>" << "\n";
	f << "\n" << "<b>Note:</b>" << "\n";
	f << "<UL>" << "\n";
	f << "  <LI> If a section name is documented, the documentation is displayed as tooltip." << "\n";
	f << "  <LI> Advanced parameter names are italic." << "\n";
	f << "</UL>" << "\n";
  f.close();
}

//**********************************************************************************
//Helper macros that can be used for easy classes
//**********************************************************************************

// For classes that have a default-constructor, simply use this macro with the
// class name
#define DOCME(class_name) \
	writeParameters(""#class_name ,class_name().getDefaults());

// For class templates and classes without default constructor use this macro
// with (1.) the class name and (2.) a class instance.
#define DOCME2(class_template_name,instantiation) \
	writeParameters(""#class_template_name,(instantiation).getDefaults());

//**********************************************************************************
//Main method - add your class here
//**********************************************************************************
int main (int argc , char** argv)
{
	// some classes require a QApplication
	QApplication app(argc,argv);

	//////////////////////////////////
	// Simple cases
	//////////////////////////////////

	DOCME(BernNorm);
	DOCME(BiGaussFitter1D);
	DOCME(BiGaussModel);
	DOCME(BinnedSharedPeakCount);
	DOCME(BinnedSpectralContrastAngle);
	DOCME(BinnedSumAgreeingIntensities);
	DOCME(ComplementFilter);
	DOCME(ComplementMarker);
	DOCME(ConsensusID);
	DOCME(DetectabilitySimulation);
	DOCME(DigestSimulation);
	DOCME(EmgFitter1D);
	DOCME(EmgModel);
	DOCME(ExtendedIsotopeFitter1D);
	DOCME(ExtendedIsotopeModel);
	DOCME(FalseDiscoveryRate);
	DOCME(FeatureDeconvolution);
	DOCME(FeatureGroupingAlgorithmLabeled);
	DOCME(FeatureGroupingAlgorithmUnlabeled);
	DOCME(GaussFilter);
	DOCME(GaussFitter1D);
	DOCME(GaussModel);
	DOCME(GoodDiffFilter);
	DOCME(IDDecoyProbability);
	DOCME(IDMapper);
	DOCME(InternalCalibration);
	DOCME(InterpolationModel);
	DOCME(IsotopeDiffFilter);
	DOCME(IsotopeFitter1D);
	DOCME(IsotopeMarker);
	DOCME(IsotopeModel);
	DOCME(ItraqChannelExtractor);
	DOCME(ItraqQuantifier);
	DOCME(LabeledPairFinder);
	DOCME(LinearResampler);
	DOCME(LmaGaussFitter1D);
	DOCME(LmaGaussModel);
	DOCME(LmaIsotopeFitter1D);
	DOCME(LmaIsotopeModel);
	DOCME(MSPFile);
	DOCME(MSSim);
	DOCME(MapAlignmentAlgorithmPoseClustering);
	DOCME(MapAlignmentAlgorithmSpectrumAlignment);
	DOCME(MapAlignmentAlgorithmApplyGivenTrafo);
	DOCME(MapAlignmentAlgorithmIdentification);
	DOCME(NLargest);
	DOCME(NeutralLossDiffFilter);
	DOCME(NeutralLossMarker);
	DOCME(Normalizer);
	DOCME(OptimizePeakDeconvolution);
	DOCME(PILISScoring);
	DOCME(ParentPeakMower);
	DOCME(PeakAlignment);
	DOCME(PeakPickerCWT);
	DOCME(PeakPickerHiRes);
	DOCME(PoseClusteringAffineSuperimposer);
	DOCME(PoseClusteringShiftSuperimposer);
	DOCME(SavitzkyGolayFilter);
	DOCME(SimplePairFinder);
	DOCME(StablePairFinder);
	DOCME(SpectrumAlignment);
	DOCME(SpectrumAlignmentScore);
	DOCME(SpectrumCheapDPCorr);
	DOCME(SpectrumPrecursorComparator);
	DOCME(SteinScottImproveScore);
	DOCME(SpectraMerger);
	DOCME(TICFilter);
	DOCME(TheoreticalSpectrumGenerator);
	DOCME(ThresholdMower);
	DOCME(TwoDOptimization);
	DOCME(WindowMower);
	DOCME(ZhangSimilarityScore);
	DOCME(CompareFouriertransform);
	DOCME(PrecursorIonSelection);
	DOCME(PrecursorIonSelectionPreprocessing);
	DOCME(MorphologicalFilter);
	DOCME(CompNovoIonScoring);
	DOCME(CompNovoIonScoringCID);
	DOCME(CompNovoIdentification);
	DOCME(CompNovoIdentificationCID);
	DOCME(MassDecompositionAlgorithm);
	DOCME(PILISModel);
	DOCME(MRMFragmentSelection);
	DOCME(PILISCrossValidation);
	DOCME(ProtonDistributionModel);
	DOCME(MascotRemoteQuery);
	DOCME(MascotGenericFile);
	DOCME(PILISNeutralLossModel);
	DOCME(PILISModelGenerator);
	DOCME(AdvancedTheoreticalSpectrumGenerator);
	DOCME(FeatureGroupingAlgorithmIdentification);
	DOCME(OfflinePrecursorIonSelection);
	DOCME(TOPPViewBase);
	DOCME(TOPPASBase);
	DOCME(Fitter1D);
	DOCME(EGHModel);
	DOCME(EGHFitter1D);
	DOCME(O18Labeler);
	DOCME(ITRAQLabeler);
	DOCME(Math::PosteriorErrorProbabilityModel);

	//////////////////////////////////
	// More complicated cases
	//////////////////////////////////

	DOCME2(FeatureFinderAlgorithmIsotopeWavelet, (FeatureFinderAlgorithmIsotopeWavelet<Peak1D,Feature>()))
	DOCME2(FeatureFinderAlgorithmPicked, (FeatureFinderAlgorithmPicked<Peak1D,Feature>()));
	DOCME2(FeatureFinderAlgorithmSimple, (FeatureFinderAlgorithmSimple<Peak1D,Feature>()));
	DOCME2(FeatureFinderAlgorithmSimplest, (FeatureFinderAlgorithmSimplest<Peak1D,Feature>()));
	DOCME2(FeatureFinderAlgorithmMRM, (FeatureFinderAlgorithmMRM<Peak1D,Feature>()))
	DOCME2(ModelFitter, (ModelFitter<Peak1D,Feature>(0,0,0)));
	DOCME2(ProductModel,ProductModel<2>());
	DOCME2(SignalToNoiseEstimatorMeanIterative,SignalToNoiseEstimatorMeanIterative<>());
	DOCME2(SignalToNoiseEstimatorMedian,SignalToNoiseEstimatorMedian<>());
	DOCME2(SimpleExtender, (SimpleExtender<Peak1D,Feature>(0,0,0)));
	DOCME2(SimpleSeeder, (SimpleSeeder<Peak1D,Feature>(0,0,0)));
	DOCME2(Spectrum1DCanvas,Spectrum1DCanvas(Param(),0));
	DOCME2(Spectrum2DCanvas,Spectrum2DCanvas(Param(),0));
	DOCME2(Spectrum3DCanvas,Spectrum3DCanvas(Param(),0));
	DOCME2(IonizationSimulation, IonizationSimulation(OpenMS::SimRandomNumberGenerator() ));
	DOCME2(RawMSSignalSimulation, RawMSSignalSimulation(OpenMS::SimRandomNumberGenerator() ));
	DOCME2(RawTandemMSSignalSimulation, RawTandemMSSignalSimulation(OpenMS::SimRandomNumberGenerator() ))
	DOCME2(RTSimulation, RTSimulation(OpenMS::SimRandomNumberGenerator() ))
	DOCME2(TraceFitter,(TraceFitter<Peak1D>()))
	DOCME2(GaussTraceFitter,(GaussTraceFitter<Peak1D>()))
	DOCME2(EGHTraceFitter,(EGHTraceFitter<Peak1D>()))

  return 0;
}
