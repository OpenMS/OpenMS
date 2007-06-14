// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <fstream>


// Documentation is in .C files:
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>
#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SumReducer.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayQRFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
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
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummyFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummySeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MarrWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PickedPeakSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>


// Documentation is in .h files:
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDecharger.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>
#include <OpenMS/FILTERING/BASELINE/MorphFilter.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>


using namespace std;
using namespace OpenMS;

//**********************************************************************************
//Helper method - do not change anything here
//**********************************************************************************
void writeParameters(std::ofstream& f, const String& class_name, const Param& param)
{
	f <<
		"/**\n"
		" @page " << class_name << "_Parameters " << class_name << " Parameters" << endl;

	// Manually create a link to the class documentation.  Doxygen just won't do this with @link or @sa.
	String class_doc = "./class";
	if ( !class_name.hasPrefix("OpenMS::") ) class_doc += "OpenMS::";
	class_doc += class_name;
	class_doc.substitute("::","_1_1");
	f << "<a href=\"" << class_doc << ".html\">" << class_name << "</a>\n";

	String type, description;
	for(map<String,DataValue>::const_iterator it = param.begin(); it != param.end();++it)
	{
		if (it->second.valueType()==DataValue::INTVALUE || it->second.valueType()==DataValue::LONVALUE || it->second.valueType()==DataValue::SHOVALUE  )
		{
			type = "int";
		}
		if (it->second.valueType()==DataValue::FLOVALUE || it->second.valueType()==DataValue::DOUVALUE )
		{
			type = "float";
		}
		if (it->second.valueType()==DataValue::STRVALUE )
		{
			type = "string";
		}
		description = param.getDescription(it->first);
		description.substitute("\n","@n ");
		f <<" - @b "<< it->first << " (" << type << "): " << description << endl;
	}
	f << "*/" << endl;
	f << endl;
}

//**********************************************************************************
//Main method - add your class here
//**********************************************************************************
int main (int, char**)
{
	ofstream f;
	f.open("DefaultParameters.doxygen");

	// for straight classes
#define DOCME(class) \
writeParameters(f,""#class ,class().getParameters());

	// For class templates (you may have to put parens around the instantiation)
#define DOCME2(class_template_name,instantiation) \
writeParameters(f,""#class_template_name,instantiation.getParameters());

	//////////////////////////////////
	//
	// Documentation is in .C files:
	//
	//////////////////////////////////

	//	DOCME(BaseSweepSeeder); // cannot instantiate (pure virtual functions)
	DOCME(BiGaussModel);
	// DOCME(BinnedRepCompareFunctor); // cannot instantiate (pure virtual functions)
	DOCME(BinnedRepMutualInformation);
	DOCME(ComplementFilter);
	DOCME(ComplementMarker);
	DOCME(ConsensusID);
	DOCME(DummyFitter);
	DOCME(DummySeeder);
	DOCME(EmgModel);
	DOCME(ExtendedModelFitter);
	DOCME(GaussModel);
	DOCME(GoodDiffFilter);
	DOCME(IsotopeDiffFilter);
	DOCME(IsotopeMarker);
	DOCME(IsotopeModel);
	DOCME(IsotopeWaveletSeeder);
	DOCME(LmaGaussModel);
	DOCME(LogNormalModel);
	DOCME(MarrWaveletSeeder);
	DOCME(MaxReducer);
	DOCME(NLargest);
	DOCME(NeutralLossDiffFilter);
	DOCME(NeutralLossMarker);
	DOCME(Normalizer);
	DOCME(OptimizePeakDeconvolution);
	DOCME(PILISIdentification);
	DOCME(PILISModel);
	DOCME(PILISModelGenerator);
	DOCME(PILISScoring);
	// DOCME(PairMatcher);  // cannot instantiate (no default constructor)
	DOCME(ParentPeakMower);
	DOCME(PeakPicker);
	DOCME(PeakPickerCWT);
	// DOCME(PeakSpectrumCompareFunctor); // cannot instantiate (pure virtual functions)
	DOCME(PickedPeakSeeder);
	DOCME(ProtonDistributionModel);
	DOCME(SavitzkyGolayQRFilter);
	DOCME(SavitzkyGolaySVDFilter);
	DOCME(SimpleExtender);
	DOCME(SimpleModelFitter);
	DOCME(SimpleSeeder);
	// DOCME(Spectrum1DCanvas);  // cannot instantiate (no default constructor)
	// DOCME(Spectrum2DCanvas);  // cannot instantiate (no default constructor)
	// DOCME(Spectrum3DCanvas);  // cannot instantiate (no default constructor)
	DOCME(SpectrumAlignment);
	DOCME(SpectrumAlignmentScore);
	DOCME(SpectrumCheapDPCorr);
	DOCME(SpectrumPrecursorComparator);
	DOCME(SumReducer);
	DOCME(TICFilter);
	// DOCME(TOPPViewBase);  // runtime error (no QApplication instantiated)
	DOCME(TheoreticalSpectrumGenerator);
	DOCME(ThresholdMower);
	DOCME(TwoDOptimization);
	DOCME(WindowMower);
	DOCME(ZhangSimilarityScore);

	//////////////////////////////////
	//
	// Documentation is in .h files:
	//
	//////////////////////////////////

	// DOCME2(BaseAlignment);
	// DOCME2(BaseModel,BaseModel<1>()); // cannot instantiate (pure virtual functions)
	// DOCME2(BasePairwiseMapMatcher);
	// DOCME2(DelaunayPairFinder);
	// DOCME2(FeatureDecharger);
	// DOCME2(GaussFilter);
	// DOCME2(HierarchicalClustering);
	// DOCME2(InterpolationModel,InterpolationModel<>()); // cannot instantiate (pure virtual functions)
	// DOCME2(LinearResampler);
	// DOCME2(MorphFilter);
	// DOCME2(PoseClusteringAffineSuperimposer);
	// DOCME2(PoseClusteringPairwiseMapMatcher);
	// DOCME2(PoseClusteringShiftSuperimposer);
	DOCME2(ProductModel,ProductModel<2>()); // YEAH!!!
	// DOCME2(SignalToNoiseEstimatorMeanIterative);
	// DOCME2(SignalToNoiseEstimatorMedian);
	// DOCME2(SimplePairFinder);

  return 0;
}


