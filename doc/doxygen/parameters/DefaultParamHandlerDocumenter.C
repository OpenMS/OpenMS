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
#include <OpenMS/KERNEL/FeatureMap.h>

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
//Helper method - use this method to generate the actual paramter documentation
//**********************************************************************************
void writeParameters(std::ofstream& f, const String& class_name, const Param& param)
{
	static vector<String> class_names;
	if (class_name=="CREATE_MAIN_PAGE")
	{
		//sort classes
		sort(class_names.begin(), class_names.end());
		f << "/**" << endl;
		f << " @page Parameters_Main_Page Class parameters summary" << endl;
		
		for (vector<String>::const_iterator it = class_names.begin(); it!= class_names.end(); ++it)
		{
			f << " - @subpage " << *it << "_Parameters" << endl;
		}
		
		f << "*/" << endl;
		f << endl;
	}
	else
	{
		class_names.push_back(class_name);

		f << "/**" << endl;
		f << " @page " << class_name << "_Parameters " << class_name << " Parameters" << endl;
	
		// Manually create a link to the class documentation.  Doxygen just won't do this with @link or @sa.
		String class_doc = "./class";
		if ( !class_name.hasPrefix("OpenMS::") ) class_doc += "OpenMS::";
		class_doc += class_name;
		class_doc.substitute("::","_1_1");
		f << "Parameters of <a href=\"" << class_doc << ".html\">" << class_name << "</a>:<BR><BR>\n";
		f << "<table border=1>" << endl;
		f <<"<tr><th>Name</th><th>Type</th><th>Default</th><th>Description</th></tr>" << endl;
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
			f <<"<tr><td><b>"<< it->first << "</b></td><td>" << type << "</td><td>" << it->second.toString() <<  "</td><td>" << description <<  "</td></tr>" << endl;
		}
		f << "</table>" << endl;
		f << "*/" << endl;
		f << endl;
	}
}

//**********************************************************************************
//Helper macros that can be used for easy classes
//**********************************************************************************

// for classes that have a default-constructor, simply use this macro with the class name
#define DOCME(class) \
writeParameters(f,""#class ,class().getParameters());

// For class templates and classes without default constructor use this macro with the macro and a class instance
// you may have to put parenteses around the instantiation
#define DOCME2(class_template_name,instantiation) \
writeParameters(f,""#class_template_name,instantiation.getParameters());

//**********************************************************************************
//Main method - add your class here
//**********************************************************************************
int main (int argc , char** argv)
{
	//some classes require a QApplication
	QApplication app(argc,argv);
	
	ofstream f;
	f.open("DefaultParameters.doxygen");

	//////////////////////////////////
	//
	// Documentation is in .C files:
	//
	//////////////////////////////////

	DOCME(BiGaussModel);
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
	
	//PairMatcher
	FeatureMap<> features;
	PairMatcher pm(features);
	writeParameters(f,"PairMatcher",pm.getParameters());
	
	DOCME(ParentPeakMower);
	DOCME(PeakPicker);
	DOCME(PeakPickerCWT);
	DOCME(PickedPeakSeeder);
	DOCME(ProtonDistributionModel);
	DOCME(SavitzkyGolayQRFilter);
	DOCME(SavitzkyGolaySVDFilter);
	DOCME(SimpleExtender);
	DOCME(SimpleModelFitter);
	DOCME(SimpleSeeder);
	DOCME2(Spectrum1DCanvas,Spectrum1DCanvas(Param(),0));
	DOCME2(Spectrum2DCanvas,Spectrum2DCanvas(Param(),0));
	DOCME2(Spectrum3DCanvas,Spectrum3DCanvas(Param(),0));
	DOCME(SpectrumAlignment);
	DOCME(SpectrumAlignmentScore);
	DOCME(SpectrumCheapDPCorr);
	DOCME(SpectrumPrecursorComparator);
	DOCME(SumReducer);
	DOCME(TICFilter);
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
	// DOCME2(BasePairwiseMapMatcher);
	// DOCME2(DelaunayPairFinder);
	// DOCME2(FeatureDecharger);
	// DOCME2(GaussFilter);
	// DOCME2(HierarchicalClustering);
	// DOCME2(LinearResampler);
	// DOCME2(MorphFilter);
	// DOCME2(PoseClusteringAffineSuperimposer);
	// DOCME2(PoseClusteringPairwiseMapMatcher);
	// DOCME2(PoseClusteringShiftSuperimposer);
	DOCME2(ProductModel,ProductModel<2>()); // YEAH!!!
	// DOCME2(SignalToNoiseEstimatorMeanIterative);
	// DOCME2(SignalToNoiseEstimatorMedian);
	// DOCME2(SimplePairFinder);

	//create main page for all parameter documentations
	writeParameters(f,"CREATE_MAIN_PAGE",Param());

  return 0;
}
