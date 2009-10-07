// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureFinder FeatureFinder
	
	@brief The feature detection application (quantitation)
	
	This module identifies "features" in a LC/MS map.
	
	By feature, we understand a peptide in a MS sample that
	reveals a characteristic isotope distribution. The algorithm
	computes positions in rt and m/z dimension and a charge estimate
	of each peptide.
	
	The algorithm identifies pronounced regions of the data around so-called <tt>seeds</tt>. 
  In the next step, we iteratively fit a model of the isotope profile and the retention time to
  these data points. Data points with a low probability under this model are removed from the
  feature region. The intensity of the feature is then given by the sum of the data points included
  in its regions.
  
  How to find suitable parameters and details of the different algorithms implemented are described 
	in the TOPP tutorial.
	
	Note that the wavelet transform is very slow on high-resolution spectra (i.e. FT, Orbitrap). We recommend 
	to use a noise or intensity filter to remove spurious points first and to speed-up the feature detection process.
  
	In the following table you, can find example values of the most important parameters for 
	different instrument types. @n These parameters are not valid for all instruments of that type,
	but can be used as a starting point for finding suitable parameters.

	<b>'centroided' algorithm</b>:
	<table>
		<tr>
			<td>&nbsp;</td>
			<td><b>Q-TOF</b></td>
			<td><b>LTQ Orbitrap</b></td>
		</tr>
		<tr>
			<td><b>intensity:bins</b></td>
			<td>10</td>
			<td>10</td>
		</tr>
		<tr>
			<td><b>mass_trace:mz_tolerance</b></td>
			<td>0.02</td>
			<td>0.004</td>
		</tr>
		<tr>
			<td><b>isotopic_pattern:mz_tolerance</b></td>
			<td>0.04</td>
			<td>0.005</td>
		</tr>
	</table>
	
	For the @em centroided algorithm centroided data is needed. In order to create centroided data from profile data use the @ref TOPP_PeakPicker.
	
	Specialized tools are available for some experimental techniques: @ref TOPP_SILACAnalyzer, @ref TOPP_ITRAQAnalyzer.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FeatureFinder.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinder
	: public TOPPBase
{
 public:
	TOPPFeatureFinder()
		: TOPPBase("FeatureFinder","Detects two-dimensional features in LC-MS data.")
	{
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file");
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output file");
		setValidFormats_("out",StringList::create("featureXML"));
		registerInputFile_("seeds","<file>","","User-specified seed list. This feature is not supported by all algorithms!", false);
		setValidFormats_("seeds",StringList::create("featureXML"));
		registerStringOption_("type","<name>","","FeatureFinder algorithm type\n",true);
		setValidStrings_("type", getToolList()[toolName_()] );
		addEmptyLine_();
		addText_("All other options of the Featurefinder depend on the algorithm type used.\n"
						 "They are set in the 'algorithm' section of the INI file.\n");	

		registerSubsection_("algorithm","Algorithm section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		return FeatureFinder().getParameters(type);
	}

	ExitCodes main_(int , const char**)
	{
		//input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");

		Param feafi_param = getParam_().copy("algorithm:",true);

		writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);
				
		String type = getStringOption_("type");
		
		//setup of FeatureFinder
		FeatureFinder ff;
		ff.setLogType(log_type_);
		
		//reading input data
		PeakMap exp;
		MzMLFile f;
		f.setLogType(log_type_);
		PeakFileOptions options;
		
		//load seeds
		FeatureMap<> seeds;
		if (getStringOption_("seeds")!="")
		{
			FeatureXMLFile().load(getStringOption_("seeds"),seeds);
		}
		
		if (type != "mrm")
		{
			//prevent loading of fragment spectra
			options.setMSLevels(vector<Int>(1,1));
			f.getOptions() = options;
		}
		f.load(in, exp);

		//prevent loading of everthing except MRM MS/MS spectra
		if (type == "mrm")
		{
			//exp.erase(remove_if(exp.begin(), exp.end(), HasScanMode<PeakMap::SpectrumType>(InstrumentSettings::SRM, true)), exp.end());
			// erase the spectra, we just need the chromatograms for the feature finder
			exp.erase(exp.begin(), exp.end());
		}
		else
		{
			exp.updateRanges();
		}

		//ouput data
		FeatureMap<> features;

		//running algorithm
		ff.run(type, exp, features, feafi_param, seeds);

    features.applyMemberFunction(&UniqueIdInterface::setUniqueId);

		//-------------------------------------------------------------
		// writing files
		//-------------------------------------------------------------

		//annotate output with data processing info
		addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));

		FeatureXMLFile map_file;
		map_file.store(out,features);			
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPFeatureFinder tool;
	return tool.main(argc,argv);
}

/// @endcond
