// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FILTERING/BASELINE/TopHatFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page BaselineFilter BaselineFilter

   @brief Executes the top-hat filter to remove the baseline of an MS experiment.

   This nonlinear filter, known as the top-hat operator in morphological mathematics
   (see Soille, ''Morphological Image Analysis''), is independent of the underlying baseline shape.
   It is able to detect an over brightness even if the environment is not uniform.
   The principle is based on the subtraction of a signal from its opening (erosion followed by a dilation).
   The size the structuring element (here a flat line) being conditioned by the width of the lineament
   (in our case the maximum width of a mass spectrometric peak) to be detected.

   @note The length (given in Thomson) of the structuring element should be wider than the
	 maximum peak width in the raw data.

   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBaselineFilter
	: public TOPPBase
{
 public:
	TOPPBaselineFilter()
		: TOPPBase("BaselineFilter","Removes the baseline from profile spectra using a top-hat filter.")
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
	  	registerInputFile_("in","<file>","","input raw data file ");
			setValidFormats_("in",StringList::create("mzData"));
			registerOutputFile_("out","<file>","","output raw data file ");
	  	setValidFormats_("out",StringList::create("mzData"));
      registerDoubleOption_("struc_elem_length","<size>",2.5,"Length of the structuring element in Th.",false);
      registerDoubleOption_("resampling","<spacing>",0.0,"Spacing for the resampling process.",false);
      addEmptyLine_();
			addText_("Note: The top-hat filter works only on uniform data (to generate equally spaced data you have to set the resampling option!)");
	}
  
 	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		double spacing = getDoubleOption_("resampling");

		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------

		MzDataFile mz_data_file;
		MSExperiment<Peak1D > exp;
		mz_data_file.setLogType(log_type_);
		mz_data_file.load(in,exp);

		//check for peak type (raw data required)
		if (exp.getProcessingMethod().getSpectrumType()==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: The file meta data claims that this is not raw data!");
		}
		if (PeakTypeEstimator().estimateType(exp[0].begin(),exp[0].end())==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: OpenMS peak type estimation indicates that this is not raw data!");
		}

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		TopHatFilter tophat;
    tophat.setLogType(log_type_);
    Param tophat_param;
    tophat_param.setValue("struc_elem_length",getDoubleOption_("struc_elem_length"));
		tophat.setParameters(tophat_param);

		// no resampling of the data
		if (spacing==0.0)
		{
			tophat.filterExperiment(exp);
		}
		else
		{
			LinearResampler lin_resampler;
			lin_resampler.setLogType(log_type_);
			Param resampler_param;
			resampler_param.setValue("spacing",spacing);
			lin_resampler.setParameters(resampler_param);
		
      tophat.startProgress(0,exp.size(),"resampling and baseline filtering of data");
			// resample and filter every scan
			for (UInt i = 0; i < exp.size(); ++i)
			{
				// temporary container for the resampled data
				MSSpectrum<Peak1D> resampled_data;
				lin_resampler.raster(exp[i],resampled_data);

				MSSpectrum<Peak1D> spectrum;
				tophat.filter(resampled_data, spectrum);
        
        exp[i].getContainer() = spectrum.getContainer();
				tophat.setProgress(i);
			}
      tophat.endProgress();
		}
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		exp.getProcessingMethod().setSpectrumType(SpectrumSettings::RAWDATA);
		mz_data_file.store(out,exp);

		return EXECUTION_OK;
	}

};




int main( int argc, const char** argv )
{
    TOPPBaselineFilter tool;
    return tool.main(argc,argv);
}

/// @endcond
