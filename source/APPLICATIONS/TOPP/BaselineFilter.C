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
   The principle is based on the subtraction of an signal from its opening (erosion followed by a dilation).
   The size the structuring element (here a flat line) being conditioned by the width of the lineament
   (in our case the maximum width of a mass spectrometric peak) to be detected.

   @note The length (given in Thomson) of the structuring element should be wider than the
	 maximal peak width in the raw data.

   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBaselineFilter
	: public TOPPBase
{
 public:
	TOPPBaselineFilter()
		: TOPPBase("BaselineFilter","top-hat filter for baseline reduction")
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
	  	registerStringOption_("in","<file>","","input mzData file (raw data)");
			registerStringOption_("out","<file>","","output mzData file (raw data)");
      registerDoubleOption_("struc_elem_length","<size>",2.5,"length of the structuring element in Th",false);
      registerDoubleOption_("resampling","<spacing>",0.0,"spacing for the resampling process",false);
      addEmptyLine_();
			addText_("Note: The top-hat filter works only on uniform data (to generate equally spaced data you have to set the resampling option!)");
	}

	ExitCodes main_(int , char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		double struc_elem_length = getDoubleOption_("struc_elem_length");
		double spacing = getDoubleOption_("resampling");

		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------

		MzDataFile mz_data_file;
		MSExperiment<DRawDataPoint<1> > ms_exp_raw;
		MSExperiment<DRawDataPoint<1> > ms_exp_filtered;
		mz_data_file.load(in,ms_exp_raw);

		//check for peak type (raw data required)
		if (ms_exp_raw.getProcessingMethod().getSpectrumType()==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: The file meta data claims that this is not raw data!");
		}
		if (PeakTypeEstimator().estimateType(ms_exp_raw[0].begin(),ms_exp_raw[0].end())==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: OpenMS peak type estimation indicates that this is not raw data!");
		}

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		TopHatFilter tophat;
		tophat.setStrucElemSize(struc_elem_length);

		LinearResampler lin_resampler;
		lin_resampler.setSpacing(spacing);

		// copy the experimental settings
		static_cast<ExperimentalSettings&>(ms_exp_filtered) = ms_exp_raw;

		// no resampling of the data
		if (spacing==0.0)
		{
			tophat.filterExperiment(ms_exp_raw,ms_exp_filtered);
		}
		else
		{
			unsigned int n = ms_exp_raw.size();
			// resample and filter every scan
			for (unsigned int i = 0; i < n; ++i)
			{
				// temporary container for the resampled data
				MSSpectrum< DRawDataPoint<1> > resampled_data;
				lin_resampler.raster(ms_exp_raw[i],resampled_data);

				MSSpectrum< DRawDataPoint<1> > spectrum;
				tophat.filter(resampled_data, spectrum);

				// if any peaks are found copy the spectrum settings
				if (spectrum.size() > 0)
				{
					// copy the spectrum settings
					static_cast<SpectrumSettings&>(spectrum) = ms_exp_raw[i];
					spectrum.setType(SpectrumSettings::RAWDATA);

					// copy the spectrum information
					spectrum.getPrecursorPeak() = ms_exp_raw[i].getPrecursorPeak();
					spectrum.setRetentionTime(ms_exp_raw[i].getRetentionTime());
					spectrum.setMSLevel(ms_exp_raw[i].getMSLevel());
					spectrum.getName() = ms_exp_raw[i].getName();

					ms_exp_filtered.push_back(spectrum);
				}
			}
		}
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		ms_exp_filtered.getProcessingMethod().setSpectrumType(SpectrumSettings::RAWDATA);
		mz_data_file.store(out,ms_exp_filtered);

		return EXECUTION_OK;
	}
};

int main( int argc, char ** argv )
{
    TOPPBaselineFilter tool;
    return tool.main(argc,argv);
}

/// @endcond
