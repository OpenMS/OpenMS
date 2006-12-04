// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

	 @todo Parameter struc_elem_length does not work as expected (Eva)

   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBaselineFilter
	: public TOPPBase
{
 public:
	TOPPBaselineFilter()
		: TOPPBase("BaselineFilter")
	{}

 protected:
	void printToolUsage_() const
	{
		cerr << endl
				 << getToolName() << " -- remove the baseline in a LC/MS experiment" << endl
				 << "Version: " << VersionInfo::getVersion() << endl
				 << endl
				 << "This application implements a morpholgical approach called top-hat filter." << endl
				 << endl
				 << "Note: The top-hat filter works only on uniform data (to generate equally spaced data you have to set the resampling flag!)." << endl
				 << endl
				 << "Usage:" << endl
				 << " " << getToolName() << " [options]" << endl
				 << "  -in <file>                input mzData file name" << endl
				 << "  -out <file>               output mzData file name" << endl
				 << "  -struc_elem_length <size> length of the structuring element (default is 2.5Th)" << endl
				 << "  -resampling <spacing>     spacing for the resampling process (default: deactivated)"<< endl
				 << endl;
	}

	void printToolHelpOpt_() const
	{
		cerr << endl
				 << getToolName() << endl
				 << endl
				 << "INI options:" << endl
				 << "  in <file>                input mzData file name" << endl
				 << "  out <file>               output mzData file name" << endl
				 << "  struc_elem_length <size> length of the structuring element (given in Th)" << endl
				 << "  resampling <spacing>     spacing for the resampling process (default: this flag is not set)" << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"struc_elem_length\" value=\"2.0\" type=\"float\"/>" << endl
				 << "  <ITEM name=\"resampling\" value=\"0.05\" type=\"float\"/>" << endl;
	}

	void setOptionsAndFlags_()
	{
		options_["-out"] = "out";
		options_["-in"] = "in";
		options_["-struc_elem_length"] = "struc_elem_length";
		options_["-resampling"] = "resampling";
	}

	ExitCodes main_(int , char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		//input file names and types
		String in = getParamAsString_("in");
		writeDebug_(String("Input file: ") + in, 1);

		//output file names and types
		String out = getParamAsString_("out");
		writeDebug_(String("Output file: ") + out, 1);

		//length of the structuring element
		String struc_elem = getParamAsString_("struc_elem_length");
		writeDebug_(String("Length of structuring element: ") + struc_elem, 1);

		//spacing for resampling process
		String resampling = getParamAsString_("resampling");
		writeDebug_(String("Resampling: ") + resampling, 1);


		float struc_elem_length = 2.5;
		float spacing = 0.;
		bool resampling_flag = false;
		try
		{
			//resampling
			if (struc_elem != "")
			{
				struc_elem_length = struc_elem.toFloat();
			}
		}
		catch(Exception::ConversionError& e)
		{
			writeLog_(String("Invalid length for the structuring element '") + struc_elem  + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		try
		{
			//resampling
			if (resampling != "")
			{
				spacing = resampling.toFloat();
				resampling = true;
			}
		}
		catch(Exception::ConversionError& e)
		{
			writeLog_(String("Invalid spacing '") + resampling + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------

		MzDataFile mz_data_file;
		MSExperiment<DRawDataPoint<1> > ms_exp_raw;
		MSExperiment<DRawDataPoint<1> > ms_exp_filtered;
		mz_data_file.load(in,ms_exp_raw);

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
		if (!resampling_flag)
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
