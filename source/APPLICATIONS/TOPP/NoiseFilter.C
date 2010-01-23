// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_NoiseFilter NoiseFilter
	
	@brief  Executes a Savitzky Golay or a Gaussian filter to reduce the noise in an MS experiment.
	
	The idea of the Savitzky Golay filter is to find filter-coefficients
	that preserve higher moments, which means to approximate the underlying
	function within the moving window by a polynomial of higher order
	(typically quadratic or quartic) (see A. Savitzky and M. J. E. Golay,
	''Smoothing and Differentiation of Data by Simplified Least Squares Procedures'').
	
	The Gaussian is a peak area preserving low-pass filter and is characterized by narrow bandwidths,
	sharp cutoffs, and low passband ripple.
	
	@note The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the @ref TOPP_Resampler tool).
	      The Gaussian filter works for uniform as well as for non-uniform data.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_NoiseFilter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPNoiseFilter
      : public TOPPBase
{
  public:
    TOPPNoiseFilter()
        : TOPPBase("NoiseFilter","Removes noise from profile spectra by using different smoothing techniques.")
    {
    }

    void registerOptionsAndFlags_()
    {
	  	registerInputFile_("in","<file>","","input raw data file ");
			setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out","<file>","","output raw data file ");
	  	setValidFormats_("out",StringList::create("mzML"));
      registerStringOption_("type","<type>","","smoothing filter type", true);
			setValidStrings_("type", StringList::create("sgolay,gaussian"));
			addEmptyLine_();
	  	addText_("Parameters for the algorithms can be given in the INI file only.");
			addEmptyLine_();
			addText_("Note: The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the Resampler tool).\n"
      				 "      The Gaussian filter works for uniform as well as for non-uniform data.");
    	registerSubsection_("algorithm","Algorithm parameters section");
    }
    
    Param getSubsectionDefaults_(const String& /*section*/) const
    {
			String type = getStringOption_("type");
			Param tmp;
			
			if (type == "sgolay")
      {
        tmp = SavitzkyGolayFilter().getDefaults();
      }
      else if (type == "gaussian")
      {
        tmp = GaussFilter().getDefaults();
      }

      return tmp;
    }

    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      String in = getStringOption_("in");
      String out = getStringOption_("out");
      String type = getStringOption_("type");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------
      MzMLFile mz_data_file;
      mz_data_file.setLogType(log_type_);
      MSExperiment<Peak1D> exp;
      mz_data_file.load(in,exp);

			//check for peak type (profile data required)
			if (PeakTypeEstimator().estimateType(exp[0].begin(),exp[0].end())==SpectrumSettings::PEAKS)
			{
				writeLog_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
			}
			
			//check if spectra are sorted
			for (Size i=0; i< exp.size(); ++i)
			{
				if (!exp[i].isSorted())
				{
					writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
					return INCOMPATIBLE_INPUT_DATA;
				}
			}

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
    	Param filter_param = getParam_().copy("algorithm:",true);
			writeDebug_("Parameters passed to filter", filter_param,3);
      if (type == "sgolay")
      {	
  			SavitzkyGolayFilter sgolay;
        sgolay.setLogType(log_type_);
  			sgolay.setParameters( filter_param );
				sgolay.filterExperiment(exp);
      }
      else if (type == "gaussian")
      {	
        GaussFilter gauss;
        gauss.setLogType(log_type_);
        gauss.setParameters(filter_param);
        try
        {
          gauss.filterExperiment(exp);
        }
        catch(Exception::IllegalArgument& e)
        {
        	writeLog_(String("Error: ") + e.getMessage()) ;
        	return INCOMPATIBLE_INPUT_DATA;
        }
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      mz_data_file.store(out,exp);

			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::SMOOTHING));

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPNoiseFilter tool;
  return tool.main(argc,argv);
}

/// @endcond
