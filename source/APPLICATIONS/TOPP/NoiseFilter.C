// -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/config.h>

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page NoiseFilter NoiseFilter
 
   @brief  Executes a Savitzky Golay or a Gaussian filter to reduce the noise in a MS experiment.
 
   The idea of the Savitzky Golay filter is to find filter-coefficients
   that preserve higher moments, which means to approximate the underlying
   function within the moving window by a polynomial of higher order
   (typically quadratic or quartic) (see A. Savitzky and M. J. E. Golay,
   ''Smoothing and Differentiation of Data by Simplified Least Squares Procedures'').
   
   The Gaussian is a peak area preserving low-pass filter and is characterized by narrow bandwidths,
   sharp cutoffs, and low passband ripple.
 
   @ingroup TOPP
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
			setValidFormats_("in",StringList::create("mzData"));
			registerOutputFile_("out","<file>","","output raw data file ");
	  	setValidFormats_("out",StringList::create("mzData"));
      registerStringOption_("type","<type>","","smoothing filter type", true);
			setValidStrings_("type", StringList::create("sgolay,gaussian"));
      registerDoubleOption_("resampling","<spacing>",0.0,"spacing for the resampling process",false);
			addEmptyLine_();
	  	addText_("Parameters for the algorithms can be given in the INI file only.");
			addEmptyLine_();
			addText_("Note: The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the resampling option).\n"
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
      float spacing = getDoubleOption_("resampling");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MzDataFile mz_data_file;
      mz_data_file.setLogType(log_type_);
      MSExperiment<Peak1D > ms_exp_raw;
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
      MSExperiment<Peak1D > ms_exp_filtered;

    	Param filter_param = getParam_().copy("algorithm:",true);
			writeDebug_("Parameters passed to filter", filter_param,3);
      if (type == "sgolay")
      {	
  			SavitzkyGolayFilter sgolay;
        sgolay.setLogType(log_type_);
  			sgolay.setParameters( filter_param );
        
        // no resampling of the data
        if (spacing==0.0)
        { 
           sgolay.filterExperiment(ms_exp_raw,ms_exp_filtered);
					 writeDebug_(String("No resampling!"), 1);
        }
        else
        {
					LinearResampler lin_resampler;
					lin_resampler.setLogType(log_type_);
					Param resampler_param;
					resampler_param.setValue("spacing",spacing);
					lin_resampler.setParameters(resampler_param);
			
          UInt n = ms_exp_raw.size();
          sgolay.startProgress(0,n,"smoothing mzData file");
          lin_resampler.startProgress(0,n,"resampling of data");
          // resample and filter every scan
          for (UInt i = 0; i < n; ++i)
          {
            // temporary container for the resampled data
            MSSpectrum<Peak1D> resampled_data;
            lin_resampler.raster(ms_exp_raw[i],resampled_data);
            lin_resampler.setProgress(i);

            MSSpectrum<Peak1D> spectrum;
						
						if (resampled_data.size() == 1)
						{
							ms_exp_filtered.push_back(resampled_data);
						}
						else
						{
							sgolay.filter(resampled_data, spectrum);
              sgolay.setProgress(i);
						}

            // if any peaks are found copy the spectrum settings
            if (spectrum.size() > 0)
            {
              // copy the spectrum settings
              static_cast<SpectrumSettings&>(spectrum) = ms_exp_raw[i];
              spectrum.setType(SpectrumSettings::RAWDATA);

              // copy the spectrum information
              spectrum.getPrecursorPeak() = ms_exp_raw[i].getPrecursorPeak();
              spectrum.setRT(ms_exp_raw[i].getRT());
              spectrum.setMSLevel(ms_exp_raw[i].getMSLevel());
              spectrum.getName() = ms_exp_raw[i].getName();

              ms_exp_filtered.push_back(spectrum);
            }
          }
          sgolay.endProgress();
          lin_resampler.endProgress();
        }
      }
      else if (type == "gaussian")
      {	
        GaussFilter gauss;
        gauss.setLogType(log_type_);
        gauss.setParameters(filter_param);
        try
        {
          gauss.filterExperiment(ms_exp_raw, ms_exp_filtered);
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
			
			ms_exp_filtered.getProcessingMethod().setSpectrumType(SpectrumSettings::RAWDATA);
      mz_data_file.store(out,ms_exp_filtered);

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPNoiseFilter tool;
  return tool.main(argc,argv);
}

/// @endcond
