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
#include <OpenMS/config.h>

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolaySVDFilter.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

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
 
   @note Use a Gaussian filter  which has approximately the same width as your mass peaks.
 
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPNoiseFilter
      : public TOPPBase
{
  public:
    TOPPNoiseFilter()
        : TOPPBase("NoiseFilter","remove the noise from LC/MS raw data")
    {
    }

    void registerOptionsAndFlags_()
    {
	  	registerStringOption_("in","<file>","","input mzData file (raw data)");
			registerStringOption_("out","<file>","","output mzData file (raw data)");
      registerStringOption_("filter_type","<type>","","smoothing filter type. Valid types are: 'sgolay' or 'gaussian'");
      registerDoubleOption_("resampling","<spacing>",0.0,"spacing for the resampling process",false);
			addEmptyLine_();
	  	addText_("Parameters for the algorithms can be given in the INI file only:\n"
								"<NODE name=\"sgolay\">\n"
								"  <ITEM name=\"frame_length\" value=\"9\" type=\"int\"/>\n"
								"  <ITEM name=\"polynomial_order\" value=\"4\" type=\"int\"/>\n"
								"</NODE>\n"
								"<NODE name=\"gaussian\">\n"
								"  <ITEM name=\"gaussian_width\" value=\"0.8\" type=\"float\" />\n"
								"</NODE>");
			addEmptyLine_();
			addText_("Note: The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the resampling option).\n"
      				 "      The Gaussian filter works for uniform as well as for non-uniform data.");
    	registerSubsection_("sgolay");
    	registerSubsection_("gaussian");
    }

    ExitCodes main_(int , char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      String in = getStringOption_("in");
      String out = getStringOption_("out");
      String filter_type = getStringOption_("filter_type");
      float spacing = getDoubleOption_("resampling");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MzDataFile mz_data_file;
      MSExperiment<DRawDataPoint<1> > ms_exp_raw;
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
      MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

      if (filter_type == "sgolay")
      {	
      	Param filter_param = getParam_().copy("sgolay:",true);
  			writeDebug_("Parameters passed to SavitzkyGolaySVDFilter", filter_param,3);
  			SavitzkyGolaySVDFilter sgolay;
  			sgolay.setParam( filter_param );
        
        LinearResampler lin_resampler;
        lin_resampler.setSpacing(spacing);

        // copy the experimental settings
        static_cast<ExperimentalSettings&>(ms_exp_filtered) = ms_exp_raw;

        // no resampling of the data
        if (spacing==0.0)
        { 
           sgolay.filterExperiment(ms_exp_raw,ms_exp_filtered);
					 writeDebug_(String("No resampling!"), 1);
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
						
						if (resampled_data.size() == 1)
						{
							ms_exp_filtered.push_back(resampled_data);
						}
						else
						{
							sgolay.filter(resampled_data, spectrum);
						}

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
      }
      else if (filter_type == "gaussian")
      {	
      	Param filter_param = getParam_().copy("gaussian:",true);
  			writeDebug_("Parameters passed to GaussFilter", filter_param,3);
        GaussFilter gauss;
        gauss.setParam(filter_param);
        gauss.filterExperiment(ms_exp_raw, ms_exp_filtered);
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
  TOPPNoiseFilter tool;
  return tool.main(argc,argv);
}

/// @endcond
