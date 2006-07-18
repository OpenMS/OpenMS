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
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/SMOOTHING/DSavitzkyGolaySVDFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/DGaussFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include "TOPPBase.h"

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
        (The default value is 0.8 Th.)

   @ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES

class TOPPNoiseFilter
      : public TOPPBase
{
public:
  TOPPNoiseFilter()
      : TOPPBase("NoiseFilter")
  {}

protected:
  void printToolUsage_()
  {
    cerr << endl
    << tool_name_ << " -- remove the noise in a LC/MS experiment" << endl
    << endl
    << "This application implements a smoothing filter. It executes a Savitzky Golay or alternatively a Gaussian filter." << endl
    << endl
    << "Note: The Savitzky Golay filter works only on uniform data (to generate equally spaced data use the resampling option)." << endl
    << "      The Gaussian filter works for uniform as well as for non-uniform data." << endl
    << endl
    << "Usage:" << endl
    << " " << tool_name_ << " [options]" << endl
    << endl
    << "Options are:" << endl
    << "  -filter_type <type>   smoothing filter type. Valid filter options are: 'sgolay' or 'gaussian'." << endl
    << "  -resampling <spacing> spacing for the resampling process (default: this flag is not set)" << endl
    << "  -in <file>            input mzData file name" << endl
    << "  -out <file>           output mzData file name" << endl
    << endl;
  }

  void printToolHelpOpt_()
  {
    cerr << endl
    << tool_name_ << endl
    << endl
    << "INI options:" << endl
    << "  in <file>            input mzData file name" << endl
    << "  out <file>           output mzData file name" << endl
    << "  filter_type <type>   smoothing filter type. Valid filter options are: 'sgolay' or 'gaussian'." << endl
    << "  resampling <spacing> spacing for the resampling process (default: deactivated)" << endl
    << endl
    << "INI File example section:" << endl
    << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"filter_type\" value=\"gaussian\" type=\"string\"/>" << endl
    << "  <ITEM name=\"resampling\" value=\"0.05\" type=\"float\"/>" << endl;
  }

  void setOptionsAndFlags_()
  {
    options_["-out"] = "out";
    options_["-in"] = "in";
    options_["-filter_type"] = "filter_type";
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

    //filter type
    String filter_type = getParamAsString_("filter_type");
    writeDebug_(String("Filter type: ") + filter_type, 1);

    //spacing for resampling process
    String resampling = getParamAsString_("resampling");
    writeDebug_(String("Resampling: ") + resampling, 1);

    float spacing = 0.;
    bool resampling_flag = false;

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
    mz_data_file.load(in,ms_exp_raw);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    MSExperiment<DRawDataPoint<1> > ms_exp_filtered;

    if (filter_type == "sgolay")
    {
      DSavitzkyGolaySVDFilter<1> sgolay(param_);

      if (!resampling_flag)
      {
        ms_exp_raw >> sgolay(ms_exp_filtered);
      }
      else
      {
        LinearResampler<DRawDataPoint<1> > lin_resampler;
        lin_resampler.setSpacing(spacing);

        MSExperiment<DRawDataPoint<1> >::const_iterator first_scan = ms_exp_raw.begin();
        MSExperiment<DRawDataPoint<1> >::const_iterator last_scan = ms_exp_raw.end();

        // copy the experimental settings
        static_cast<ExperimentalSettings>(ms_exp_filtered) = ms_exp_raw;

        while (first_scan != last_scan)
        {
          MSSpectrum<DRawDataPoint<1> >::const_iterator first_data_point = first_scan->begin();
          MSSpectrum<DRawDataPoint<1> >::const_iterator last_data_point = first_scan->end();

          // create the resampled spectrum
          int number_resampled_points = (int)(ceil(((last_data_point-1)->getPos() - first_data_point->getPos()) / spacing + 1));
          DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > resampled_data(number_resampled_points);

          lin_resampler.start(first_data_point,last_data_point,resampled_data.begin());

          // create the baseline filtered spectrum
          DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > filtered_data(number_resampled_points);
          sgolay.filter(resampled_data.begin(),resampled_data.end(),filtered_data.begin());

          MSSpectrum<DRawDataPoint<1> > spectrum;
          spectrum.setContainer(filtered_data);

          spectrum.setRetentionTime(first_scan->getRetentionTime(), first_scan->getRetentionTimeStart(), first_scan->getRetentionTimeStop());
          spectrum.setMSLevel(first_scan->getMSLevel());
          spectrum.setName(first_scan->getName());

          ms_exp_filtered.std::vector< MSSpectrum< DRawDataPoint<1> > >::push_back(spectrum);
          ++first_scan;
        }
      }
    }
    else
      if (filter_type == "gaussian")
      {
        DGaussFilter<1> gauss(param_);
        ms_exp_raw >> gauss(ms_exp_filtered);
      }
      else
      {
        writeLog_(String("Unknown filter type '") + filter_type + "' given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }



    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    mz_data_file.store(out,ms_exp_filtered);

    return OK;
  }
};

/// @endcond

int main( int argc, char ** argv )
{
  TOPPNoiseFilter tool;
  return tool.main(argc,argv);
}


