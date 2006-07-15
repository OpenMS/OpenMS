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
#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include "TOPPBase.h"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page PeakPicker PeakPicker

   @brief Executes the peak picking algorithm
   as described by Lange et al. (2006) Proc. PSB-06.

   The conversion of the ''raw'' ion count data acquired
   by the machine into peak lists for further processing
   is usually called peak picking. Our algorithm is independent
   of the underlying machine or ionization method, and is able
   to resolve highly convoluted and asymmetric signals.
   The method uses the multiscale nature of spectrometric data by
   first detecting the mass peaks in the wavelet-transformed signal
   before a given asymmetric peak function is fitted to the raw data.
   In an optional third stage, the resulting fit can be further improved using
   techniques from nonlinear optimization.

   @ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES 

class TOPPPeakPicker
      : public TOPPBase
{
public:
  TOPPPeakPicker()
      : TOPPBase("PeakPicker")
  {}

protected:
  void printToolUsage_()
  {

    cerr << endl
    << tool_name_ << " -- find mass spectrometric peaks in LC/MC experiments." << endl
    << "This application implements an algorithm for peak picking as " << endl
    << "described in Lange et al. (2006) Proc. PSB-06. "<< endl
    << endl
    << "Usage:" << endl
    << " " << tool_name_ << " [options]" << endl
    << "Options are:" << endl
    << "  -optimize_peaks   switch for the optimization of peak parameters. Valid options are 'on' or 'off' (default: 'off')" << endl
    << "  -in <file>        input mzData file name" << endl
    << "  -out <file>       output mzData file name" << endl
    << endl;
  }

  void printToolHelpOpt_()
  {
    cerr << endl
    << tool_name_ << endl
    << endl
    << "INI options:" << endl
    << "  optimize_peaks   switch for the optimization of peak parameters. Valid options are 'on' or 'off' (default: 'off')" << endl
    << "  in               input mzData file name" << endl
    << "  out              output mzData file name" << endl
    << endl
    << "INI File example section:" << endl
    << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"optimize_peaks\" value=\"on\" type=\"string\"/>" << endl;
  }

  void setOptionsAndFlags_()
  {
    options_["-out"] = "out";
    options_["-in"] = "in";
    options_["-optimize_peaks"] = "optimize_peaks";
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

    //optimze flag
    String optimize_peaks = getParamAsString_("optimize_peaks");
    writeDebug_(String("Optimization of peaks: ") + optimize_peaks, 1);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzDataFile mz_data_file;
    MSExperiment<DRawDataPoint<1> > ms_exp_raw;
    mz_data_file.load(in,ms_exp_raw);


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String ini_location = String(tool_name_) + ":" + String(instance_number_) + ":";
    Param pepi_param = getParamCopy_(ini_location);

    if (optimize_peaks != "")
    {
      //optimization
      if (optimize_peaks == "on")
      {
        pepi_param.setValue("Optimization:SkipOptimization","no");
      }
      else if (optimize_peaks == "off")
      {
        pepi_param.setValue("Optimization:SkipOptimization","yes");
      }
      else
      {
        writeLog_(String("Invalid option '") + optimize_peaks + "' given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      pepi_param.setValue("Optimization:SkipOptimization","yes");
    }

    std::cout << pepi_param << std::endl;

    PeakPickerCWT peak_picker(pepi_param);
    
    MSExperiment<DPickedPeak<1> > ms_exp_peaks;
    // copy the experimental settings
    static_cast<ExperimentalSettings&>(ms_exp_peaks) = ms_exp_raw;    

		// pick peaks on each scan
    for (unsigned int i = 0; i < ms_exp_raw.size(); ++i)
    {
    	MSSpectrum<DPickedPeak<1> > spectrum;
    	 
    	peak_picker.pick(ms_exp_raw[i],spectrum);
    	
    	// if any peaks are found copy the spectrum settings
    	if (spectrum.size() > 0)
    		{
    			// copy the spectrum settings
    			static_cast<SpectrumSettings&>(spectrum) = ms_exp_raw[i];  
    			spectrum.setType(SpectrumSettings::PEAKS);
    				
    			// copy the spectrum information
    			spectrum.getPrecursorPeak() = ms_exp_raw[i].getPrecursorPeak();
					spectrum.setRetentionTime(ms_exp_raw[i].getRetentionTime());
					spectrum.setMSLevel(ms_exp_raw[i].getMSLevel());
					spectrum.getName() = ms_exp_raw[i].getName();
				
    			ms_exp_peaks.push_back(spectrum);
    	}
    }
   

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    mz_data_file.store(out,ms_exp_peaks);

    return OK;
  }
};

/// @endcond


int main( int argc, char ** argv )
{
  TOPPPeakPicker tool;
  return tool.main(argc,argv);
}

