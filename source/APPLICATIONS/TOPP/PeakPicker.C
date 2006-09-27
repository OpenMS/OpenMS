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
 
   @brief Executes the peak picking algorithm as described by Lange et al. (2006) Proc. PSB-06.
 
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
 	 
 	 @todo Add the three main parameters to the the command line / ini file (Eva)
 	 
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
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
				 << endl
				 << "Options are:" << endl
				 << "  -optimize_peaks   flag that turns on for the optimization of peak parameters" << endl
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
				 << "  optimize_peaks   flag that turns on for the optimization of peak parameters" << endl
				 << "  in <file>        input mzData file name" << endl
				 << "  out <file>       output mzData file name" << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"optimize_peaks\" value=\"\" type=\"string\"/>" << endl;
  }

  void setOptionsAndFlags_()
  {
    options_["-out"] = "out";
    options_["-in"] = "in";
    flags_["-optimize_peaks"] = "optimize_peaks";
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
    bool optimize_peaks = getParamAsBool_("optimize_peaks");
    if (optimize_peaks)
    {
    	writeDebug_(String("Optimization of peaks: ON"), 1);
		}
		else
		{
			writeDebug_(String("Optimization of peaks: OFF"), 1);
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
    String ini_location = String(tool_name_) + ":" + String(instance_number_) + ":";
    Param pepi_param = getParamCopy_(ini_location);

    //optimization
    if (optimize_peaks)
    {
      pepi_param.setValue("Optimization:SkipOptimization","no");
    }
    else
    {
      pepi_param.setValue("Optimization:SkipOptimization","yes");
    }

    PeakPickerCWT peak_picker(pepi_param);

    MSExperiment<DPickedPeak<1> > ms_exp_peaks;
    peak_picker.pickExperiment(ms_exp_raw,ms_exp_peaks);
    
  
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------


		mz_data_file.store(out,ms_exp_peaks);

		return OK;
	}
};


int main( int argc, char ** argv )
{
  TOPPPeakPicker tool;
  return tool.main(argc,argv);
}

/// @endcond
