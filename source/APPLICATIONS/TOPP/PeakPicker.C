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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/APPLICATIONS/TOPPBase2.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

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
      : public TOPPBase2
{
 public:
  TOPPPeakPicker()
		: TOPPBase2("PeakPicker","find mass spectrometric peaks in LC/MS raw data")
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
  	registerStringOption_("in","<file>","","input mzData file (raw data)");
		registerStringOption_("out","<file>","","output mzData file (peak data)");
    registerFlag_("optimize_peaks","flag that turns on the optimization of peak parameters");
		addEmptyLine_();
  	addText_("Parameters for the peak picker algorithm can be given in the 'algorithm' part of INI file:\n"
							"<NODE algorithm>\n"
							"  <NODE name=\"signal_to_noise_estimation\">\n"
							"    <ITEM name=\"bucket\" value=\"10\" type=\"int\" />\n"
							"    <ITEM name=\"window\" value=\"700\" type=\"int\" />\n"
							"  </NODE>\n"
							"  <ITEM name=\"optimize_peaks\" value=\"off\" type=\"string\"/>\n"
							"  <NODE name=\"wavelet_transform\">\n"
							"    <ITEM name=\"scale\" value=\"0.2\" type=\"float\" />\n"
							"  </NODE>\n"
							"  <NODE name=\"thresholds\">\n"
							"    <ITEM name=\"peak_bound\" value=\"100\" type=\"float\" />\n"
							"    <ITEM name=\"peak_bound_ms2_level\" value=\"30\" type=\"float\" />\n"
							"    <ITEM name=\"fwhm_bound\" value=\"0.1\" type=\"float\"/>\n"
							"  </NODE>\n"
							"</NODE>");
		addEmptyLine_();
  	addText_("This application implements an algorithm for peak picking as\n"
				     "described in Lange et al. (2006) Proc. PSB-06. ");
  	registerSubsection_("algorithm");
  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    bool optimize_peaks = getFlag_("optimize_peaks");
		
    //-------------------------------------------------------------
    // Init peak picker
    //-------------------------------------------------------------
		Param pepi_param = getParam_().copy("algorithm:",true);
		
    //optimization
    if (optimize_peaks)
    {
      pepi_param.setValue("Optimization:skip_optimization","no");
    }
    else
    {
      pepi_param.setValue("Optimization:skip_optimization","yes");
    }
		
		writeDebug_("Parameters passed to PeakPickerCWT", pepi_param,3);
    PeakPickerCWT peak_picker;
    peak_picker.setParam(pepi_param);

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
    // pick
    //-------------------------------------------------------------

    MSExperiment<DPickedPeak<1> > ms_exp_peaks;
    peak_picker.pickExperiment(ms_exp_raw,ms_exp_peaks);
  
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		ms_exp_peaks.getProcessingMethod().setSpectrumType(SpectrumSettings::PEAKS);
		mz_data_file.store(out,ms_exp_peaks);

		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
  TOPPPeakPicker tool;
  return tool.main(argc,argv);
}

/// @endcond
