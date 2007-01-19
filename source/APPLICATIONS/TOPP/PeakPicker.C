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
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
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
	In case of low-resoluted data an optional step for the separation of
	overlapping peaks can be added.
	In an optional third stage, the resulting fit can be further improved using
	techniques from nonlinear optimization.
	
	<b>Finding the right parameters:</b><br>
	Finding the right parameters for the peak picking is not trivial. The default parameters
	will not work on most datasets. In order to find the standard parameters we propose the
	following procedure:
	-# Extract a single scan from your data
		- You can do that using the \ref DTAExtractor tool
		- You should extract a scan from the middle of the gradient. To get information about the
	    gradient length, you can use the \ref FileInfo tool.
	-# Load the scan in TOPPView
	-# Experiment with the peak picking parameters until you have found the proper settings
		- You can find the peak picker in \ref TOPPView in the menu 'Tools'
		- The most important parameters are <tt>peak height</tt> and <tt>full width at half max</tt>. 
			If you set them too high, not many peaks will be found!
	-# Set the parameters <tt>peak_bound</tt>, <tt>peak_bound_ms2_level</tt>, <tt>fwhm_bound</tt> 
		and <tt>scale</tt> in an INI file.
		- <tt>scale</tt> can be set to the equal to the <tt>peak_bound</tt> parameter.
	
	@todo INI file for NoiseFilter, PeakPicker and MapAlignment
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPicker
      : public TOPPBase
{
 public:
  TOPPPeakPicker()
		: TOPPBase("PeakPicker","find mass spectrometric peaks in LC/MS raw data")
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
  	registerStringOption_("in","<file>","","input mzData file (raw data)");
		registerStringOption_("out","<file>","","output mzData file (peak data)");
		addEmptyLine_();
  	addText_("Parameters for the peak picker algorithm can be given in the 'algorithm' part of INI file:\n"
							"<NODE algorithm>\n"
						 	" <NODE name=\"SignalToNoiseEstimationParameter\">\n "
						  "   <ITEM name=\"WindowLength\" value=\"200\" type=\"int\" />\n "
				      "   <ITEM name=\"BinCount\" value=\"30\" type=\"int\" />\n "
        			"   <ITEM name=\"MinReqElementsInWindow\" value=\"10\" type=\"int\" />\n"
						  "   <ITEM name=\"NoiseEmptyWindow\" value=\"2\" type=\"int\" /> \n"
						  " </NODE> \n"
						  " <NODE name=\"Optimization\"> \n"
							"  <ITEM name=\"optimization\" value=\"no\" type=\"string\"/>\n"
						  " </NODE> "
							"  <NODE name=\"wavelet_transform\">\n"
							"    <ITEM name=\"scale\" value=\"0.2\" type=\"float\" />\n"
							"  </NODE>\n"
							"  <NODE name=\"thresholds\">\n"
							"    <ITEM name=\"peak_bound\" value=\"100\" type=\"float\" />\n"
							"    <ITEM name=\"peak_bound_ms2_level\" value=\"30\" type=\"float\" />\n"
							"    <ITEM name=\"fwhm_bound\" value=\"0.1\" type=\"float\"/>\n"
							"  </NODE>\n"
						  " <NODE name=\"deconvolution\"> \n"
						  "  <ITEM name=\"skip_deconvolution\" value=\"yes\" type=\"string\"/>\n"
						  " </NODE>\n"
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
    //-------------------------------------------------------------
    // Init peak picker
    //-------------------------------------------------------------
		Param pepi_param = getParam_().copy("algorithm:",true);

		
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
