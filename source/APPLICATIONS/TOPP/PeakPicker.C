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
	
	@brief A tool for peak detection in raw data
	
	Executes the peak picking algorithm as described by Lange et al. (2006) Proc. PSB-06.
	
	The conversion of the ''raw'' ion count data acquired
	by the machine into peak lists for further processing
	is usually called peak picking. Our algorithm is independent
	of the underlying machine or ionization method, and is able
	to resolve highly convoluted and asymmetric signals.
	The method uses the multiscale nature of spectrometric data by
	first detecting the mass peaks in the wavelet-transformed signal
	before a given asymmetric peak function is fitted to the raw data.
	In case of low-resoluted data, an optional step for the separation of
	overlapping peaks can be added.
	In an optional third stage, the resulting fit can be further improved using
	techniques from nonlinear optimization.
	
	How to find @ref TOPP_example_signalprocessing_parameters is explained in the TOPP tutorial. 

	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPicker
      : public TOPPBase
{
 public:
  TOPPPeakPicker()
		: TOPPBase("PeakPicker","Finds mass spectrometric peaks in profile mass spectra.")
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
  	registerInputFile_("in","<file>","","input raw data file ");
		setValidFormats_("in",StringList::create("mzData"));
		registerOutputFile_("out","<file>","","output peak file ");
	  setValidFormats_("out",StringList::create("mzData"));
		addEmptyLine_();
  	addText_("Parameters for the peak picker algorithm can be given in the 'algorithm' part of INI file.");
		addEmptyLine_();
  	addText_("This application implements an algorithm for peak picking as\n"
				     "described in Lange et al. (2006) Proc. PSB-06. ");
  	registerSubsection_("algorithm","Algorithm parameters section");
  }
  
  Param getSubsectionDefaults_(const String& /* section*/) const
  {
      return PeakPickerCWT().getDefaults();
  }

  ExitCodes main_(int , const char**)
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
    peak_picker.setLogType(log_type_);
		peak_picker.setParameters(pepi_param);
		
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
    // pick
    //-------------------------------------------------------------

    MSExperiment<> ms_exp_peaks;
    peak_picker.pickExperiment(ms_exp_raw,ms_exp_peaks);
  
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		ms_exp_peaks.getProcessingMethod().setSpectrumType(SpectrumSettings::PEAKS);
		mz_data_file.store(out,ms_exp_peaks);


		
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPPeakPicker tool;
  return tool.main(argc,argv);
}

/// @endcond
