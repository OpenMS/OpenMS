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
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_PeakPicker PeakPicker
	
	@brief A tool for peak detection in profile data
	
	Executes the peak picking algorithm as described by Lange et al. (2006) Proc. PSB-06.
	
	The conversion of the ''raw'' ion count data acquired
	by the machine into peak lists for further processing
	is usually called peak picking. Our algorithm is independent
	of the underlying machine or ionization method, and is able
	to resolve highly convoluted and asymmetric signals.
	The method uses the multiscale nature of spectrometric data by
	first detecting the mass peaks in the wavelet-transformed signal
	before a given asymmetric peak function is fitted to the profile data.
	In case of low-resoluted data, an optional step for the separation of
	overlapping peaks can be added.
	In an optional third stage, the resulting fit can be further improved using
	techniques from nonlinear optimization.
	
	How to find @ref TOPP_example_signalprocessing_parameters is explained in the TOPP tutorial. 
	
	In the following table you, can find example values of the most important parameters for 
	different instrument types. @n These parameters are not valid for all instruments of that type,
	but can be used as a starting point for finding suitable parameters.
	<table>
		<tr>
			<td>&nbsp;</td>
			<td><b>Q-TOF</b></td>
			<td><b>LTQ Orbitrap</b></td>
		</tr>
		<tr>
			<td><b>signal_to_noise</b></td>
			<td>2</td>
			<td>0</td>
		</tr>
		<tr>
		<td><b>peak_width</b></td>
			<td>0.1</td>
			<td>0.012</td>
		</tr>
	</table>
	
	In order to impove the results of the peak detection on low resolution data @ref TOPP_NoiseFilter and @ref TOPP_BaselineFilter can be applied.
	For high resolution data this is not necessary.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PeakPicker.cli
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
  	registerInputFile_("in","<file>","","input profile data file ");
		setValidFormats_("in",StringList::create("mzML"));
		registerOutputFile_("out","<file>","","output peak file ");
	  setValidFormats_("out",StringList::create("mzML"));
		registerStringOption_("type","<name>","","peak detection algorithm type",true);
		setValidStrings_("type", getToolList()[toolName_()] );
		addEmptyLine_();
  	addText_("Parameters for the peak picker algorithm can be given in the 'algorithm' part of INI file.");
  	registerSubsection_("algorithm","Algorithm parameters section");
  }
  
	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		Param tmp;
		
		if (type == "wavelet")
    {
      tmp = PeakPickerCWT().getDefaults();
    }
    else if (type == "high_res")
    {
      tmp = PeakPickerHiRes().getDefaults();
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
    MSExperiment<Peak1D > ms_exp_raw;
    mz_data_file.load(in,ms_exp_raw);
		
		//check for peak type (profile data required)
		if (PeakTypeEstimator().estimateType(ms_exp_raw[0].begin(),ms_exp_raw[0].end())==SpectrumSettings::PEAKS)
		{
			writeLog_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
		}

		//check if spectra are sorted
		for (Size i=0; i<ms_exp_raw.size(); ++i)
		{
			if (!ms_exp_raw[i].isSorted())
			{
				writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
				return INCOMPATIBLE_INPUT_DATA;
			}
		}

    //-------------------------------------------------------------
    // pick
    //-------------------------------------------------------------
    MSExperiment<> ms_exp_peaks;
    
		Param pepi_param = getParam_().copy("algorithm:",true);		
		writeDebug_("Parameters passed to PeakPicker", pepi_param,3);
		
    if (type == "wavelet")
    {	
    	PeakPickerCWT pp;
      pp.setLogType(log_type_);
			pp.setParameters(pepi_param);
			pp.pickExperiment(ms_exp_raw,ms_exp_peaks);
    }
    else if (type == "high_res")
    {	
    	PeakPickerHiRes pp;
      pp.setLogType(log_type_);
			pp.setParameters(pepi_param);
			pp.pickExperiment(ms_exp_raw,ms_exp_peaks);
    }

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		//annotate output with data processing info
		addDataProcessing_(ms_exp_peaks, getProcessingInfo_(DataProcessing::PEAK_PICKING));

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
