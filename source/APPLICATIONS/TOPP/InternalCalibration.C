// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_InternalCalibration InternalCalibration
	
	@brief Performs an internal calibration on an MS experiment.
	
	This is a simple calibration method: given a list of reference masses and an MS experiment, 
	the relative errors of the peaks in the data are approximated by linear interpolation and
	subtracted from the data. This is done scanwise, i.e. at least two reference masses need to
	be present in each scan, otherwise the scan can't be calibrated. If the input file contains
	raw data, an additional peak picking step is performed.
	
	@note The default input is raw data, if you have peak data, please use the flag peak_data.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_InternalCalibration.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPInternalCalibration
      : public TOPPBase
 {
 public:
  TOPPInternalCalibration()
    : TOPPBase("InternalCalibration","Applies an internal calibration.")
  {
  }

 protected:

	 void registerOptionsAndFlags_()
	 {
		 registerInputFile_("in","<file>","","input raw data or peak file ");
		 setValidFormats_("in",StringList::create("mzML"));
		 registerOutputFile_("out","<file>","","output file ");
	   setValidFormats_("out",StringList::create("mzML"));
		 registerInputFile_("ref_masses","<file>","","input file containing reference masses (one per line)",true);
		 registerFlag_("peak_data","set this flag, if you have peak data, not raw data");
		 addEmptyLine_();
		 addText_("If you want to calibrate raw data, it is necessary to perform a peak picking step before the "
							"actual calibration is done. The parameters for the peak picking step can be given "
							"given in the 'algorithm' part of INI file in the subsection PeakPicker.");
		 addEmptyLine_();
		 registerSubsection_("algorithm","Settings for the peak picking step.");
	 }
	 
	 Param getSubsectionDefaults_(const String& /* section*/) const
	 {
		 Param tmp;
		 tmp.insert("PeakPicker:",PeakPickerCWT().getDefaults());
		 return tmp;
	 }
	 
	 
	 ExitCodes main_(int , const char**)
	 {
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		String ref = getStringOption_("ref_masses");
		//-------------------------------------------------------------
		// init InternalCalibration
		//-------------------------------------------------------------
		
		InternalCalibration calib;
		Param param = getParam_().copy("algorithm:",true);
		calib.setParameters(param);
		
		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		MSExperiment<Peak1D > ms_exp_raw;
		
		MzMLFile mz_data_file;
		mz_data_file.setLogType(log_type_);
		mz_data_file.load(in,ms_exp_raw);
		
		
		
		vector<double> ref_masses;
		TextFile ref_file;
		
		
		ref_file.load(ref,true);
		
		for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
		{
			ref_masses.push_back(atof(iter->c_str()));
		}
		
		//-------------------------------------------------------------
		// perform calibration
		//-------------------------------------------------------------
		
		calib.calibrate(ms_exp_raw,ref_masses,getFlag_("peak_data"));
		
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		//annotate output with data processing info
		addDataProcessing_(ms_exp_raw, getProcessingInfo_(DataProcessing::CALIBRATION));
		
		mz_data_file.store(out,ms_exp_raw);
		
		return EXECUTION_OK;
  }

};


int main( int argc, const char** argv )
{
  TOPPInternalCalibration tool;
  return tool.main(argc,argv);
}

/// @endcond
