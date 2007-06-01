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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/CALIBRATION/ExternalCalibration.h>

using namespace OpenMS;
using namespace std;



//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	 @page ExternalCalibration ExternalCalibration
	 
   @brief Performs an external calibration for tof spectra.

	 Given one or more calibrant spectra containing flight times, the instrument's calibration constants and the
	 expected masses the quadratic function y_i = a + b*x_i + c*x_i^2 is fitted, where x_i is the ith flight time.
	 If there are more than one calibrant spectra the coefficients a, b and c are averaged. The fitted function is
	 then used to convert the flight times of the given experiment to m/z-values.

	 You can choose to calibrate picked or raw data. If you use picked data, set the flag peak_data. If you have
	 raw data an additional peak picking step for the calibrant spectra is needed, the parameters for the
	 peak picker can be set in the ini-file.

	
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPExternalCalibration
      : public TOPPBase
 {
 public:
  TOPPExternalCalibration()
    : TOPPBase("ExternalCalibration","Apply external calibration.")
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
    registerStringOption_("in","<input file>","","input mzData file (peak or raw data)");
    registerStringOption_("out","<output file>","","output mzData file (peak or raw data)");
		registerStringOption_("ext_calibrants","<input file>","","mzData file containing the external calibrant spectra (peak or raw data)");
    registerStringOption_("ref_masses","<reference file>","","file containing reference masses of the external calibrant spectra (one per line)",true);
		registerStringOption_("tof_const","<file>","","File containing TOF conversion constants."
													" These can be either two or three constants\n" 
													"per set, depending on the conversion type. Either one set for all calibrant spectra \n"
													"(tab separated), or one for each spectrum.\n"
													"For a detailed description, please have a look at the doxygen documentation."
													"(one set, tab separated, per line)",true);
		registerFlag_("peak_data","set this flag, if you have peak data, not raw data");
		addText_("\nIf you want to calibrate raw data, it is necessary to perform a peak picking step before the "
							"actual calibration is done. \nThe parameters for the peak picking step can be given "
							"given in the 'algorithm' part of INI file in the subsection PeakPicker, e.g.:\n"
							"<NODE name=\"algorithm\">\n"
						  " <NODE name=\"PeakPicker\">\n"
							"  <NODE name=\"wavelet_transform\">\n"
							"    <ITEM name=\"scale\" value=\"0.2\" type=\"float\" />\n"
							"  </NODE>\n"
							"  <NODE name=\"thresholds\">\n"
							"    <ITEM name=\"peak_bound\" value=\"100\" type=\"float\" />\n"
							"    <ITEM name=\"correlation\" value=\"0.5\" type=\"float\" />\n"
							"    <ITEM name=\"fwhm_bound\" value=\"0.1\" type=\"float\"/>\n"
							"  </NODE>\n"
						  " </NODE>\n"
						  "</NODE>");
		 addEmptyLine_();
		 registerSubsection_("algorithm","Algorithm section for peak picking");
  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
		String in_calib = getStringOption_("ext_calibrants");
    String ref = getStringOption_("ref_masses");
		String conv = getStringOption_("tof_const");
		bool peak_data = getFlag_("peak_data");
    //-------------------------------------------------------------
    // init ExternalCalibration
    //-------------------------------------------------------------

    ExternalCalibration calib;
		calib.setLogType(log_type_);
		Param param = getParam_().copy("algorithm:",true);
		calib.setParameters(param);
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MSExperiment<RawDataPoint1D > ms_exp_calib,ms_exp_raw;
		MSExperiment<PickedPeak1D > ms_exp_p,ms_exp_calib_p;
		MzDataFile mz_data_file;
		mz_data_file.setLogType(log_type_);
		if(peak_data)
			{
				mz_data_file.load(in_calib,ms_exp_calib_p);
				mz_data_file.load(in,ms_exp_p);
			}
		else
			{
				mz_data_file.load(in_calib,ms_exp_calib);
				mz_data_file.load(in,ms_exp_raw);
			}
		vector<double> ref_masses;
		TextFile ref_file;
		ref_file.load(ref,true);
		
		for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
			{
				ref_masses.push_back(atof(iter->c_str()));
			}
		TextFile const_file;
		const_file.load(conv,true);
		std::vector<String> vec;
		TextFile::Iterator iter = const_file.begin();
		iter->split('\t',vec);

		std::vector<double> ml1,ml2,ml3;
		ml1.push_back(atof(vec[0].c_str()));
		ml2.push_back(atof(vec[1].c_str()));
		if(vec.size()==3)
			{
				ml3.push_back(atof(vec[2].c_str()));				
			}
		++iter;
		
		for(; iter != const_file.end(); ++iter)
			{
				iter->split('\t',vec);
				ml1.push_back(atof(vec[0].c_str()));
				ml2.push_back(atof(vec[1].c_str()));
				if(vec.size()==3)
					{
						ml3.push_back(atof(vec[2].c_str()));				
					}
			}

		if(ml1.size() != 1 && (!(peak_data && ml1.size() == ms_exp_calib_p.size()) && ml1.size() != ms_exp_calib.size()))
			{
				writeLog_("Incorrect number of calibration constants given. Aborting!");
				return INPUT_FILE_CORRUPT;
			}
		calib.setML1s(ml1);
		calib.setML2s(ml2);
		if(!ml3.empty()) calib.setML3s(ml3);
    //-------------------------------------------------------------
    // perform calibration
    //-------------------------------------------------------------
		if(peak_data)
			{
				calib.calibrate(ms_exp_calib_p,ms_exp_p,ref_masses);
			}
		else calib.calibrate(ms_exp_calib,ms_exp_raw,ref_masses);
    
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
		if(peak_data)
			{
				mz_data_file.store(out,ms_exp_p);
			}
		else mz_data_file.store(out,ms_exp_raw);

    return EXECUTION_OK;
  }
};

Param TOPPBase::getSubsectionDefaults_(const String& /*section*/) const
{
  // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the PeakPicker param
  Param tmp;
  tmp.insert("PeakPicker:",PeakPickerCWT().getDefaults());
  return tmp;
}

int main( int argc, char ** argv )
{
  TOPPExternalCalibration tool;
  return tool.main(argc,argv);
}

/// @endcond
