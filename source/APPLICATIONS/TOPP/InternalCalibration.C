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
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>

using namespace OpenMS;
using namespace std;



//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**	
   @brief Performs an internal calibration.

	 This is a simle calibration method: given a list of reference masses, 
	 the relative errors of the peaks in the data are approximated by linear interpolation and
	 subtracted from the data. This is done scanwise, i.e. at least two reference masses need to
	 be present in each scan, otherwise the scan can't be calibrated.

	
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPInternalCalibration
      : public TOPPBase
 {
 public:
  TOPPInternalCalibration()
    : TOPPBase("InternalCalibration","Apply internal calibration.")
  {
  }

 protected:

  void registerOptionsAndFlags_()
  {
    registerStringOption_("in","<input file>","","input file (mzData or mzXML)");
    registerStringOption_("out","<output file>","","output file (mzData or mzXML)");
    registerStringOption_("ref_masses","<reference file>","","file containing reference masses(one per line)",true);
		registerDoubleOption_("peak_bound","<peak bound>",1000,"minimal intensity for the peak picking step",false);
  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String ref = getStringOption_("ref_masses");
		double peak_bound = getDoubleOption_("peak_bound");
    //-------------------------------------------------------------
    // init InternalCalibration
    //-------------------------------------------------------------

    InternalCalibration calib;
    calib.setPeakBound(peak_bound);
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MSExperiment<RawDataPoint1D > ms_exp_raw;
    
    if(in.hasSuffix("mzXML"))
      {
				MzXMLFile mz_xml_file;
				mz_xml_file.load(in,ms_exp_raw);
				
      }
    else if(in.hasSuffix("mzData"))
      {
				MzDataFile mz_data_file;
				mz_data_file.load(in,ms_exp_raw);
				
      }
    else
      {
				std::cout << "Unknown input format. Aborting!"<<std::endl;
				return INPUT_FILE_NOT_READABLE;
      }
    
    float mz;
    vector<double> ref_masses;
    FILE *stream;		
    if((stream = fopen( ref.c_str(), "rt"))!= NULL)
      {
				// read ref masses
				while (fscanf(stream,"%f\n",&mz)!=EOF)
					{
						ref_masses.push_back(mz);
					}
      }
    else
      {
				std::cout << "Reference file not found. Aborting!"<<std::endl;
				return INPUT_FILE_NOT_FOUND;
      }
		
		
    //-------------------------------------------------------------
    // perform calibration
    //-------------------------------------------------------------
		
    calib.calibrate(ms_exp_raw,ref_masses);

    
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if(out.hasSuffix("mzXML"))
      {
				MzXMLFile mz_xml_file;
				mz_xml_file.store(out,ms_exp_raw);
				
      }
    else if(out.hasSuffix("mzData"))
      {
				MzDataFile mz_data_file;
				mz_data_file.store(out,ms_exp_raw);
				
      }
    else
      {
				std::cout << "Unknown output format. Aborting!"<<std::endl;
				return CANNOT_WRITE_OUTPUT_FILE;
      }
    
    return EXECUTION_OK;
  }
};


int main( int argc, char ** argv )
{
  TOPPInternalCalibration tool;
  return tool.main(argc,argv);
}

/// @endcond
