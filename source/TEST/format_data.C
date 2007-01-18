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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/ANDIFile.h>
#include <OpenMS/FORMAT/HANDLERS/ANDIHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>


// This program created the MzData, MzXML and ANDIFile data used in format tests.
// It is capable of creating data of variable size to test visualization.

using namespace OpenMS;
using namespace std;

int main(int argc, const char* argv[])
{
	const String NAME = (argc>2)? argv[2] : "tmp";
	const Size SPEC_NUM = (argc>1)? QString(argv[1]).toInt() : 100;

	cerr << "\nBuilding test data for MzXML, MzData and ANDIFile with " << SPEC_NUM << " scans.\n"
			 <<   "------------------------------------------------------------------------------------\n";

	vector< vector<float> > mz, intens;
	// create 2D-map looking like a half-pyramid
	for (UnsignedInt j=0, k=1; j<SPEC_NUM; j++, k+=2)
	{
		vector<float> tmp_mz, tmp_int;
		for (int i=0, m=-j; i<(int)k; i++, m++)
		{
			tmp_mz.push_back(SPEC_NUM*40-j*10+10*i);		
			tmp_int.push_back((j+1)*100-abs(m)*100);
		}
		mz.push_back(tmp_mz);
		intens.push_back(tmp_int);
	}

	//---------------------------------------------------------------------------
	// Create MzData test file
	//---------------------------------------------------------------------------
	
	String testfilename="data/" + NAME + ".mzData";
	cerr << "Creating file: " << testfilename << endl;

	Base64 b64;
	
	ofstream mzdata(testfilename.c_str());

	mzdata << "<!-- -*- Mode: XML; tab-width: 2; -*- -->\n<mzData>\n"
					<< "\t<spectrumList count=\"" << SPEC_NUM << "\">\n";

	for (UnsignedInt spec=0; spec<SPEC_NUM; spec++)
	{
		mzdata << "\t\t<spectrum id=\"" << (spec+1) << "\">\n"
						<< "\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n"
						<< "\t\t\t\t<spectrumInstrument MSLevel=\"1\" mzRangeStart=\"300.0\" mzRangeStop=\"1500.0\">\n"
						<< "\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000038\" name=\"TimeInMinutes\" value=\""
						<< (spec+1) << "\"/>\n"
						<< "\t\t\t\t</spectrumInstrument>\n\t\t\t</spectrumSettings>\n\t\t</spectrumDesc>\n";

		float* tmp = b64.getFloatBuffer(mz[spec].size());
		
		for (UnsignedInt i=0; i<mz[spec].size(); i++)
		{
			tmp[i] = mz[spec][i];
		}
		mzdata << "\t\t\t<mzArrayBinary>\n\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
						<< mz[spec].size()
						<< "\">"
						<< b64.encodeFloat()
						<< "</data>\n\t\t\t</mzArrayBinary>\n";

		tmp = b64.getFloatBuffer(intens[spec].size());
		for (UnsignedInt i=0; i<intens[spec].size(); i++)
		{
			tmp[i] = intens[spec][i];
		}
		mzdata << "\t\t\t<intenArrayBinary>\n\t\t\t\t<data precision=\"32\" endian=\"little\" length=\""
						<< intens[spec].size()
						<< "\">"
						<< b64.encodeFloat()
						<< "</data>\n\t\t\t</intenArrayBinary>\n\t\t</spectrum>\n";
	}
	mzdata << "\t</spectrumList>\n</mzData>\n";

	//---------------------------------------------------------------------------
	// Create MzXML test file
	//---------------------------------------------------------------------------
	

	testfilename="data/" + NAME + ".mzXML";
	cerr << "Creating file: " << testfilename << endl;

	ofstream mzxml(testfilename.c_str());


	mzxml << "<!-- -*- Mode: XML; tab-width: 2; -*- -->\n<msRun scanCount=\"" 
					<< SPEC_NUM << "\" startTime=\"PT0.220000S\" endTime=\"PT3180.090000S\">\n";

	for (UnsignedInt spec=0; spec<SPEC_NUM; spec++)
	{
		mzxml << "\t<scan num=\"" << (spec+1) << "\"MSLevel=\"1\" peaksCount=\"" << mz[spec].size()
						<< "\" polarity=\"+\" scanType=\"full\" centroided=\"1\" retentionTime=\"PT"
						<< 60*(spec+1) << "S\">\n\t\t<peaks precision=\"32\">";
			
		float* tmp = b64.getFloatBuffer(mz[spec].size()*2);
		
		for (UnsignedInt i=0; i<mz[spec].size(); i++)
		{
			tmp[2*i]   = mz[spec][i];
			tmp[2*i+1] = intens[spec][i];
		}
		mzxml << b64.encodeFloatCorrected() << "</peaks>\n\t</scan>\n";
	}

	mzxml << "</msRun>\n";


	//---------------------------------------------------------------------------
	// Create ANDI test file
	//---------------------------------------------------------------------------
 
	testfilename="data/" + NAME + ".cdf";
	cerr << "Creating file: " << testfilename << endl << endl;

	unsigned long nscans = SPEC_NUM;
	unsigned long ninst = 3;
	MS_Sample_Data fill_sd = {"30", "31","32","33","34","35","36","37","38","39","40","41","42",state_solid};
	MS_Test_Data fill_td = {	separation_glc,inlet_membrane,2.7,ionization_ei,
														polarity_plus,23.56,56.23,"43",12.3,"44","45",
														1.2,2.3,3.4,4.5,detector_em,5.6,6.7,resolution_constant,
														"46",function_scan,direction_up,law_linear,12.2,"47","48","49","50"};
	MS_Raw_Data_Global fill_rdg = {	nscans,true,false,1.0,1.0,1.0,0.0,mass_m_z,time_seconds,
																	intensity_counts,intensity_volts,data_float,
																	data_float,data_float,"51","52","53",0,0,0,
																	0,0,0,0,0,0,0,0,0,"54"};
	MS_Instrument_Data fill_id[3] = {{0, "i1", "i2","i3","i4","i5","i6","i7","i8","i9","i10"},
																	 {1, "i11", "i12","i13","i14","i15","i16","i17","i18","i19","i20"},
																	 {2, "i21", "i22","i23","i24","i25","i26","i27","i28","i29","i30"}};
	// First 4 values and the value after "28" are set by NETCDF regardless of values set here
	MS_Admin_Data fill_ad = {"C1+C2", "1.0.1","2.3.2","English","5","6","7","20031211093000+0001","9",
													 "20021201091000+0002","11","12","13","14","15","16","17","18","19",
													 "20011201081000+0003","21","22","23","24","25","26","27","28","",expt_centroid ,123,456,ninst};

	MS_Admin_Data admin_data;
  MS_Sample_Data sample_data;
  MS_Test_Data test_data;
  MS_Raw_Data_Global raw_global_data; 
  MS_Instrument_Data inst_data;
  MS_Raw_Library lib_data;
	MS_Raw_Per_Scan raw_data;

	unsigned long index;
	int file_id;
	int err_code;
  ncopts = 0; //NC_VERBOSE;
	
	char* name = (char*)testfilename.c_str();
	file_id = ms_open_write(name ,expt_centroid,nscans,ninst,
													 data_float,data_float,data_long,TRUE,FALSE);
	if (MS_ERROR==file_id) return 1;
	ms_init_global( 0, &admin_data, &sample_data, &test_data, &raw_global_data);

	admin_data = fill_ad;
	sample_data = fill_sd;
	test_data = fill_td;
	raw_global_data = fill_rdg;

	err_code =  ms_write_global( file_id, &admin_data, &sample_data, &test_data, &raw_global_data);
	if (MS_ERROR==err_code) return 1;

	ms_init_instrument(0,&inst_data);
	for ( index = 0; index < ninst; index++)
  {
		inst_data = fill_id[index];
		err_code = ms_write_instrument( file_id, &inst_data);
		if (MS_ERROR==err_code) return 1;
		ms_init_instrument(0,&inst_data);
	}

	ms_init_per_scan(0, &raw_data, &lib_data);
	for (index = 0; index < nscans; index++) {
		const UnsignedInt SIZE = mz[index].size();
		float* ms_m = (float*) new float[SIZE];
		float* ms_i = (float*) new float[SIZE];

		for (UnsignedInt i=0; i<SIZE; i++)
		{
			ms_m[i] = mz[index][i];
			ms_i[i] = intens[index][i];
		}

		MS_Raw_Per_Scan fill_rps = {index,SIZE,0,index,0,0,0,(index+1)*60,0,0,0,0,0,0,0,ms_m, NULL, ms_i,	NULL, NULL};

		raw_data = fill_rps;
		err_code = ms_write_per_scan(file_id, &raw_data, NULL);
		if (MS_ERROR==err_code) return 1;

		//ms_init_per_scan( 0, &raw_data, &lib_data);
		delete [] ms_m;
		delete [] ms_i;
	}

	ms_init_global( 0, &admin_data, &sample_data, &test_data, &raw_global_data);
	ms_close( file_id);
}
