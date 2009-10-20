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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/TextFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TOFCalibration, "$Id$")

/////////////////////////////////////////////////////////////

TOFCalibration* ptr = 0;
START_SECTION((TOFCalibration()))
  ptr = new TOFCalibration;
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~TOFCalibration()))
  delete ptr;
END_SECTION

TOFCalibration tc;

START_SECTION((const std::vector<double>& getML1s() const))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML1s(vec);
  TEST_EQUAL(tc.getML1s()== vec,true)
END_SECTION

START_SECTION((const std::vector<double>& getML2s() const))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML2s(vec);
  TEST_EQUAL(tc.getML2s()== vec,true)
END_SECTION

START_SECTION((const std::vector<double>& getML3s() const))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML3s(vec);
  TEST_EQUAL(tc.getML3s()== vec,true)
END_SECTION  

  
START_SECTION((void setML1s(const std::vector< double > &ml1s)))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML1s(vec);
  TEST_EQUAL(tc.getML1s()== vec,true)
END_SECTION

START_SECTION((void setML2s(const std::vector< double > &ml2s)))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML2s(vec);
  TEST_EQUAL(tc.getML2s()== vec,true)
END_SECTION

START_SECTION((void setML3s(const std::vector< double > &ml3s)))
  std::vector<double> vec;
  vec.push_back(0.1);
  vec.push_back(0.3);
  tc.setML3s(vec);
  TEST_EQUAL(tc.getML3s()== vec,true)
END_SECTION
  


START_SECTION((template<typename PeakType> void pickAndCalibrate(MSExperiment< Peak1D > &calib_spectra, MSExperiment< PeakType > &exp, std::vector< double > &exp_masses)))
  TOLERANCE_ABSOLUTE(0.000001)
  std::cout.precision(writtenDigits<DoubleReal>());
  MSExperiment<Peak1D> calib_exp;
  MSExperiment<Peak1D> exp,res_exp;
  MzDataFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_calibrants.mzData"),calib_exp);
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test.mzData"),exp);
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_output.mzData"),res_exp);
  vector<double> ref_masses;
  TextFile ref_file;

  ref_file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_calibrant_masses.txt"),true);

  for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
	{
		ref_masses.push_back(String(iter->c_str()).toDouble());
	}

  std::vector<double> ml1;
                                 
  ml1.push_back(418327.924993827);
                
  std::vector<double> ml2;
                
  ml2.push_back(253.645187196031);
  std::vector<double> ml3;
                
  ml3.push_back(-0.0414243465397252);
  tc.setML1s(ml1);
  tc.setML2s(ml2);
  tc.setML3s(ml3);

  Param param;
  param.setValue("PeakPicker:thresholds:peak_bound",400.0);
  param.setValue("PeakPicker:thresholds:correlation",0.0);
  param.setValue("PeakPicker:signal_to_noise",3.0);
  param.setValue("PeakPicker:centroid_percentage",0.6);
  tc.setParameters(param);
  tc.pickAndCalibrate(calib_exp,exp,ref_masses);
	
TOLERANCE_ABSOLUTE(0.01)
  TEST_EQUAL(exp.size()==res_exp.size(),true)
	for (Size i=0; i<exp.size(); ++i)
	{
		for (Size j=0; j<exp[i].size(); ++j)
		{
			TEST_REAL_SIMILAR(exp[i][j].getPos(),res_exp[i][j].getPos())
			TEST_REAL_SIMILAR(exp[i][j].getIntensity(),res_exp[i][j].getIntensity())
		}
	}
END_SECTION

tc = TOFCalibration();

START_SECTION((template<typename PeakType> void calibrate(MSExperiment<Peak1D> &calib_spectra, MSExperiment< PeakType > &exp, std::vector< double > &exp_masses)))
  TOLERANCE_ABSOLUTE(0.000001)
	std::cout.precision(writtenDigits<>(DoubleReal()));
  MSExperiment<> calib_exp;
  MSExperiment<> exp,res_exp;
  MzDataFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_calibrants2.mzData"),calib_exp);
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test2.mzData"),exp);
  file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_output2.mzData"),res_exp);
  vector<double> ref_masses;
  TextFile ref_file;

  ref_file.load(OPENMS_GET_TEST_DATA_PATH("TOFCalibration_test_calibrant_masses.txt"),true);

  for(TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
   {
     ref_masses.push_back(atof(iter->c_str()));
   }

  std::vector<double> ml1;
                                 
  ml1.push_back(418327.924993827);
  ml1.push_back(418257.238180361);
  ml1.push_back(418295.348979040);
  std::vector<double> ml2;
                
  ml2.push_back(253.645187196031);
  ml2.push_back(250.532666867861);
  ml2.push_back(251.878402283764);
  std::vector<double> ml3;
                
  ml3.push_back(-0.0414243465397252);
  ml3.push_back(-0.0428127107041497);
  ml3.push_back(-0.0419329877166861);
  tc.setML1s(ml1);
  tc.setML2s(ml2);
  tc.setML3s(ml3);

  tc.calibrate(calib_exp,exp,ref_masses);

	TOLERANCE_ABSOLUTE(0.01)
  TEST_EQUAL(exp.size()==res_exp.size(),true)
	for (Size i=0; i<exp.size(); ++i)
	{
		for (Size j=0; j<exp[i].size(); ++j)
		{
			TEST_REAL_SIMILAR(res_exp[i][j].getPos(),exp[i][j].getPos())
			TEST_REAL_SIMILAR(res_exp[i][j].getIntensity(),exp[i][j].getIntensity())
		}
	}

END_SECTION	

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
