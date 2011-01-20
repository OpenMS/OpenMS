// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InternalCalibration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InternalCalibration* ptr = 0;
START_SECTION(InternalCalibration())
{
	ptr = new InternalCalibration();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~InternalCalibration())
{
	delete ptr;
}
END_SECTION

ptr = new InternalCalibration();

MSExperiment<> exp;
MzDataFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_test.mzData"),exp);
std::vector<double> ref_masses;
ref_masses.push_back(1296.68476942);
ref_masses.push_back(2465.19833942);
Param param;
param.setValue("mz_tolerance",100.);
param.setValue("mz_tolerance_unit","ppm");

START_SECTION((template < typename InputPeakType > void calibrateMapSpectrumwise(const MSExperiment< InputPeakType > &exp, MSExperiment< InputPeakType > &calibrated_exp, std::vector< DoubleReal > &ref_masses)))
{
  TOLERANCE_ABSOLUTE(0.000001)
  MSExperiment<> calibrated_exp;
  
  ptr->setParameters(param);
  ptr->calibrateMapSpectrumwise(exp,calibrated_exp,ref_masses);
	
  TEST_REAL_SIMILAR(calibrated_exp[0][14].getMZ(),1296.68476942)
  TEST_REAL_SIMILAR(calibrated_exp[0][77].getMZ(),2465.19833942)

}
END_SECTION

START_SECTION((template < typename InputPeakType > void calibrateMapGlobally(const MSExperiment< InputPeakType > &exp, MSExperiment< InputPeakType > &calibrated_exp, std::vector< DoubleReal > &ref_masses, String trafo_file_name="")))
{
  TOLERANCE_ABSOLUTE(0.000001)
  MSExperiment<> calibrated_exp;
  ptr->setParameters(param);
  ptr->calibrateMapGlobally(exp,calibrated_exp,ref_masses);
  
  TEST_REAL_SIMILAR(calibrated_exp[0][14].getMZ(),1296.68476942)
  TEST_REAL_SIMILAR(calibrated_exp[1][40].getMZ(),1296.68476942)
  TEST_REAL_SIMILAR(calibrated_exp[0][77].getMZ(),2465.19833942)
  TEST_REAL_SIMILAR(calibrated_exp[1][90].getMZ(),2465.19833942)
}
END_SECTION
IdXMLFile id_file;
std::vector<ProteinIdentification> prot_ids;
std::vector<PeptideIdentification> pep_ids;
id_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_1.IdXML"),prot_ids,pep_ids);
START_SECTION((template < typename InputPeakType > void calibrateMapGlobally(const MSExperiment< InputPeakType > &exp, MSExperiment< InputPeakType > &calibrated_exp, std::vector< PeptideIdentification > &ref_ids, String trafo_file_name="")))
{
  TOLERANCE_ABSOLUTE(0.000001)
  MSExperiment<> calibrated_exp;
  ptr->setParameters(param);
  ptr->calibrateMapGlobally(exp,calibrated_exp,pep_ids);
  

	TEST_REAL_SIMILAR(calibrated_exp[0][14].getMZ(),1296.68476942)
  TEST_REAL_SIMILAR(calibrated_exp[1][40].getMZ(),1296.68476942)
  TEST_REAL_SIMILAR(calibrated_exp[0][77].getMZ(),2465.19833942)
  TEST_REAL_SIMILAR(calibrated_exp[1][90].getMZ(),2465.19833942)
}
END_SECTION

FeatureMap<> f_map;
FeatureXMLFile f_file;
f_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_annotated.featureXML"),f_map);
START_SECTION((void calibrateMapGlobally(const FeatureMap<> &feature_map, FeatureMap<> &calibrated_feature_map, String trafo_file_name="")))
{
  FeatureMap<> calibrated_f_map;
  ptr->calibrateMapGlobally(f_map,calibrated_f_map);
  TEST_REAL_SIMILAR(calibrated_f_map[0].getMZ(),687.841430243171)
  TEST_REAL_SIMILAR(calibrated_f_map[1].getMZ(),720.005082366204)
  TEST_REAL_SIMILAR(calibrated_f_map[2].getMZ(),927.493444113771)
  TEST_REAL_SIMILAR(calibrated_f_map[3].getMZ(),1052.06529617992)
  TEST_REAL_SIMILAR(calibrated_f_map[4].getMZ(),1224.59976809287)
  TEST_REAL_SIMILAR(calibrated_f_map[5].getMZ(),998.486309862771)
}
END_SECTION
id_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_2.IdXML"),prot_ids,pep_ids);
START_SECTION((void calibrateMapGlobally(const FeatureMap<> &feature_map, FeatureMap<> &calibrated_feature_map, std::vector< PeptideIdentification > &ref_ids, String trafo_file_name="")))
{
  FeatureMap<> calibrated_f_map;
  ptr->calibrateMapGlobally(f_map,calibrated_f_map,pep_ids);
  TEST_REAL_SIMILAR(calibrated_f_map[0].getMZ(),687.841430243171)
  TEST_REAL_SIMILAR(calibrated_f_map[1].getMZ(),720.005082366204)
  TEST_REAL_SIMILAR(calibrated_f_map[2].getMZ(),927.493444113771)
  TEST_REAL_SIMILAR(calibrated_f_map[3].getMZ(),1052.06529617992)
  TEST_REAL_SIMILAR(calibrated_f_map[4].getMZ(),1224.59976809287)
  TEST_REAL_SIMILAR(calibrated_f_map[5].getMZ(),998.486309862771)
}
END_SECTION

START_SECTION((template < typename InputPeakType > void calibrateMapList(std::vector< MSExperiment< InputPeakType > > &exp_list, std::vector< MSExperiment< InputPeakType > > &calibrated_exp_list, std::vector< DoubleReal > &ref_masses, std::vector< DoubleReal > &detected_background_masses)))
{
  NOT_TESTABLE  // not yet existing
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


