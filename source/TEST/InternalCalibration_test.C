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
// $Maintainer: Alexandra Zerck$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzDataFile.h>
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


//START_SECTION((template <typename InputPeakType> void calibrate(MSExperiment< InputPeakType > &exp, std::vector< double > &ref_masses)))
//{
//     TOLERANCE_ABSOLUTE(0.000001)
//   MSExperiment<Peak1D> exp;
//   MSExperiment<> exp_peaks;
//   MzDataFile file;
//   file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_test.mzData"),exp);
//   std::vector<double> ref_masses;
//   ref_masses.push_back(1296.68476942);
//   ref_masses.push_back(2465.19833942);
	
//   Param param;
// 	param.setValue("PeakPicker:peak_width",0.15);
//   param.setValue("PeakPicker:thresholds:peak_bound",800.0);
//   param.setValue("PeakPicker:fwhm_bound_factor",0.0);
//   param.setValue("PeakPicker:thresholds:correlation",0.0);
// //  param.setValue("PeakPicker:centroid_percentage",0.6);
//   ptr->setParameters(param);
//   ptr->calibrate(exp,ref_masses);
	
//   PeakPickerCWT pp;
//   Param pp_param;
// 	param.setValue("peak_width",0.15);
// 	pp_param.setValue("thresholds:correlation",0.0);
//   pp_param.setValue("thresholds:peak_bound",800.0);
//   pp_param.setValue("fwhm_bound_factor",0.0);
//   pp.setParameters(pp_param);
//   pp.pickExperiment(exp,exp_peaks);
//   Peak1D peak;
//   peak.setMZ(1296.68476942);
//   MSExperiment<>::SpectrumType::Iterator it = lower_bound(exp_peaks[0].begin(),exp_peaks[0].end(),peak,Peak1D::PositionLess());
//   --it;
//   TEST_REAL_SIMILAR(it->getMZ(),1296.68476942)
//   peak.setMZ(2465.19833942);
//   it = lower_bound(exp_peaks[0].begin(),exp_peaks[0].end(),peak,Peak1D::PositionLess());
//   --it;
//   TEST_REAL_SIMILAR(it->getMZ(),2465.19833942)
// }
//END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


