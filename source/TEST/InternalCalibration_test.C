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
// $Maintainer: Alexandra Zerck$
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
CHECK(InternalCalibration())
{
	ptr = new InternalCalibration();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~InternalCalibration())
{
	delete ptr;
}
RESULT

ptr = new InternalCalibration();

CHECK((InternalCalibration(InternalCalibration &obj)))
{
  InternalCalibration copy(*ptr);
  TEST_EQUAL(copy.getPeaks()== ptr->getPeaks(),true)
  TEST_EQUAL(copy.getMonoisotopicPeaks()==ptr->getMonoisotopicPeaks(),true )
}
RESULT

CHECK((InternalCalibration& operator=(const InternalCalibration &obj)))
{
  InternalCalibration copy;
  copy = *ptr;
  TEST_EQUAL(copy.getPeaks()== ptr->getPeaks(),true)
  TEST_EQUAL(copy.getMonoisotopicPeaks()==ptr->getMonoisotopicPeaks(),true )

}
RESULT

CHECK((template <typename InputPeakType> void calibrate(MSExperiment< InputPeakType > &exp, std::vector< double > &ref_masses, bool peak_data=false)))
{
    PRECISION(0.000001)
  MSExperiment<RawDataPoint1D> exp;
  MSExperiment<PickedPeak1D> exp_peaks;
  MzDataFile file;
  file.load("data/InternalCalibration_test.mzData",exp);
  std::vector<double> ref_masses;
  ref_masses.push_back(1296.68476942);
  ref_masses.push_back(2465.19833942);
	
  Param param;
  param.setValue("PeakPicker:thresholds:peak_bound",800);
  param.setValue("PeakPicker:thresholds:fwhm_bound",0.0);
  param.setValue("PeakPicker:thresholds:correlation",0.0);
  ptr->setParameters(param);
  ptr->calibrate(exp,ref_masses,false);
	
  PeakPickerCWT pp;
  pp.setPeakCorrBound(0.0);
  pp.setPeakBound(800);
  pp.setFwhmBound(0.0);
  pp.pickExperiment(exp,exp_peaks);
  PickedPeak1D peak;
  peak.setMZ(1296.68476942);
  MSSpectrum<PickedPeak1D>::Iterator it = lower_bound(exp_peaks[0].begin(),exp_peaks[0].end(),peak,RawDataPoint1D::PositionLess());
  --it;
  TEST_REAL_EQUAL(it->getMZ(),1296.68476942)
  peak.setMZ(2465.19833942);
  it = lower_bound(exp_peaks[0].begin(),exp_peaks[0].end(),peak,RawDataPoint1D::PositionLess());
  --it;
  TEST_REAL_EQUAL(it->getMZ(),2465.19833942)
}
RESULT

CHECK((const MSExperiment<PickedPeakType>& getPeaks() const))
{
  MSExperiment<PickedPeak1D> exp;
  MSSpectrum<PickedPeak1D> spec;
  PickedPeak1D peak;
  peak.setMZ(100.1);
  spec.push_back(peak);
  peak.setMZ(102.1);
  spec.push_back(peak);
  exp.push_back(spec);
  ptr->setPeaks(exp);
  TEST_EQUAL(ptr->getPeaks()== exp,true)
}
RESULT

CHECK((void setPeaks(const MSExperiment< PickedPeakType > &exp_peaks)))
{
  MSExperiment<PickedPeak1D> exp;
  MSSpectrum<PickedPeak1D> spec;
  PickedPeak1D peak;
  peak.setMZ(100.1);
  spec.push_back(peak);
  peak.setMZ(102.1);
  spec.push_back(peak);
  exp.push_back(spec);
  
  ptr->setPeaks(exp);
  TEST_EQUAL(ptr->getPeaks()==exp,true)

}
RESULT

CHECK((const std::vector<std::vector<UInt> >& getMonoisotopicPeaks() const))
{
  std::vector<std::vector<UInt> > p;
  std::vector<UInt> vec;
  vec.push_back(1);
  vec.push_back(2);
  vec.push_back(3);
  p.push_back(vec);
  
  p.push_back(vec);
  ptr->setMonoisotopicPeaks(p);
  TEST_EQUAL(ptr->getMonoisotopicPeaks()== p,true)
}
RESULT

CHECK((void setMonoisotopicPeaks(const std::vector< std::vector< UInt > > &monoiso_peaks)))
{
  std::vector<std::vector<UInt> > p;
  std::vector<UInt> vec;
  vec.push_back(1);
  vec.push_back(2);
  vec.push_back(3);
  p.push_back(vec);
  ptr->setMonoisotopicPeaks(p);
  TEST_EQUAL(ptr->getMonoisotopicPeaks()==p,true)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



