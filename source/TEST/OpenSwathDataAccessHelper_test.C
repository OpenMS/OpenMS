// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest, $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathDataAccessHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathDataAccessHelper* ptr = 0;
OpenSwathDataAccessHelper* nullPointer = 0;

START_SECTION(OpenSwathDataAccessHelper())
{
  ptr = new OpenSwathDataAccessHelper();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathDataAccessHelper())
{
  delete ptr;
}
END_SECTION

START_SECTION(OpenSwathDataAccessHelper::convertToSpectrumPtr(sptr))
{

  MSSpectrum<> sptr,omsptr;
  Peak1D p1;
  p1.setIntensity(1.0f);
  p1.setMZ(2.0);

  Peak1D p2;
  p2.setIntensity(2.0f);
  p2.setMZ(10.0);

  Peak1D p3;
  p3.setIntensity(3.0f);
  p3.setMZ(30.0);

  TEST_STRING_EQUAL(sptr.getName(), "")
  sptr.setName("my_fancy_name");
  sptr.push_back(p1);
  sptr.push_back(p2);
  sptr.push_back(p3);
  OpenSwath::SpectrumPtr p = OpenSwathDataAccessHelper::convertToSpectrumPtr(sptr);
  OpenSwathDataAccessHelper::convertToOpenMSSpectrum(p,omsptr);
  TEST_REAL_SIMILAR(p->getMZArray()->data[0],2.0);
  TEST_REAL_SIMILAR(p->getMZArray()->data[1],10.0);
  TEST_REAL_SIMILAR(p->getMZArray()->data[2],30.0);


  TEST_REAL_SIMILAR(p->getIntensityArray()->data[0],1.0f);
  TEST_REAL_SIMILAR(p->getIntensityArray()->data[1],2.0f);
  TEST_REAL_SIMILAR(p->getIntensityArray()->data[2],3.0f);
}
END_SECTION

START_SECTION(OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram,cptr))
{
  //void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(OpenMS::MSChromatogram<> & chromatogram,
  //                                                          const OpenSwath::ChromatogramPtr cptr)
  OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram());
  cptr->getTimeArray()->data.push_back(1.0);
  cptr->getTimeArray()->data.push_back(2.0);
  cptr->getTimeArray()->data.push_back(3.0);
  cptr->getTimeArray()->data.push_back(4.0);

  cptr->getIntensityArray()->data.push_back(4.0);
  cptr->getIntensityArray()->data.push_back(3.0);
  cptr->getIntensityArray()->data.push_back(2.0);
  cptr->getIntensityArray()->data.push_back(1.0);

  MSChromatogram<> chromatogram;
  OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram,cptr);

  TEST_REAL_SIMILAR(chromatogram[0].getRT(),1.);
  TEST_REAL_SIMILAR(chromatogram[0].getIntensity(),4.);
  TEST_REAL_SIMILAR(chromatogram[1].getRT(),2.);
  TEST_REAL_SIMILAR(chromatogram[1].getIntensity(),3.);
  TEST_REAL_SIMILAR(chromatogram[2].getRT(),3.);
  TEST_REAL_SIMILAR(chromatogram[2].getIntensity(),2.);

}
END_SECTION

START_SECTION(convertToOpenMSSpectrum(spectrum,sptr))
{
  OpenSwath::SpectrumPtr cptr(new OpenSwath::Spectrum());
  cptr->getMZArray()->data.push_back(1.0);
  cptr->getMZArray()->data.push_back(2.0);
  cptr->getMZArray()->data.push_back(3.0);
  cptr->getMZArray()->data.push_back(4.0);

  cptr->getIntensityArray()->data.push_back(4.0);
  cptr->getIntensityArray()->data.push_back(3.0);
  cptr->getIntensityArray()->data.push_back(2.0);
  cptr->getIntensityArray()->data.push_back(1.0);

  MSSpectrum<> spectrum;
  OpenSwathDataAccessHelper::convertToOpenMSSpectrum(cptr, spectrum);

  TEST_REAL_SIMILAR(spectrum[0].getMZ(),1.);
  TEST_REAL_SIMILAR(spectrum[0].getIntensity(),4.);
  TEST_REAL_SIMILAR(spectrum[1].getMZ(),2.);
  TEST_REAL_SIMILAR(spectrum[1].getIntensity(),3.);
  TEST_REAL_SIMILAR(spectrum[2].getMZ(),3.);
  TEST_REAL_SIMILAR(spectrum[2].getIntensity(),2.);
}
END_SECTION


START_SECTION(convertTargetedExp(transition_exp_, transition_exp))
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



