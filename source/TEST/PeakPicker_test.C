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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPicker, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPicker* ptr = 0;
CHECK((PeakPicker()))
  ptr = new PeakPicker();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~PeakPicker()))
  delete ptr;
RESULT

CHECK((PeakPicker& operator=(const PeakPicker& pp)))
  Param param;
  param.setValue("thresholds:signal_to_noise",7);
  param.setValue("thresholds:peak_bound",100);
  param.setValue("thresholds:peak_bound_ms2_level",10);
  param.setValue("thresholds:fwhm_bound",0.5);
  PeakPicker p;
  p.setParam(param);
  
  PeakPicker p_copy;
  p_copy = p;
  TEST_REAL_EQUAL(p_copy.getPeakBound(),100)
  TEST_REAL_EQUAL(p_copy.getPeakBoundMs2Level(),10)
  TEST_REAL_EQUAL(p_copy.getSignalToNoiseLevel(),7)
  TEST_REAL_EQUAL(p_copy.getFwhmBound(),0.5)
  TEST_EQUAL(p_copy.getParam() == param, true)
RESULT

CHECK((void setParam(Param param)))
  Param param;
  param.setValue("thresholds:signal_to_noise",7);
  param.setValue("thresholds:peak_bound",100);
  param.setValue("thresholds:peak_bound_ms2_level",10);
  param.setValue("thresholds:fwhm_bound",0.5);
  
  PeakPicker pp;
  pp.setParam(param);
  
  TEST_REAL_EQUAL(pp.getPeakBound(),100)
  TEST_REAL_EQUAL(pp.getPeakBoundMs2Level(),10)
  TEST_REAL_EQUAL(pp.getSignalToNoiseLevel(),7)
  TEST_REAL_EQUAL(pp.getFwhmBound(),0.5)
  TEST_EQUAL(pp.getParam() == param, true)
RESULT

CHECK((PeakPicker(const PeakPicker& pp)))
  Param param;
  param.setValue("thresholds:signal_to_noise",7);
  param.setValue("thresholds:peak_bound",100);
  param.setValue("thresholds:peak_bound_ms2_level",10);
  param.setValue("thresholds:fwhm_bound",0.5);
  PeakPicker p;
  p.setParam(param);
  
  PeakPicker p_copy(p);
  TEST_REAL_EQUAL(p_copy.getPeakBound(),100)
  TEST_REAL_EQUAL(p_copy.getPeakBoundMs2Level(),10)
  TEST_REAL_EQUAL(p_copy.getSignalToNoiseLevel(),7)
  TEST_REAL_EQUAL(p_copy.getFwhmBound(),0.5)
  TEST_EQUAL(p_copy.getParam() == param, true)
RESULT

CHECK((const Param& getParam() const))
  Param param;
  param.setValue("thresholds:signal_to_noise",7);
  param.setValue("thresholds:peak_bound",100);
  param.setValue("thresholds:peak_bound_ms2_level",10);
  param.setValue("thresholds:fwhm_bound",0.5);
  
  PeakPicker pp;
  pp.setParam(param);

  TEST_REAL_EQUAL((double)pp.getParam().getValue("thresholds:signal_to_noise"),7.0);
  TEST_REAL_EQUAL((double)pp.getParam().getValue("thresholds:peak_bound"),100.0);
  TEST_REAL_EQUAL((double)pp.getParam().getValue("thresholds:peak_bound_ms2_level"),10.0);
  TEST_REAL_EQUAL((double)pp.getParam().getValue("thresholds:fwhm_bound"),0.5);
RESULT

CHECK((const float& getFwhmBound() const))
  PeakPicker pp;
  pp.setParam(Param());
  TEST_REAL_EQUAL(pp.getFwhmBound(),0.2)
RESULT

CHECK((const float& getPeakBound() const))
  PeakPicker pp;
  pp.setParam(Param());
  TEST_REAL_EQUAL(pp.getPeakBound(),200)
RESULT

CHECK((const float& getPeakBoundMs2Level() const))
  PeakPicker pp;
  pp.setParam(Param());
  TEST_REAL_EQUAL(pp.getPeakBoundMs2Level(),50)
RESULT

CHECK((const float& getSignalToNoiseLevel() const))
  PeakPicker pp;
  pp.setParam(Param());
  TEST_REAL_EQUAL(pp.getSignalToNoiseLevel(),3)
RESULT

CHECK((void setFwhmBound(const float& fwhm)))
  PeakPicker pp;
  
  pp.setFwhmBound(.3);
  TEST_REAL_EQUAL(pp.getFwhmBound(),0.3)
RESULT

CHECK((void setParam(const Param& param)))
  // ???
RESULT

CHECK((void setPeakBound(const float& peak_bound)))
  PeakPicker pp;
  
  pp.setPeakBound(1000);
  TEST_REAL_EQUAL(pp.getPeakBound(),1000)
RESULT

CHECK((void setPeakBoundMs2Level(const float& peak_bound_ms2_level)))
  PeakPicker pp;
  
  pp.setPeakBoundMs2Level(10);
  TEST_REAL_EQUAL(pp.getPeakBoundMs2Level(),10)
RESULT

CHECK((void setSignalToNoiseLevel(const float& signal_to_noise)))
  PeakPicker pp;
    
  pp.setSignalToNoiseLevel(10);
  TEST_REAL_EQUAL(pp.getSignalToNoiseLevel(),10)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



