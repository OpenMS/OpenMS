// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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


CHECK((PeakPicker(const PeakPicker& pp)))
  Param param;
  param.setValue("thresholds:signal_to_noise",7.0);
  param.setValue("thresholds:peak_bound",100.0);
  param.setValue("thresholds:peak_bound_ms2_level",10.0);
  param.setValue("thresholds:fwhm_bound",0.5);
  PeakPicker p;
  p.setParameters(param);
  
  PeakPicker p_copy(p);
  TEST_EQUAL(p_copy.getParameters() == param, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



