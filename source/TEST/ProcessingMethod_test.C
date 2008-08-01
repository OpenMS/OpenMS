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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/ProcessingMethod.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProcessingMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProcessingMethod* ptr = 0;
CHECK(ProcessingMethod())
	ptr = new ProcessingMethod();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ProcessingMethod())
	delete ptr;
RESULT

CHECK(SpectrumSettings::SpectrumType getSpectrumType() const)
  ProcessingMethod tmp;
  TEST_EQUAL(tmp.getSpectrumType(),SpectrumSettings::UNKNOWN);
RESULT

CHECK(void setSpectrumType(SpectrumSettings::SpectrumType method))
  ProcessingMethod tmp;
  tmp.setSpectrumType(SpectrumSettings::PEAKS);
  TEST_EQUAL(tmp.getSpectrumType(),SpectrumSettings::PEAKS);
RESULT

CHECK(bool getChargeDeconvolution() const)
  ProcessingMethod tmp;
  TEST_EQUAL(tmp.getChargeDeconvolution(),false);  
RESULT

CHECK(void setChargeDeconvolution(bool charge_deconvolution))
  ProcessingMethod tmp;
  tmp.setChargeDeconvolution(true);
  TEST_EQUAL(tmp.getChargeDeconvolution(),true);  
RESULT

CHECK(bool getDeisotoping() const)
  ProcessingMethod tmp;
  TEST_EQUAL(tmp.getDeisotoping(),false);  
RESULT

CHECK(void setDeisotoping(bool deisotoping))
  ProcessingMethod tmp;
  tmp.setDeisotoping(true);
  TEST_EQUAL(tmp.getDeisotoping(),true);  
RESULT

CHECK(float getIntensityCutoff() const)
  ProcessingMethod tmp;
  TEST_REAL_EQUAL(tmp.getIntensityCutoff(), 0);
RESULT

CHECK(void setIntensityCutoff(float cutoff))
  ProcessingMethod tmp;
  tmp.setIntensityCutoff(22.6);
  TEST_REAL_EQUAL(tmp.getIntensityCutoff(), 22.6);
RESULT

CHECK(ProcessingMethod& operator= (const ProcessingMethod& source))
  ProcessingMethod tmp;
  tmp.setChargeDeconvolution(true);
  tmp.setDeisotoping(true);
  tmp.setSpectrumType(SpectrumSettings::PEAKS);
  tmp.setIntensityCutoff(3.4);
  tmp.setMetaValue("label",String("label"));
  
  ProcessingMethod tmp2(tmp);
  TEST_EQUAL(tmp2.getChargeDeconvolution(),true);
  TEST_EQUAL(tmp2.getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSpectrumType(),SpectrumSettings::PEAKS);
  TEST_REAL_EQUAL(tmp2.getIntensityCutoff(), 3.4);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
RESULT

CHECK(ProcessingMethod(const ProcessingMethod& source))
  ProcessingMethod tmp;
  tmp.setChargeDeconvolution(true);
  tmp.setDeisotoping(true);
  tmp.setSpectrumType(SpectrumSettings::PEAKS);
  tmp.setIntensityCutoff(2.8);
  tmp.setMetaValue("label",String("label"));
  
  ProcessingMethod tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getChargeDeconvolution(),true);
  TEST_EQUAL(tmp2.getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSpectrumType(),SpectrumSettings::PEAKS);
  TEST_REAL_EQUAL(tmp2.getIntensityCutoff(), 2.8);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  
  tmp2 = ProcessingMethod();
  TEST_EQUAL(tmp2.getChargeDeconvolution(),false);
  TEST_EQUAL(tmp2.getDeisotoping(),false);
  TEST_EQUAL(tmp2.getSpectrumType(),SpectrumSettings::UNKNOWN);
  TEST_REAL_EQUAL(tmp2.getIntensityCutoff(), 0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const ProcessingMethod& rhs) const)
  ProcessingMethod edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
  edit.setChargeDeconvolution(true);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setDeisotoping(true);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setSpectrumType(SpectrumSettings::PEAKS);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setIntensityCutoff(99.24);
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK(bool operator!= (const ProcessingMethod& rhs) const)
  ProcessingMethod edit, empty;
  
  TEST_EQUAL(edit!=empty, false);
  
  edit.setChargeDeconvolution(true);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.setDeisotoping(true);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.setSpectrumType(SpectrumSettings::PEAKS);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.setIntensityCutoff(99.24);
  TEST_EQUAL(edit!=empty, true);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty, true);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



