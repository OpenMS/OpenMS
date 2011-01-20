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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/HPLC.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HPLC, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HPLC* ptr = 0;
START_SECTION(HPLC())
	ptr = new HPLC();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~HPLC())
	delete ptr;
END_SECTION

START_SECTION(Gradient& getGradient())
  HPLC tmp;
  TEST_EQUAL(tmp.getGradient().getEluents().size(),0);
END_SECTION

START_SECTION(void setGradient(const Gradient& gradient))
  HPLC tmp;
  Gradient g;
  g.addEluent("A");
  tmp.setGradient(g);
  TEST_EQUAL(tmp.getGradient().getEluents().size(),1);
  TEST_EQUAL(tmp.getGradient().getEluents()[0],"A");
END_SECTION

START_SECTION(const Gradient& getGradient() const)
  HPLC tmp;
  tmp.getGradient().addEluent("A");
  TEST_EQUAL(tmp.getGradient().getEluents().size(),1);
  TEST_EQUAL(tmp.getGradient().getEluents()[0],"A");
END_SECTION

START_SECTION(UInt getFlux() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getFlux(),0);  
END_SECTION

START_SECTION(void setFlux(UInt flux))
  HPLC tmp;
  tmp.setFlux(5);
  TEST_EQUAL(tmp.getFlux(),5);  
END_SECTION

START_SECTION(UInt getPressure() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getPressure(),0);  
END_SECTION

START_SECTION(void setPressure(UInt pressure))
  HPLC tmp;
  tmp.setPressure(5);
  TEST_EQUAL(tmp.getPressure(),5);  
END_SECTION

START_SECTION(Int getTemperature() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getTemperature(),21);  
END_SECTION

START_SECTION(void setTemperature(Int temperature))
  HPLC tmp;
  tmp.setTemperature(5);
  TEST_EQUAL(tmp.getTemperature(),5);  
END_SECTION

START_SECTION(String getComment() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getComment(),"");  
END_SECTION

START_SECTION(void setComment(String comment))
  HPLC tmp;
  tmp.setComment("comment");
  TEST_EQUAL(tmp.getComment(),"comment");  
END_SECTION

START_SECTION(const String& getInstrument() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getInstrument(),"");  
END_SECTION

START_SECTION(void setInstrument(const String& instrument))
  HPLC tmp;
  tmp.setInstrument("instrument");
  TEST_EQUAL(tmp.getInstrument(),"instrument");  
END_SECTION

START_SECTION(const String& getColumn() const)
  HPLC tmp;
  TEST_EQUAL(tmp.getColumn(),"");  
END_SECTION

START_SECTION(void setColumn(const String& column))
  HPLC tmp;
  tmp.setColumn("column");
  TEST_EQUAL(tmp.getColumn(),"column");  
END_SECTION

START_SECTION(HPLC(const HPLC& source))
  HPLC tmp;
  tmp.setInstrument("instrument");
  tmp.setComment("comment");
  tmp.setColumn("column");
  tmp.setPressure(5);
  tmp.setFlux(6);
  tmp.setTemperature(7);
  
  HPLC tmp2(tmp);
	TEST_EQUAL(tmp2.getInstrument(),"instrument");
  TEST_EQUAL(tmp2.getComment(),"comment");
  TEST_EQUAL(tmp2.getColumn(),"column");
  TEST_EQUAL(tmp2.getPressure(),5);
  TEST_EQUAL(tmp2.getFlux(),6);
  TEST_EQUAL(tmp2.getTemperature(),7);
END_SECTION

START_SECTION(HPLC& operator = (const HPLC& source))
  HPLC tmp;
  tmp.setInstrument("instrument");
  tmp.setComment("comment");
  tmp.setColumn("column");
  tmp.setPressure(5);
  tmp.setFlux(6);
  tmp.setTemperature(7);
  
  HPLC tmp2(tmp);
	TEST_EQUAL(tmp2.getInstrument(),"instrument");
  TEST_EQUAL(tmp2.getComment(),"comment");
  TEST_EQUAL(tmp2.getColumn(),"column");
  TEST_EQUAL(tmp2.getPressure(),5);
  TEST_EQUAL(tmp2.getFlux(),6);
  TEST_EQUAL(tmp2.getTemperature(),7);
  
  tmp2 = HPLC();
	TEST_EQUAL(tmp2.getInstrument(),"");
  TEST_EQUAL(tmp2.getComment(),"");
  TEST_EQUAL(tmp2.getColumn(),"");
  TEST_EQUAL(tmp2.getPressure(),0);
  TEST_EQUAL(tmp2.getFlux(),0);
  TEST_EQUAL(tmp2.getTemperature(),21);  
END_SECTION

START_SECTION(bool operator == (const HPLC& source) const)
  HPLC edit, empty;
  
  TEST_EQUAL(edit==empty,true);  
  
  edit.setInstrument("instrument");
  TEST_EQUAL(edit==empty,false); 
  
  edit = empty;
  edit.setComment("comment");
  TEST_EQUAL(edit==empty,false); 
  
  edit = empty;
  edit.setColumn("column");
  TEST_EQUAL(edit==empty,false); 
  
  edit = empty;
  edit.setPressure(5);
  TEST_EQUAL(edit==empty,false); 
  
  edit = empty;
  edit.setFlux(6);
  TEST_EQUAL(edit==empty,false); 
  
  edit = empty;
  edit.setTemperature(7);
  TEST_EQUAL(edit==empty,false); 
END_SECTION

START_SECTION(bool operator != (const HPLC& source) const)
  HPLC edit, empty;
  
  TEST_EQUAL(edit!=empty,false);  
  
  edit.setInstrument("instrument");
  TEST_EQUAL(edit!=empty,true); 
  
  edit = empty;
  edit.setComment("comment");
  TEST_EQUAL(edit!=empty,true); 
  
  edit = empty;
  edit.setColumn("column");
  TEST_EQUAL(edit!=empty,true); 
  
  edit = empty;
  edit.setPressure(5);
  TEST_EQUAL(edit!=empty,true); 
  
  edit = empty;
  edit.setFlux(6);
  TEST_EQUAL(edit!=empty,true); 
  
  edit = empty;
  edit.setTemperature(7);
  TEST_EQUAL(edit!=empty,true); 
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



