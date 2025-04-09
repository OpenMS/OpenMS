// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/Gradient.h>
#include <OpenMS/DATASTRUCTURES/String.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Gradient, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Gradient* ptr = nullptr;
Gradient* nullPointer = nullptr;
START_SECTION((Gradient()))
	ptr = new Gradient();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~Gradient()))
	delete ptr;
END_SECTION

START_SECTION((const std::vector<String>& getEluents() const))
  Gradient tmp;
  TEST_EQUAL(tmp.getEluents().size(),0);
END_SECTION

START_SECTION((void addEluent(const String& eluent) ))
  Gradient tmp;

  tmp.addEluent("A");
  TEST_EQUAL(tmp.getEluents().size(),1);
  TEST_EQUAL(tmp.getEluents()[0],"A");

  tmp.addEluent("B");
  TEST_EQUAL(tmp.getEluents().size(),2);
  TEST_EQUAL(tmp.getEluents()[0],"A");
  TEST_EQUAL(tmp.getEluents()[1],"B");
  
  TEST_EXCEPTION(Exception::InvalidValue,tmp.addEluent("B"))
END_SECTION

START_SECTION((void clearEluents()))
  Gradient tmp;

  tmp.addEluent("A");
  tmp.clearEluents();
  TEST_EQUAL(tmp.getEluents().size(),0);
END_SECTION

START_SECTION((const std::vector<Int>& getTimepoints() const))
  Gradient tmp;
  TEST_EQUAL(tmp.getTimepoints().size(),0);
END_SECTION

START_SECTION((void addTimepoint(Int timepoint) ))
  Gradient tmp;

  tmp.addTimepoint(5);
  TEST_EQUAL(tmp.getTimepoints().size(),1);
  TEST_EQUAL(tmp.getTimepoints()[0],5);

  tmp.addTimepoint(7);
  TEST_EQUAL(tmp.getTimepoints().size(),2);
  TEST_EQUAL(tmp.getTimepoints()[0],5);
  TEST_EQUAL(tmp.getTimepoints()[1],7);

  TEST_EXCEPTION(Exception::OutOfRange,tmp.addTimepoint(6));
  TEST_EXCEPTION(Exception::OutOfRange,tmp.addTimepoint(7));
  tmp.addTimepoint(8);
END_SECTION

START_SECTION((void clearTimepoints()))
  Gradient tmp;

  tmp.addTimepoint(5);
  tmp.clearTimepoints();
  TEST_EQUAL(tmp.getTimepoints().size(),0);
END_SECTION

START_SECTION((const std::vector< std::vector< UInt > > & getPercentages() const))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addTimepoint(7);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.addEluent("C");
  TEST_EQUAL(tmp.getPercentages().size(),3);
  TEST_EQUAL(tmp.getPercentages()[0].size(),2);
  TEST_EQUAL(tmp.getPercentages()[1].size(),2);
  TEST_EQUAL(tmp.getPercentages()[2].size(),2);
  TEST_EQUAL(tmp.getPercentages()[0][0],0);
  TEST_EQUAL(tmp.getPercentages()[0][1],0);
  TEST_EQUAL(tmp.getPercentages()[1][0],0);
  TEST_EQUAL(tmp.getPercentages()[1][1],0);
  TEST_EQUAL(tmp.getPercentages()[2][0],0);
  TEST_EQUAL(tmp.getPercentages()[2][1],0);
END_SECTION

START_SECTION((void setPercentage(const String& eluent, Int timepoint, UInt percentage) ))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addTimepoint(7);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.addEluent("C");
  
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  tmp.setPercentage("C",5,0);
  tmp.setPercentage("A",7,30);
  tmp.setPercentage("B",7,50);
  tmp.setPercentage("C",7,20);
  
  TEST_EQUAL(tmp.getPercentages().size(),3);
  TEST_EQUAL(tmp.getPercentages()[0].size(),2);
  TEST_EQUAL(tmp.getPercentages()[1].size(),2);
  TEST_EQUAL(tmp.getPercentages()[2].size(),2);
  TEST_EQUAL(tmp.getPercentages()[0][0],90);
  TEST_EQUAL(tmp.getPercentages()[0][1],30);
  TEST_EQUAL(tmp.getPercentages()[1][0],10);
  TEST_EQUAL(tmp.getPercentages()[1][1],50);
  TEST_EQUAL(tmp.getPercentages()[2][0],0);
  TEST_EQUAL(tmp.getPercentages()[2][1],20);
  
  TEST_EXCEPTION(Exception::InvalidValue,tmp.setPercentage("D",7,20));
  TEST_EXCEPTION(Exception::InvalidValue,tmp.setPercentage("C",9,20));
  TEST_EXCEPTION(Exception::InvalidValue,tmp.setPercentage("C",7,101));
END_SECTION

START_SECTION((UInt getPercentage(const String& eluent, Int timepoint) const ))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addTimepoint(7);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.addEluent("C");
  
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  tmp.setPercentage("C",5,0);
  tmp.setPercentage("A",7,30);
  tmp.setPercentage("B",7,50);
  tmp.setPercentage("C",7,20);
  
  TEST_EQUAL(tmp.getPercentage("A",5),90);
  TEST_EQUAL(tmp.getPercentage("A",7),30);
  TEST_EQUAL(tmp.getPercentage("B",5),10);
  TEST_EQUAL(tmp.getPercentage("B",7),50);
  TEST_EQUAL(tmp.getPercentage("C",5),0);
  TEST_EQUAL(tmp.getPercentage("C",7),20);
  
  TEST_EXCEPTION(Exception::InvalidValue,tmp.getPercentage("D",7));
  TEST_EXCEPTION(Exception::InvalidValue,tmp.getPercentage("C",9));
END_SECTION

START_SECTION((void clearPercentages()))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addTimepoint(7);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.addEluent("C");
  
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  tmp.setPercentage("C",5,0);
  tmp.setPercentage("A",7,30);
  tmp.setPercentage("B",7,50);
  tmp.setPercentage("C",7,20);
  
  tmp.clearPercentages();
  
  TEST_EQUAL(tmp.getPercentages().size(),3);
  TEST_EQUAL(tmp.getPercentages()[0].size(),2);
  TEST_EQUAL(tmp.getPercentages()[1].size(),2);
  TEST_EQUAL(tmp.getPercentages()[2].size(),2);
  TEST_EQUAL(tmp.getPercentages()[0][0],0);
  TEST_EQUAL(tmp.getPercentages()[0][1],0);
  TEST_EQUAL(tmp.getPercentages()[1][0],0);
  TEST_EQUAL(tmp.getPercentages()[1][1],0);
  TEST_EQUAL(tmp.getPercentages()[2][0],0);
  TEST_EQUAL(tmp.getPercentages()[2][1],0);
END_SECTION

START_SECTION((bool isValid() const))
  Gradient tmp;
  TEST_EQUAL(tmp.isValid(),true);
  tmp.addTimepoint(5);
  tmp.addTimepoint(7);
  TEST_EQUAL(tmp.isValid(),false);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.addEluent("C");
  TEST_EQUAL(tmp.isValid(),false);
  
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  tmp.setPercentage("C",5,0);
  tmp.setPercentage("A",7,30);
  tmp.setPercentage("B",7,50);
  tmp.setPercentage("C",7,20);
  TEST_EQUAL(tmp.isValid(),true);
  
  tmp.setPercentage("A",5,91);
  TEST_EQUAL(tmp.isValid(),false);
  tmp.setPercentage("B",5,9);
  TEST_EQUAL(tmp.isValid(),true);

  tmp.setPercentage("A",7,29);
  TEST_EQUAL(tmp.isValid(),false);
  tmp.setPercentage("B",7,51);
  TEST_EQUAL(tmp.isValid(),true);

  tmp.setPercentage("C",5,1);
  TEST_EQUAL(tmp.isValid(),false);
  tmp.setPercentage("A",5,90);
  TEST_EQUAL(tmp.isValid(),true);

  tmp.setPercentage("C",7,19);
  TEST_EQUAL(tmp.isValid(),false);
  tmp.setPercentage("A",7,30);
  TEST_EQUAL(tmp.isValid(),true);
END_SECTION

START_SECTION((Gradient(const Gradient& source)))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  
  Gradient tmp2(tmp);
  TEST_EQUAL(tmp2.getPercentages().size(),2);
  TEST_EQUAL(tmp2.getPercentages()[0].size(),1);
  TEST_EQUAL(tmp2.getPercentages()[1].size(),1);
  TEST_EQUAL(tmp2.getPercentage("A",5),90);
  TEST_EQUAL(tmp2.getPercentage("B",5),10);
END_SECTION

START_SECTION((Gradient& operator = (const Gradient& source)))
  Gradient tmp;
  tmp.addTimepoint(5);
  tmp.addEluent("A");
  tmp.addEluent("B");
  tmp.setPercentage("A",5,90);
  tmp.setPercentage("B",5,10);
  
  Gradient tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getPercentages().size(),2);
  TEST_EQUAL(tmp2.getPercentages()[0].size(),1);
  TEST_EQUAL(tmp2.getPercentages()[1].size(),1);
  TEST_EQUAL(tmp2.getPercentage("A",5),90);
  TEST_EQUAL(tmp2.getPercentage("B",5),10);

  tmp2 = Gradient();
  TEST_EQUAL(tmp2.getPercentages().size(),0);
	TEST_EQUAL(tmp2.getTimepoints().size(),0);
	TEST_EQUAL(tmp2.getEluents().size(),0);
END_SECTION

START_SECTION((bool operator == (const Gradient& source) const))
  Gradient base;
  base.addTimepoint(5);
  base.addEluent("A");
  base.addEluent("B");
  base.setPercentage("A",5,90);
  base.setPercentage("B",5,10);
	
	Gradient edit(base);
	TEST_EQUAL(edit==base,true);
	
	edit.addEluent("C");
	TEST_EQUAL(edit==base,false);
	
	edit = base;
	edit.addTimepoint(7);
	TEST_EQUAL(edit==base,false);

	edit = base;
	edit.setPercentage("B",5,11);
	TEST_EQUAL(edit==base,false);
END_SECTION

START_SECTION((bool operator != (const Gradient& source) const))
  Gradient base;
  base.addTimepoint(5);
  base.addEluent("A");
  base.addEluent("B");
  base.setPercentage("A",5,90);
  base.setPercentage("B",5,10);
	
	Gradient edit(base);
	TEST_EQUAL(edit!=base,false);
	
	edit.addEluent("C");
	TEST_EQUAL(edit!=base,true);
	
	edit = base;
	edit.addTimepoint(7);
	TEST_EQUAL(edit!=base,true);

	edit = base;
	edit.setPercentage("B",5,11);
	TEST_EQUAL(edit!=base,true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



