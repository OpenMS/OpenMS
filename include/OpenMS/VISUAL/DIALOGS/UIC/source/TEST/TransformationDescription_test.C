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
// $Maintainer: Clemens Groepl $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <vector>

///////////////////////////

START_TEST(TransformationDescription, "$Id: TransformationDescription_test.C 6446 2009-11-20 16:21:41Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


TransformationDescription* ptr = 0;
START_SECTION((TransformationDescription()))
	ptr = new TransformationDescription;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~TransformationDescription()))
	delete ptr;
END_SECTION

START_SECTION((const String& getName() const))
	TransformationDescription td;
	TEST_STRING_EQUAL(td.getName(),"")
END_SECTION

START_SECTION((void setName(const String& name)))
	TransformationDescription td;
	td.setName("bla");
	TEST_STRING_EQUAL(td.getName(),"bla")

END_SECTION

START_SECTION((const Param& getParameters() const))
	TransformationDescription td;
	TEST_EQUAL(td.getParameters(),Param())
END_SECTION

START_SECTION((const DataValue& getParam(const String &name) const))
	TransformationDescription td;
	TEST_EXCEPTION(Exception::ElementNotFound, td.getParam("bla"))
END_SECTION

START_SECTION((void setParam(const String& name, DoubleReal value)))
  TransformationDescription td;
  td.setParam("bla",4.5);
  TEST_REAL_SIMILAR(td.getParam("bla"),4.5)
END_SECTION

START_SECTION((void setParam(const String& name, Int value)))
  TransformationDescription td;
  td.setParam("bla",17);
  TEST_EQUAL(td.getParam("bla"),17)
END_SECTION

START_SECTION((void setParam(const String& name, const String& value)))
  TransformationDescription td;
  td.setParam("bla","yummyummmmyummmmmy");
  TEST_EQUAL(td.getParam("bla"),"yummyummmmyummmmmy")
END_SECTION


START_SECTION((void setParameters(const Param& param)))
	TransformationDescription td;
	Param p;
	p.setValue("int",5);
	td.setParameters(p);
	TEST_EQUAL((Int)td.getParameters().size(),1)
	TEST_EQUAL((Int)td.getParameters().getValue("int"),5)
END_SECTION


TransformationDescription::PairVector pairs;
pairs.push_back(make_pair(1.2,5.2));
pairs.push_back(make_pair(3.2,7.3));
pairs.push_back(make_pair(2.2,6.25));

START_SECTION((const PairVector& getPairs() const))
	TransformationDescription td;
	TEST_EQUAL(td.getPairs().size(),0)
END_SECTION

START_SECTION((PairVector& getPairs()))
{
  TransformationDescription td;
  TEST_EQUAL(td.getPairs().size(),0)
  td.getPairs().push_back(make_pair(12.34,56.78));
  TEST_EQUAL(td.getPairs().size(),1)
  TEST_EQUAL(td.getPairs()[0].first,12.34); // not TEST_REAL_SIMILAR!  I want to see an alert when someone messes up float vs. double!
  TEST_EQUAL(td.getPairs()[0].second,56.78);
  TEST_NOT_EQUAL(td.getPairs().empty(),true);
  td.getPairs().clear();
  TEST_EQUAL(td.getPairs().empty(),true);
}
END_SECTION

START_SECTION((void setPairs(const PairVector& pairs)))
{
	TransformationDescription td;
	td.setPairs(pairs);
	TEST_EQUAL(td.getPairs().size(),3)

	TransformationDescription::PairVector pairs_empty;
  pairs_empty.clear();
	td.setPairs(pairs_empty);
	TEST_EQUAL(td.getPairs().size(),0)
}
END_SECTION

START_SECTION((TransformationDescription(const TransformationDescription& rhs)))
{
	TransformationDescription td;
	td.setName("dummy");
	td.setParam("int",5);
	td.setPairs(pairs);

	TEST_EQUAL(td.getName()==td.getName(),true)
	TEST_EQUAL(td.getParameters()==td.getParameters(),true)
	TEST_EQUAL(td.getPairs().size(),3)
}
END_SECTION

START_SECTION((TransformationDescription& operator = (const TransformationDescription& rhs)))
{
	TransformationDescription td;
	td.setName("dummy");
	td.setParam("int",5);
	td.setPairs(pairs);
	TransformationDescription td2;
	td2 = td;

 	TEST_STRING_EQUAL(td2.getName(),td.getName());
	TEST_EQUAL(td2.getParameters()==td.getParameters(),true);
	TEST_EQUAL(td2.getPairs()==td.getPairs(),true);
}
END_SECTION

START_SECTION((void clear()))
{
	TransformationDescription td;

	td.setName("linear");
	td.setParam("slope",2.0);
	td.setParam("intercept",47.12);
	td.setPairs(pairs);

	DoubleReal value = 5.0;
	td.apply(value);
	TEST_REAL_SIMILAR(value,57.12);

 	TEST_STRING_EQUAL(td.getName(),"linear");
	TEST_EQUAL((DoubleReal)td.getParameters().getValue("slope"),2.0);
	TEST_EQUAL((DoubleReal)td.getParameters().getValue("intercept"),47.12);
	TEST_EQUAL(td.getPairs()==pairs,true);
	TEST_EQUAL(td.getPairs().size(),3);

	td.clear();

 	TEST_STRING_EQUAL(td.getName(),"");
	TEST_EQUAL(td.getParameters().empty(),true);
	TEST_EQUAL(td.getPairs()==pairs,false);
	TEST_EQUAL(td.getPairs().size(),0);
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));
}
END_SECTION

START_SECTION((void apply(DoubleReal &value) const))
{
	DoubleReal value = 5.0;
	TransformationDescription td;

	//test missing name and parameters
 	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));

	td.setName("bla");
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));

	//test with identity
	td.setName("none");
	td.apply(value);
	TEST_REAL_SIMILAR(value,5.0);

	//test for missing parameter
	td.setName("linear");
	td.setParam("slope",1.0);
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));

	//real test (linear, identity)
	td.setParam("intercept",0.0);
	TEST_REAL_SIMILAR(value,5.0);

	//real test (linear, no identity)
	td.setParam("slope",2.0);
	td.setParam("intercept",47.12);
	td.apply(value);
	TEST_REAL_SIMILAR(value,57.12);

	td.clear();
	td.setName("interpolated_linear");
	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));
	td.setPairs(pairs);
	td.apply(value);

	// VALUES FROM ABOVE:
	// pairs.push_back(make_pair(1.2,5.2));
	// pairs.push_back(make_pair(2.2,6.25));
	// pairs.push_back(make_pair(3.2,7.3));

	value = 0.2;
	td.apply(value);
	TEST_REAL_SIMILAR(value,4.15);

	value = 0.7;
	td.apply(value);
	TEST_REAL_SIMILAR(value,4.675);

	value = 1.2;
	td.apply(value);
	TEST_REAL_SIMILAR(value,5.2);

	value = 1.45;
	td.apply(value);
	TEST_REAL_SIMILAR(value,5.4625);

	value = 1.7;
	td.apply(value);
	TEST_REAL_SIMILAR(value,5.725);

	value = 2.2;
	td.apply(value);
	TEST_REAL_SIMILAR(value,6.25);

	value = 2.45;
	td.apply(value);
	TEST_REAL_SIMILAR(value, 6.5125);

	value = 2.7;
	td.apply(value);
	TEST_REAL_SIMILAR(value,6.775);

	value = 3.2;
	td.apply(value);
	TEST_REAL_SIMILAR(value,7.3);

	value = 4.2;
	td.apply(value);
	TEST_REAL_SIMILAR(value,8.35);

	//--------------------------

	td.clear();
	td.setName("b_spline");
	td.setParam("num_breakpoints", 4);

	TEST_EXCEPTION(Exception::IllegalArgument,td.apply(value));

	pairs.clear();
  pairs.push_back(make_pair(1.2,5.2));
  pairs.push_back(make_pair(3.2,7.3));
  pairs.push_back(make_pair(2.2,6.25));
  pairs.push_back(make_pair(2.2,3.1));
  pairs.push_back(make_pair(2.2,7.25));
  pairs.push_back(make_pair(3.0,8.5));
  pairs.push_back(make_pair(3.1,4.7));
  pairs.push_back(make_pair(1.7,6.0));
  pairs.push_back(make_pair(2.9,4.7));
  pairs.push_back(make_pair(4.2,5.0));
  pairs.push_back(make_pair(3.7,-2.4));

  td.setPairs(pairs);


#if 0
    // Since the numbers in this test were verified by manual (in fact, visual) inspection...
    // Here is a recipe how this was done:
    //
    // To grep for output, "pairs:" and "spline:", you might use a command line like this:
    // make TransformationDescription_test && ctest -V -R TransformationDescription_test &&  ctest -V -R TransformationDescription_test | grep pairs: > points.dat && ctest -V -R TransformationDescription_test | grep spline: > bla.dat
    // To have a look at the results using gnuplot:
    // gnuplot -  # (to start gnuplot)
    // plot './bla.dat' u 5:6, './points.dat' u 5:6

    for ( UInt i = 0; i < pairs.size(); ++i )
    {
      STATUS("pairs: " << pairs[i].first << " " << pairs[i].second);
    }

    for ( Int i = -10; i <= 60; i+=5 )
    {
      DoubleReal value = i;
      value /= 10;
      DoubleReal image = value;
      td.apply(image);
      STATUS("spline: " << value << " " << image);
    }
#endif

    DoubleReal sample_values[] =
        { -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6 };
    DoubleReal sample_images[] =
        { -14.3519415123977, -9.91557518507088, -5.4792088577441, -1.04284253041731, 3.39352379690948, 6.4561466812738, 5.4858954730427,
            6.14659387774751, 6.77299727168147, 0.646024122587505, -1.13062259235381, 18.3842099268184, 40.7826815802615, 63.1811532337045,
            85.5796248871476 };
    for ( Size i = 0; i < sizeof(sample_values)/sizeof(*sample_values); ++i)
    {
      DoubleReal x = sample_values[i];
      td.apply(x);
      TEST_REAL_SIMILAR(x,sample_images[i]);
    }

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
