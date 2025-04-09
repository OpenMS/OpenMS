// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Compomer.h>

using namespace OpenMS;
using namespace std;

START_TEST(MassExplainer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassExplainer* ptr = nullptr;
MassExplainer* nullPointer = nullptr;
START_SECTION(MassExplainer())
{
  ptr = new MassExplainer();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MassExplainer())
{
	delete ptr;
}
END_SECTION

START_SECTION((MassExplainer(AdductsType adduct_base)))
{
	Adduct a(2, 1, 123.12, "Na", -0.5,0);
	MassExplainer::AdductsType va;
	va.push_back(a);  
	MassExplainer me(va);
	TEST_EQUAL(me.getAdductBase().size(), 1);
}
END_SECTION

START_SECTION((MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp)))
{
	MassExplainer me(5,10,2,-10.3);
	TEST_EQUAL(me.getAdductBase().size(), 4);
}
END_SECTION

START_SECTION((MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals)))
{
	MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me(va,5,10,2,-10.3,0);
	TEST_EQUAL(me.getAdductBase().size(), 2);
}
END_SECTION

START_SECTION((MassExplainer& operator=(const MassExplainer &rhs)))
{
  MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me(va,5,10,2,-10.3,0);
	MassExplainer me2;
	me2 = me;
	TEST_EQUAL(me.getAdductBase().size(), 2);

}
END_SECTION

START_SECTION((void setAdductBase(AdductsType adduct_base)))
{
	MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me;
	me.setAdductBase(va);
	TEST_EQUAL(me.getAdductBase().size(),2);
}
END_SECTION

START_SECTION((AdductsType getAdductBase() const))
{
  NOT_TESTABLE; // tested above
}
END_SECTION

START_SECTION((const Compomer& getCompomerById(Size id) const))
{
  MassExplainer me;
	me.compute();
	TEST_EQUAL(me.getCompomerById(0).getID(), 0)
}
END_SECTION

START_SECTION((void compute()))
{
  NOT_TESTABLE; // tested above
}
END_SECTION

START_SECTION((SignedSize query(const Int net_charge, const float mass_to_explain, const float mass_delta, const float thresh_log_p, std::vector< Compomer >::const_iterator &firstExplanation, std::vector< Compomer >::const_iterator &lastExplanation) const))
{
  MassExplainer me;
	me.compute();
	
	
	MassExplainer::CompomerIterator s,e;
	SignedSize hits = me.query(2, 45.0, 13.0, -100000,s,e);
	
	std::cout << "hits: " << hits << std::endl;
	for (; s!=e; ++s)
	{
			std::cout << *s << std::endl;
	}
	
	TEST_EQUAL(hits, 5);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
