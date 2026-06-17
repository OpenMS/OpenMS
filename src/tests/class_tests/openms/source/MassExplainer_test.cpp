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

START_SECTION((SignedSize queryMultimer(const Int net_charge, const double m1, const double m2, const Int n1, const Int n2, const double tolerance, const float thresh_log_p, std::vector<CompomerIterator>& hits) const))
{
  // Use default MassExplainer (H+, Na+, NH4+, K+ adducts) which produces
  // net_charge=0 compomers with adducts on both LEFT and RIGHT sides.
  MassExplainer me;
  me.compute();

  // Simulate a monomer-dimer pair with matching adducts.
  // Feature1 = [M+K]+ at some m/z (monomer, n1=1)
  // Feature2 = [2M+H]+ at some m/z (dimer, n2=2)
  // M = 99.0
  // K mass ~ 38.963158, H mass ~ 1.00728
  // m1 = mz1 * |q1| = M + K_mass = 99.0 + 38.963158 = 137.963158
  // m2 = mz2 * |q2| = 2*M + H_mass = 2*99.0 + 1.00728 = 199.00728
  // Cross-mult: n2*m1 - n1*m2 = 2*137.963158 - 1*199.00728 = 76.919036
  // Expected from compomer (K1 left, H1 right): n2*K_mass - n1*H_mass = 2*38.963158 - 1*1.00728 = 76.919036
  double K_mass = 38.963158;
  double H_mass = 1.00728;
  double M = 99.0;
  double m1_val = M + K_mass;       // 137.963158
  double m2_val = 2.0 * M + H_mass; // 199.00728

  std::vector<MassExplainer::CompomerIterator> multimer_hits;
  SignedSize n_hits = me.queryMultimer(0, m1_val, m2_val, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits > 0, true);

  // Wrong multiplier combo should find no match
  multimer_hits.clear();
  n_hits = me.queryMultimer(0, m1_val, m2_val, 1, 3, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits, 0);

  // Test with non-zero net_charge: cross-charge cross-adduct scenario.
  // Feature1 = [M+K]+ (q1=1, monomer), Feature2 = [2M+2H]2+ (q2=2, dimer)
  // m1 = M + K_mass = 99.0 + 38.963158 = 137.963158
  // m2 = 2*M + 2*H_mass = 198.0 + 2*1.00728 = 200.01456
  // net_charge = q2 - q1 = 1
  // observed = 2*137.963158 - 200.01456 = 75.91176
  // Compomer K1_LEFT + H2_RIGHT has net_charge = -1+2 = +1
  //   left_mass = K_mass, right_mass = 2*H_mass
  //   expected = 2*K_mass - 1*(2*H_mass) = 77.926316 - 2.01456 = 75.91176
  double m1_cross = M + K_mass;                        // 137.963158
  double m2_cross = 2.0 * M + 2.0 * H_mass;           // 200.01456
  multimer_hits.clear();
  n_hits = me.queryMultimer(1, m1_cross, m2_cross, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits > 0, true);

  // Wrong net_charge should find no match for these masses
  multimer_hits.clear();
  n_hits = me.queryMultimer(-1, m1_cross, m2_cross, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits, 0);
}
END_SECTION

START_SECTION((void compute(bool include_identity)))
{
  // Without identity compomers (default), no mass=0 net_charge=0 compomers exist
  MassExplainer me_no_id;
  me_no_id.compute();  // include_identity defaults to false

  MassExplainer::CompomerIterator s, e;
  // mass=0, net_charge=0 should have no hits without identity
  SignedSize hits_no_id = me_no_id.query(0, 0.0, 0.5, -100000, s, e);
  TEST_EQUAL(hits_no_id, 0);

  // With identity compomers, should find compomers at mass=0, net_charge=0
  MassExplainer me_id;
  me_id.compute(true);  // include_identity = true

  SignedSize hits_id = me_id.query(0, 0.0, 0.5, -100000, s, e);
  TEST_EQUAL(hits_id > 0, true);

  // Each identity compomer should have getSideMass(LEFT) == getSideMass(RIGHT)
  for (auto it = s; it != e; ++it)
  {
    TEST_REAL_SIMILAR(it->getSideMass(Compomer::LEFT), it->getSideMass(Compomer::RIGHT));
    TEST_EQUAL(it->getNetCharge(), 0);
  }

  // Same-adduct queryMultimer should now work
  // H+ identity: expected = 2*H_mass - 1*H_mass = H_mass
  // Use default adducts (H+, Na+, NH4+, K+)
  std::vector<MassExplainer::CompomerIterator> multimer_hits;
  // Query with approximate H+ mass (proton mass ~ 1.007)
  double H_mass = 1.00728;
  double m1 = 100.0 + H_mass;       // [M+H]+ with M=100
  double m2 = 200.0 + H_mass;       // [2M+H]+ with M=100
  SignedSize n_hits = me_id.queryMultimer(0, m1, m2, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits > 0, true);

  // Without identity, same query should fail
  multimer_hits.clear();
  n_hits = me_no_id.queryMultimer(0, m1, m2, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits, 0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
