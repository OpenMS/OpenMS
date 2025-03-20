// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(FineIsotopePatternGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FineIsotopePatternGenerator* ptr = nullptr;
FineIsotopePatternGenerator* nullPointer = nullptr;
START_SECTION((FineIsotopePatternGenerator()))
  ptr = new FineIsotopePatternGenerator();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~FineIsotopePatternGenerator()))
  delete ptr;
END_SECTION

START_SECTION(( IsotopeDistribution run(const EmpiricalFormula&) const ))
{
  EmpiricalFormula ef ("C6H12O6");

  // simple way of getting an IsotopeDistribution
  IsotopeDistribution test_id = ef.getIsotopeDistribution(FineIsotopePatternGenerator(0.01, false, false));
  TEST_EQUAL(test_id.size(), 3)

  // simple way of getting an IsotopeDistribution using absolute tol
  test_id = ef.getIsotopeDistribution(FineIsotopePatternGenerator(0.01, false, true));
  TEST_EQUAL(test_id.size(), 3)

  // simple way of getting an IsotopeDistribution using total probability
  test_id = ef.getIsotopeDistribution(FineIsotopePatternGenerator(0.01, true, false));
  TEST_EQUAL(test_id.size(), 3)

  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922633) // 0.922119)

    TEST_REAL_SIMILAR(id[2].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.0113774 )
  }

  {
    const double threshold = 1e-5;
    FineIsotopePatternGenerator gen(threshold, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 14)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922633)

    TEST_REAL_SIMILAR(id[4].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[4].getIntensity(), 0.0113774 )

    TEST_REAL_SIMILAR(id[13].getMZ(), 184.07434277234)
    TEST_REAL_SIMILAR(id[13].getIntensity(), 2.02975552383577e-05)
  }

  {
    FineIsotopePatternGenerator gen(1e-12, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 104)

    gen.setThreshold(1e-25);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 634)

    gen.setThreshold(1e-50);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 1883)

    gen.setThreshold(1e-100);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 2548)

    gen.setThreshold(0.0);
    TEST_EQUAL(gen.run(EmpiricalFormula(ef)).size(), 2548)
  }

  // For a C100 molecule
  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setThreshold(1e-2);
    IsotopeDistribution id = gen.run(EmpiricalFormula("C100"));
    TEST_EQUAL(id.size(), 6)

    // for (auto i : id.getContainer())
    //   std::cout << i << std::endl;

    gen.setThreshold(1e-5);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 9)

    gen.setThreshold(1e-10);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 14)

    gen.setThreshold(1e-20);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 21)

    gen.setThreshold(1e-40);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 34)

    gen.setThreshold(1e-60);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 46)

    gen.setThreshold(1e-100);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 65)

    gen.setThreshold(1e-150);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 86)

    gen.setThreshold(1e-196);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 100)

    gen.setThreshold(1e-198);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 101)

    gen.setThreshold(1e-250);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 101)

    gen.setThreshold(0.0);
    TEST_EQUAL(gen.run(EmpiricalFormula("C100")).size(), 101)

    TEST_REAL_SIMILAR(gen.run(EmpiricalFormula("C100"))[100].getIntensity(), 8.67e-198) // note: Intensity is only float, so nothing beyond 1e-38
    TEST_REAL_SIMILAR(gen.run(EmpiricalFormula("C100"))[100].getMZ(), 1300.3355000000001)
  }

  // {
  //   std::string formula = "C100H202"; // add 202 hydrogen
  //   FineIsotopePatternGenerator gen(0.99, true, true);
  //   IsotopeDistribution id = gen.run(EmpiricalFormula(formula));
  //   TEST_EQUAL(id.size(), 9)

  //   for (auto i : id.getContainer())
  //     std::cout << i << std::endl;
  // }

  {
    std::string formula = "C100H202"; // add 202 hydrogen
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setThreshold(1e-2);
    IsotopeDistribution id = gen.run(EmpiricalFormula(formula));
    TEST_EQUAL(id.size(), 9)

    gen.setThreshold(1e-5);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 21)
    id = gen.run(EmpiricalFormula(formula));

    gen.setThreshold(1e-10);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 50)

    gen.setThreshold(1e-20);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 131)

    gen.setThreshold(1e-40);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 368)

    gen.setThreshold(1e-60);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 677)

    gen.setThreshold(1e-100);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 1474)

    gen.setThreshold(1e-150);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 2743)

    gen.setThreshold(1e-250);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 5726)

    gen.setThreshold(1e-320);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 7687)

    gen.setThreshold(0.0);
    TEST_EQUAL(gen.run(EmpiricalFormula(formula)).size(), 101* 203)
  }

  // Also test a molecule with 2048 atoms (a value that does not fit into the
  // lookup table any more, it should still work).
  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setThreshold(1e-2);
    IsotopeDistribution id = gen.run(EmpiricalFormula("C2048"));
    TEST_EQUAL(id.size(), 28)

    gen.setThreshold(1e-5);
    TEST_EQUAL(gen.run(EmpiricalFormula("C2048")).size(), 44)
  }
}
END_SECTION

START_SECTION(( [EXTRA CH]IsotopeDistribution run(const EmpiricalFormula&) const ))
{
  EmpiricalFormula ef ("C6H12O6");

  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922633) // 0.922119)

    TEST_REAL_SIMILAR(id[2].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.0113774 )
  }

  ef.setCharge(2);
  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.922633) // 0.922119)

    TEST_REAL_SIMILAR(id[2].getMZ(), 182.068 ) 
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.0113774 )
  }
}
END_SECTION

START_SECTION(( [EXTRA]IsotopeDistribution run(const EmpiricalFormula&) const ))
{
  {
    // human insulin
    EmpiricalFormula ef ("C520H817N139O147S8");

    FineIsotopePatternGenerator gen(0.01, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 267)

    gen.setThreshold(1e-5);
    IsotopeDistribution id2 = gen.run(ef);
    TEST_EQUAL(id2.size(), 5513)

    IsotopeDistribution id3 = ef.getIsotopeDistribution(FineIsotopePatternGenerator(0.01, false, false));
    TEST_EQUAL(id3.size(), 267)

    IsotopeDistribution id4 = ef.getIsotopeDistribution(FineIsotopePatternGenerator(1e-5, false, false));
    TEST_EQUAL(id4.size(), 5513)
  }

  {
    EmpiricalFormula ef("C222N190O110");
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setThreshold(1e-3);
    IsotopeDistribution id = gen.run(ef);

    TEST_EQUAL(id.size(), 154)

    // int idx = 0; for (const auto ele : id ) std::cout << idx++ << " : " << ele << std::endl;

    TEST_REAL_SIMILAR(id[0].getMZ(), 7084.02466902)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.0348636) // cmp with 0.0349429

    TEST_REAL_SIMILAR(id[1].getMZ(), 7085.0217039152)
    TEST_REAL_SIMILAR(id[2].getMZ(), 7085.0280238552)
    TEST_REAL_SIMILAR(id[3].getMZ(), 7085.0288861574)

    TEST_REAL_SIMILAR(id[1].getIntensity() + id[2].getIntensity() + id[3].getIntensity(), 0.109638) // cmp with 0.109888

    TEST_REAL_SIMILAR(id[4].getMZ(), 7086.0187388104)
    TEST_REAL_SIMILAR(id[9].getMZ(), 7086.0322409926)
    TEST_REAL_SIMILAR(id[4].getIntensity() +
                      id[5].getIntensity() +
                      id[6].getIntensity() +
                      id[7].getIntensity() +
                      id[8].getIntensity() +
                      id[9].getIntensity(), 0.179746) // cmp with 0.180185 -- difference of 0.24%

    TEST_REAL_SIMILAR(id[10].getMZ(), 7087.0157737056)
    TEST_REAL_SIMILAR(id[19].getMZ(), 7087.0355958278)
    TEST_REAL_SIMILAR(id[10].getIntensity() +
                      id[11].getIntensity() +
                      id[12].getIntensity() +
                      id[13].getIntensity() +
                      id[14].getIntensity() +
                      id[15].getIntensity() +
                      id[16].getIntensity() +
                      id[17].getIntensity() +
                      id[18].getIntensity() +
                      id[19].getIntensity(),
                      0.203836) // cmp with 0.204395 -- difference of 0.27%

    // Cmp with CoarseIsotopePatternGenerator:
    // container.push_back(IsotopeDistribution::MassAbundance(7084, 0.0349429));
    // container.push_back(IsotopeDistribution::MassAbundance(7085, 0.109888));
    // container.push_back(IsotopeDistribution::MassAbundance(7086, 0.180185));
    // container.push_back(IsotopeDistribution::MassAbundance(7087, 0.204395));
    // container.push_back(IsotopeDistribution::MassAbundance(7088, 0.179765));
    // container.push_back(IsotopeDistribution::MassAbundance(7089, 0.130358));
    // container.push_back(IsotopeDistribution::MassAbundance(7090, 0.0809864));
    // container.push_back(IsotopeDistribution::MassAbundance(7091, 0.0442441));
    // container.push_back(IsotopeDistribution::MassAbundance(7092, 0.0216593));
    // container.push_back(IsotopeDistribution::MassAbundance(7093, 0.00963707));
    // container.push_back(IsotopeDistribution::MassAbundance(7094, 0.0039406));
  }

  {
    // test gapped isotope distributions, e.g. bromide 79,81 (missing 80)

    EmpiricalFormula ef("CBr2");
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setThreshold(1e-3);
    IsotopeDistribution id = gen.run(ef);

    // int idx = 0; for (const auto ele : id ) std::cout << idx++ << " : " << ele << std::endl;

    TEST_REAL_SIMILAR(id[0].getMZ(), 169.8366742)
    TEST_REAL_SIMILAR(id[1].getMZ(), 170.8400292)

    IsotopeDistribution::ContainerType container;
    container.push_back(IsotopeDistribution::MassAbundance(170, 0.254198270573));
    container.push_back(IsotopeDistribution::MassAbundance(171, 0.002749339427));
    container.push_back(IsotopeDistribution::MassAbundance(172, 0.494555798854));
    container.push_back(IsotopeDistribution::MassAbundance(173, 0.005348981146));
    container.push_back(IsotopeDistribution::MassAbundance(174, 0.240545930573));
    container.push_back(IsotopeDistribution::MassAbundance(175, 0.002601679427));
    for (Size i = 0; i != id.size(); ++i)
    {
      TEST_EQUAL(round(id.getContainer()[i].getMZ()), container[i].getMZ())
      TEST_REAL_SIMILAR(id.getContainer()[i].getIntensity(), container[i].getIntensity())
    }
  }

#if 0
  // Do some stress testing of the library...
  // Stress test takes about 20 seconds
  // there is a significant drop in speed due to copying (and sorting) of data
  int sum = 0;
  for (Size k = 0; k < 2e5; k++)
  {
    EmpiricalFormula ef ("C520H817N139O147");
    FineIsotopePatternGenerator gen(1e-2, false, false);
    IsotopeDistribution id = gen.run(ef);
    sum += id.size();
  }
  TEST_EQUAL(sum, 139*2*1e5) // we use OpenMS isotopic tables, we get 139 instead of 140

  int calculated_masses = 0;
  for (Size k = 0; k < 100; k++)
  {
    // human insulin
    EmpiricalFormula ef (String("C") + (520 + k) +
                         String("H") + (817 + k) + 
                         String("N") + (139 + k) +
                         String("O") + (147 + k) +
                         String("S") + ( 8  + int(k/5)) ); // Sulfur is hard to do because of the abundant isotope 34

    std::cout << " Working on stress test " << k << " " << ef.toString() << std::endl;

    {
      FineIsotopePatternGenerator gen(0.01, false, false);
      IsotopeDistribution id = gen.run(ef);
      calculated_masses += id.size();

      gen.setThreshold(1e-5);
      id = gen.run(ef);
      calculated_masses += id.size();
    }
  }
  TEST_EQUAL(calculated_masses, 1592882)
  for (Size k = 0; k < 100; k++)
  {
    // human insulin
    EmpiricalFormula ef (String("C") + (520 + k) +
                         String("H") + (817 + k) + 
                         String("N") + (139 + k) +
                         String("O") + (147 + k) +
                         String("S") + ( 8  + int(k/5)) ); // Sulfur is hard to do because of the abundant isotope 34

    std::cout << " Working on stress test " << k << " " << ef.toString() << std::endl;

    {
      FineIsotopePatternGenerator gen(0.01, false, false);
      IsotopeDistribution id = gen.run(ef);
      calculated_masses += id.size();

      gen.setThreshold(1e-5);
      id = gen.run(ef);
      calculated_masses += id.size();
    }
  }
  TEST_EQUAL(calculated_masses, 1592882*2) // repeat the test, we should get the same result
#endif
}
END_SECTION

START_SECTION(( void setAbsolute(bool absolute) ))
{
  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    gen.setAbsolute(true);
    TEST_EQUAL(gen.getAbsolute(), true);
    gen.setAbsolute(false);
    TEST_EQUAL(gen.getAbsolute(), false);
  }
  // human insulin
  EmpiricalFormula ef ("C520H817N139O147S8");

  {
    FineIsotopePatternGenerator gen(0.01, false, false);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 267)

    gen.setAbsolute(true);
    id = gen.run(ef);
    TEST_EQUAL(id.size(), 21)

    gen.setThreshold(1e-3);
    id = gen.run(ef);
    TEST_EQUAL(id.size(), 151)
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
