// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

///////////////////////////

// This one is going to be tested.
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

START_TEST(CoarseIsotopePatternGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

IsotopeDistribution* nullPointer = nullptr;

START_SECTION(CoarseIsotopePatternGenerator())
{
    CoarseIsotopePatternGenerator* ptr = new CoarseIsotopePatternGenerator();
    Size max_isotope = ptr->getMaxIsotope();
    TEST_EQUAL(max_isotope, 0)
    TEST_EQUAL(ptr->getRoundMasses(), false)
    TEST_NOT_EQUAL(ptr, nullPointer)
    delete ptr;
}
END_SECTION

START_SECTION(CoarseIsotopePatternGenerator(Size max_isotope))
{
    CoarseIsotopePatternGenerator* ptr = new CoarseIsotopePatternGenerator(117);
    Size max_isotope = ptr->getMaxIsotope();
    TEST_EQUAL(max_isotope, 117)
    TEST_EQUAL(ptr->getRoundMasses(), false)
    TEST_NOT_EQUAL(ptr, nullPointer)
    delete ptr;
}
END_SECTION

START_SECTION(CoarseIsotopePatternGenerator(Size max_isotope, bool calc_mass))
{
    CoarseIsotopePatternGenerator* ptr = new CoarseIsotopePatternGenerator(117, true);
    Size max_isotope = ptr->getMaxIsotope();
    TEST_EQUAL(max_isotope, 117)
    TEST_EQUAL(ptr->getRoundMasses(), true)
    TEST_NOT_EQUAL(ptr, nullPointer)
    delete ptr;
}
END_SECTION

CoarseIsotopePatternGenerator* solver = new CoarseIsotopePatternGenerator();

START_SECTION(~CoarseIsotopePatternGenerator())
    CoarseIsotopePatternGenerator* ptr = new CoarseIsotopePatternGenerator(117);
    delete ptr;
END_SECTION

START_SECTION(void setRoundMasses(bool round_masses))
{
    CoarseIsotopePatternGenerator solver2 = CoarseIsotopePatternGenerator();
    TEST_EQUAL(solver2.getRoundMasses(), false)
    solver2.setRoundMasses(true);
    TEST_EQUAL(solver2.getRoundMasses(), true)
}
END_SECTION

START_SECTION(bool getRoundMasses() const)
    NOT_TESTABLE
END_SECTION

START_SECTION(void setMaxIsotope(Size max_isotope))
{
    IsotopeDistribution iso = solver->estimateFromPeptideWeight(1234.2);
    TEST_EQUAL(solver->getMaxIsotope(), 0)
    TEST_EQUAL(iso.getContainer().size(), 317)
    solver->setMaxIsotope(117);
    TEST_EQUAL(solver->getMaxIsotope(), 117)
}
END_SECTION

START_SECTION(Size getMaxIsotope() const)
    NOT_TESTABLE
END_SECTION

START_SECTION(IsotopeDistribution convolve_(const CoarseIsotopePatternGenerator& isotope_distribution) const)
{
    IsotopeDistribution iso1, iso2;
    solver->setMaxIsotope(1);
    IsotopeDistribution::ContainerType result = solver->convolve(iso1.getContainer(), iso2.getContainer());
    TEST_EQUAL(result.size(), 1)
    TEST_EQUAL(result[0].getMZ(), 0)
    TEST_EQUAL(result[0].getIntensity(), 1)
}
END_SECTION

START_SECTION(( [EXTRA CH]IsotopeDistribution run(const EmpiricalFormula&) const ))
{
  EmpiricalFormula ef ("C6H12O6");

  {
    CoarseIsotopePatternGenerator gen(3);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.923456)

    TEST_REAL_SIMILAR(id[2].getMZ(), 182.0701)
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.013232)
  }

  // TODO: is that a good idea?
  ef.setCharge(2);
  {
    CoarseIsotopePatternGenerator gen(3);
    IsotopeDistribution id = gen.run(ef);
    TEST_EQUAL(id.size(), 3)

    // TEST_REAL_SIMILAR(id[0].getMZ(), 180.063)
    TEST_REAL_SIMILAR(id[0].getMZ(), 182.077943)
    TEST_REAL_SIMILAR(id[0].getIntensity(), 0.923456)

    TEST_REAL_SIMILAR(id[2].getMZ(), 184.0846529)
    TEST_REAL_SIMILAR(id[2].getIntensity(), 0.013232)
  }
}
END_SECTION

START_SECTION(CoarseIsotopePatternGenerator& convolvePow_(Size factor))
{
  // IsotopeDistribution iso = EmpiricalFormula("C60H97N15O19").getIsotopeDistribution(CoarseIsotopePatternGenerator());
  // for(auto elem : iso.getContainer())
  // {
    // std::cout << elem.getMZ() << " " << elem.getIntensity() << std::endl;
  // }

    EmpiricalFormula ef("C222N190O110");
    IsotopeDistribution id = ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(11));
    IsotopeDistribution::ContainerType container;
    container.push_back(IsotopeDistribution::MassAbundance(7084, 0.0349429));
    container.push_back(IsotopeDistribution::MassAbundance(7085, 0.109888));
    container.push_back(IsotopeDistribution::MassAbundance(7086, 0.180185));
    container.push_back(IsotopeDistribution::MassAbundance(7087, 0.204395));
    container.push_back(IsotopeDistribution::MassAbundance(7088, 0.179765));
    container.push_back(IsotopeDistribution::MassAbundance(7089, 0.130358));
    container.push_back(IsotopeDistribution::MassAbundance(7090, 0.0809864));
    container.push_back(IsotopeDistribution::MassAbundance(7091, 0.0442441));
    container.push_back(IsotopeDistribution::MassAbundance(7092, 0.0216593));
    container.push_back(IsotopeDistribution::MassAbundance(7093, 0.00963707));
    container.push_back(IsotopeDistribution::MassAbundance(7094, 0.0039406));

    for (Size i = 0; i != id.size(); ++i)
    {
        TEST_EQUAL(round(id.getContainer()[i].getMZ()), container[i].getMZ())
        TEST_REAL_SIMILAR(id.getContainer()[i].getIntensity(), container[i].getIntensity())
    }

  // test gapped isotope distributions, e.g. bromide 79,81 (missing 80)
  {
    EmpiricalFormula ef("Br2");
    IsotopeDistribution id = ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(5));
    container.clear();
    // the expected results as pairs of
    // [nominal mass, probability]
    // derived via convolution of elemental probabilities; the sum of all probabilities is 1
    // For Br2, this is simply the product of Bromine x Bromine, which
    // has a light isotope (79 Da, ~50% probability) and a heavy isotope (81 Da, ~50% probability)
    container.push_back(IsotopeDistribution::MassAbundance(158, 0.2569476));  // 79+79, ~ 0.5 * 0.5
    container.push_back(IsotopeDistribution::MassAbundance(159, 0.0));        // this mass cannot be explained by two Br atoms
    container.push_back(IsotopeDistribution::MassAbundance(160, 0.49990478)); // 79+81 (or 81+79), ~ 0.5 * 0.5 + 0.5 * 0.5
    container.push_back(IsotopeDistribution::MassAbundance(161, 0.0));        // same as mass 159
    container.push_back(IsotopeDistribution::MassAbundance(162, 0.24314761)); // 81+81, ~ 0.5 * 0.5
    for (Size i = 0; i != id.size(); ++i)
    {
      TEST_EQUAL(round(id.getContainer()[i].getMZ()), container[i].getMZ())
      TEST_REAL_SIMILAR(id.getContainer()[i].getIntensity(), container[i].getIntensity())
    }
  }
  {
    // testing a formula which has more than one element (here: C and Br), since the internal computation is different
    // The convolution is similar to the one above, but add another convolution step with Carbon (hence the lightest mass is 12 Da heavier)
    EmpiricalFormula ef("CBr2");
    IsotopeDistribution id = ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(7));
    container.clear();
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
}
END_SECTION

START_SECTION(IsotopeDistribution estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P))
{
    // We are testing that the parameterized version matches the hardcoded version.
    IsotopeDistribution iso;
    IsotopeDistribution iso2;
    solver->setMaxIsotope(3);
    iso = solver->estimateFromWeightAndComp(1000.0, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0.0);
    iso2 = solver->estimateFromPeptideWeight(1000.0);
    TEST_EQUAL(iso.begin()->getIntensity(),iso2.begin()->getIntensity());
    TEST_EQUAL(iso.begin()->getMZ(),iso2.begin()->getMZ());
}
END_SECTION

START_SECTION(IsotopeDitribution CoarseIsotopePatternGenerator::estimateFromPeptideWeight(double average_weight))
{
    // hard to test as this is an rough estimate
    IsotopeDistribution iso;
    solver->setMaxIsotope(3);
    iso = solver->estimateFromPeptideWeight(100.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.949735)
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 100.170);

    iso = solver->estimateFromPeptideWeight(1000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.586906)
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 999.714);

    iso = solver->estimateFromPeptideWeight(10000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.046495)
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 9994.041);

    solver->setRoundMasses(true);
    iso = solver->estimateFromPeptideWeight(100.0);
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 100);

    iso = solver->estimateFromPeptideWeight(1000.0);
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 1000);

    iso = solver->estimateFromPeptideWeight(10000.0);
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 9994);

    solver->setRoundMasses(false);
}
END_SECTION

START_SECTION(IsotopeDitribution CoarseIsotopePatternGenerator::approximateFromPeptideWeight(double mass, int num_peaks))
{
  std::vector<float> masses_to_test = {20, 300, 1000, 2500};
  for (auto mass = masses_to_test.begin(); mass != masses_to_test.end(); ++mass)
  {
    IsotopeDistribution approximation = CoarseIsotopePatternGenerator::approximateFromPeptideWeight(*mass);
    solver->setMaxIsotope(approximation.size());
    IsotopeDistribution coarse_truth = solver->estimateFromPeptideWeight(*mass);

    // compute KL divergence (Sum over all x: P(x) * log(P(x) / Q(x)), where P is a distribution and Q its approximation
    double KL = 0;
    double sum = 0.0;
    for (UInt peak = 0; peak != approximation.size() && peak != coarse_truth.size(); ++peak) // coarse_truth.size() is 18, although approximation.size() = 20 for mass = 20
    {
      double Px = coarse_truth[peak].getIntensity();
      double Qx = approximation[peak].getIntensity();
      sum += Qx;
      if (Px != 0) // KL is 0 for Px = 0
      {
        KL += Px * log(Px / Qx);
      }
    }

    TEST_REAL_SIMILAR(sum, 1.0);
    TEST_EQUAL(0 < KL && KL < 0.05, true);
  }
}
END_SECTION

START_SECTION(IsotopeDitribution CoarseIsotopePatternGenerator::approximateIntensities(double mass, int num_peaks))
{
  std::vector<Int> masses_to_test = {20, 300, 1000, 2500};
  for (auto mass = masses_to_test.begin(); mass != masses_to_test.end(); ++mass)
  {
    std::vector<double> approximation = CoarseIsotopePatternGenerator::approximateIntensities(*mass);
    solver->setMaxIsotope(approximation.size());
    IsotopeDistribution coarse_truth = solver->estimateFromPeptideWeight(*mass);

    // compute KL divergence (Sum over all x: P(x) * log(P(x) / Q(x)), where P is a distribution and Q its approximation
    double KL = 0;
    double sum = 0.0;
    for (UInt peak = 0; peak != approximation.size() && peak != coarse_truth.size(); ++peak)// coarse_truth.size() is 18, although approximation.size() = 20 for mass = 20
    {
      double Px = coarse_truth[peak].getIntensity();
      double Qx = approximation[peak];
      sum += Qx;
      if (Px != 0)// KL is 0 for Px = 0
      {
        KL += Px * log(Px / Qx);
      }
    }

    TEST_REAL_SIMILAR(sum, 1.0);
    TEST_EQUAL(0 < KL && KL < 0.05, true);
  }
}
END_SECTION

START_SECTION(IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::vector<UInt>& precursor_isotopes))
{
    IsotopeDistribution iso;
    IsotopeDistribution iso2;
    std::set<UInt> precursor_isotopes;
    solver->setMaxIsotope(0);
    // We're isolating the M+2 precursor isotopes
    precursor_isotopes.insert(2);
    // These are regression tests, but the results also follow an expected pattern.

    // With 0 sulfurs, it should be somewhat unlikely for the fragment to be M+2.
    iso = solver->estimateForFragmentFromPeptideWeightAndS(200.0, 0, 100.0, 0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.355445559123552)

    // At such a small size, the regular averagine method should also result in 0 sulfurs.
    // The approximate EmpiricalFormulas should be the same, and therefore so should
    // their isotopic distributions.
    iso2 = solver->estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
    iso2.renormalize();

    IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso2.begin());
    for (; it1 != iso.end(); ++it1, ++it2)
    {
        TEST_EQUAL(it1->getMZ(), it2->getMZ())
        TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
    }

    // With the only sulfur being in the fragment, it's much more likely that the fragment
    // is M+2.
    iso = solver->estimateForFragmentFromPeptideWeightAndS(200.0, 1, 100.0, 1, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.900804974056174)
    // Both sulfurs are in the fragment, so it's even more likely for the fragment to be M+2
    iso = solver->estimateForFragmentFromPeptideWeightAndS(200.0, 2, 100.0, 2, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.947862830751023)
    // All 3 sulfurs are in the fragment
    iso = solver->estimateForFragmentFromPeptideWeightAndS(200.0, 3, 100.0, 3, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.969454586761089)
    // Any more sulfurs while keeping the masses constant would violate the function preconditions
}
END_SECTION

START_SECTION(IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromPeptideWeightAndS(double average_weight_precursor, UInt S))
{
    IsotopeDistribution iso;
    IsotopeDistribution iso2;
    // These are regression tests, but the results also follow an expected pattern.

    // With 0 sulfurs, it should be very unlikely for this tiny peptide to be M+2.
    solver->setMaxIsotope(3);
    iso = solver->estimateFromPeptideWeightAndS(100.0, 0);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.00290370998965918)

    // At such a small size, the regular averagine method should also result in 0 sulfurs.
    // The approximate EmpiricalFormulas should be the same, and therefore so should
    // their isotopic distributions.
    iso2 = solver->estimateFromPeptideWeightAndS(100.0, 0);
    iso2.renormalize();

    IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso2.begin());
    for (; it1 != iso.end(); ++it1, ++it2)
    {
        TEST_EQUAL(it1->getMZ(), it2->getMZ());
        TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
    }

    // With one sulfur, it's more likely that the precursor is M+2 compared to having 0 sulfurs.
    iso = solver->estimateFromPeptideWeightAndS(100.0, 1);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.0439547771832361)
    // With two sulfurs, the M+2 isotope is more likely
    iso = solver->estimateFromPeptideWeightAndS(100.0, 2);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.0804989104418586)
    // With three sulfurs, the M+2 isotope is even more likely
    iso = solver->estimateFromPeptideWeightAndS(100.0, 3);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.117023432503842)
    // Any more sulfurs while keeping the masses constant would violate the function preconditions
}
END_SECTION

START_SECTION(IsotopeDistribution estimateFromRNAWeight(double average_weight))
{
  // hard to test as this is an rough estimate
  IsotopeDistribution iso;
  solver->setMaxIsotope(3);
  iso = solver->estimateFromRNAWeight(100.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.958166)

  iso = solver->estimateFromRNAWeight(1000.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.668538)

  iso = solver->estimateFromRNAWeight(10000.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.080505)
}
END_SECTION

START_SECTION(IsotopeDistribution CoarseIsotopePatternGenerator::estimateFromDNAWeight(double average_weight))
{
  // hard to test as this is an rough estimate
  IsotopeDistribution iso;
  solver->setMaxIsotope(3);

  iso = solver->estimateFromDNAWeight(100.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.958166)

  iso = solver->estimateFromDNAWeight(1000.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.657083)

  iso = solver->estimateFromDNAWeight(10000.0);
  TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.075138)
}
END_SECTION

START_SECTION(IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
{
    IsotopeDistribution iso;
    solver->setMaxIsotope(0);
    std::set<UInt> precursor_isotopes;
    // We're isolating the M0 and M+1 precursor isotopes
    precursor_isotopes.insert(0);
    precursor_isotopes.insert(1);
    // These are regression tests, but the results also follow an expected pattern.

    // All the fragments from the M0 precursor will also be monoisotopic, while a fragment
    // that is half the mass of the precursor and coming from the M+1 precursor will be
    // roughly 50/50 monoisotopic/M+1.
    // For such a small peptide, the M0 precursor isotope is much more abundant than M+1.
    // Therefore, it's much more likely this fragment will be monoisotopic.
    iso = solver->estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.954654801320083)
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 100.170);

    // This peptide is large enough that the M0 and M+1 precursors are similar in abundance.
    // However, since the fragment is only 1/20th the mass of the precursor, it's
    // much more likely for the extra neutron to be on the complementary fragment.
    // Therefore, it's still much more likely this fragment will be monoisotopic.
    iso = solver->estimateForFragmentFromPeptideWeight(2000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.975984866212216)
    // The size of the precursor should have no effect on the mass of the fragment
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 100.170);

    // Same explanation as the previous example.
    iso = solver->estimateForFragmentFromPeptideWeight(20000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995783521351781)

    // Like the first example, the fragment is half the mass of the precursor so
    // the fragments from the M+1 precursor will be roughly 50/50 monoisotopic/M+1.
    // However, this time the peptide is larger than the first example, so the M0
    // and M+1 precursors are also roughly 50/50 in abundance. Therefore, the
    // probability of the M0 fragment should be in the 75% range.
    //                  i.e. (100% * 50%) + (50% * 50%) = 75%
    // M0 frags due to M0 precursor^             ^M0 frags due to M+1 precursor
    iso = solver->estimateForFragmentFromPeptideWeight(2000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.741290977639283)
    TEST_REAL_SIMILAR(iso.begin()->getMZ(), 999.714);

    // Same explanation as the second example, except now the M+1 precursor is
    // more abundant than the M0 precursor. But, the fragment is so small that
    // it's still more likely for the fragment to be monoisotopic.
    iso = solver->estimateForFragmentFromPeptideWeight(20000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.95467154987681)

    // Same explanation as above.
    iso = solver->estimateForFragmentFromPeptideWeight(20000.0, 10000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.542260764523188)

    // If the fragment is identical to the precursor, then the distribution
    // should be the same as if it was just a precursor that wasn't isolated.
    iso = solver->estimateForFragmentFromPeptideWeight(200.0, 200.0, precursor_isotopes);
    IsotopeDistribution iso_precursor;
    iso_precursor = solver->estimateFromPeptideWeight(200.0);
    IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

    for (; it1 != iso.end(); ++it1, ++it2)
    {
        TEST_EQUAL(it1->getMZ(), it2->getMZ())
        TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
    }

    solver->setRoundMasses(true);

    // Rounded masses
    iso = solver->estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
    TEST_EQUAL(iso.begin()->getMZ(), 100);

    solver->setRoundMasses(false);
}
END_SECTION

START_SECTION(IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
{
    IsotopeDistribution iso;
    solver->setMaxIsotope(0);
    std::set<UInt> precursor_isotopes;
    // We're isolating the M0 and M+1 precursor isotopes
    precursor_isotopes.insert(0);
    precursor_isotopes.insert(1);

    // These are regression tests, but the results also follow an expected pattern.
    // See the comments in estimateForFragmentFromPeptideWeight for an explanation.
    iso = solver->estimateForFragmentFromDNAWeight(200.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.963845242419331)

    iso = solver->estimateForFragmentFromDNAWeight(2000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.978300783455351)

    iso = solver->estimateForFragmentFromDNAWeight(20000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995652529512413)

    iso = solver->estimateForFragmentFromDNAWeight(2000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.776727852910751)

    iso = solver->estimateForFragmentFromDNAWeight(20000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.95504592203456)

    iso = solver->estimateForFragmentFromDNAWeight(20000.0, 10000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.555730613643729)

    iso = solver->estimateForFragmentFromDNAWeight(200.0, 200.0, precursor_isotopes);
    IsotopeDistribution iso_precursor;
    solver->setMaxIsotope(2);
    iso_precursor = solver->estimateFromDNAWeight(200.0);
    IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

    for (; it1 != iso.end(); ++it1, ++it2)
    {
        TEST_EQUAL(it1->getMZ(), it2->getMZ())
        TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
    }
}
END_SECTION

START_SECTION(IsotopeDistribution CoarseIsotopePatternGenerator::estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
{
    IsotopeDistribution iso;
    solver->setMaxIsotope(0);
    std::set<UInt> precursor_isotopes;
    // We're isolating the M0 and M+1 precursor isotopes
    precursor_isotopes.insert(0);
    precursor_isotopes.insert(1);

    // These are regression tests, but the results also follow an expected pattern.
    // See the comments in estimateForFragmentFromPeptideWeight for an explanation.
    iso = solver->estimateForFragmentFromRNAWeight(200.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.963845242419331)

    iso = solver->estimateForFragmentFromRNAWeight(2000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.977854088814216)

    iso = solver->estimateForFragmentFromRNAWeight(20000.0, 100.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995465661923629)

    iso = solver->estimateForFragmentFromRNAWeight(2000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.784037437107401)

    iso = solver->estimateForFragmentFromRNAWeight(20000.0, 1000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.955768644474843)

    iso = solver->estimateForFragmentFromRNAWeight(20000.0, 10000.0, precursor_isotopes);
    iso.renormalize();
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.558201381343203)

    iso = solver->estimateForFragmentFromRNAWeight(200.0, 200.0, precursor_isotopes);
    IsotopeDistribution iso_precursor;
    solver->setMaxIsotope(2);
    iso_precursor = solver->estimateFromRNAWeight(200.0);
    IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

    for (; it1 != iso.end(); ++it1, ++it2)
    {
        TEST_EQUAL(it1->getMZ(), it2->getMZ())
        TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
    }
}
END_SECTION

START_SECTION(IsotopeDistribution calcFragmentIsotopeDist(const CoarseIsotopePatternGenerator & comp_fragment_isotope_distribution, const std::set<UInt>& precursor_isotopes, const double fragment_mono_mass))
{
  EmpiricalFormula ef_complementary_fragment = EmpiricalFormula("C2");
  EmpiricalFormula ef_fragment = EmpiricalFormula("C1");
  // The input to calcFragmentIsotopeDist should be isotope distributions
  // where the solver used atomic numbers for the mass field
  IsotopeDistribution iso1(ef_fragment.getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true))); // fragment
  IsotopeDistribution iso2(ef_complementary_fragment.getIsotopeDistribution(CoarseIsotopePatternGenerator(11, true))); // complementary fragment

  std::set<UInt> precursor_isotopes;
  precursor_isotopes.insert(0);
  precursor_isotopes.insert(1);

  precursor_isotopes.insert(2);
  IsotopeDistribution iso3;
  solver->setMaxIsotope(0);
  iso3 = solver->calcFragmentIsotopeDist(iso1, iso2, precursor_isotopes, ef_fragment.getMonoWeight());
  iso3.renormalize();

  // Need the distribution with accurate masses for the next comparison
  // because that's what the solver used for the fragment distribution
  IsotopeDistribution iso1_calc_mass(ef_fragment.getIsotopeDistribution(CoarseIsotopePatternGenerator(11))); // fragment

  IsotopeDistribution::ConstIterator it1(iso1_calc_mass.begin()), it2(iso3.begin());
  // By isolating all the precursor isotopes, the fragment isotopic distribution of a fragment molecule
  // should be the same as if it was the precursor. The probabilities can be slightly different due to
  // numerical issues.
  for (; it1 != iso1_calc_mass.end(); ++it1, ++it2)
  {
    TEST_EQUAL(it1->getMZ(), it2->getMZ())
    TEST_REAL_SIMILAR(it1->getIntensity(), it2->getIntensity())
  }

  precursor_isotopes.erase(precursor_isotopes.find(2));
  solver->setMaxIsotope(0);
  IsotopeDistribution iso4 = solver->calcFragmentIsotopeDist(iso1, iso2, precursor_isotopes, ef_fragment.getMonoWeight());
  iso4.renormalize();


  TEST_EQUAL(iso1_calc_mass.getContainer()[0].getMZ(), iso4.getContainer()[0].getMZ())
  TEST_EQUAL(iso1_calc_mass.getContainer()[1].getMZ(), iso4.getContainer()[1].getMZ())
  // Now that we're not isolating every precursor isotope, the probabilities should NOT be similar.
  // Since there's no TEST_REAL_NOT_SIMILAR, we test their similarity to the values they should be
  TEST_REAL_SIMILAR(iso1.getContainer()[0].getIntensity(), 0.989300)
  TEST_REAL_SIMILAR(iso1.getContainer()[1].getIntensity(), 0.010700)

  TEST_REAL_SIMILAR(iso4.getContainer()[0].getIntensity(), 0.989524)
  TEST_REAL_SIMILAR(iso4.getContainer()[1].getIntensity(), 0.010479)

  solver->setRoundMasses(true);
  IsotopeDistribution iso5 = solver->calcFragmentIsotopeDist(iso1, iso2, precursor_isotopes, ef_fragment.getMonoWeight());
  double result_mass[] = { 12.0, 13.0033548378 };
  double result_rounded_mass[] = { 12, 13 };
  Size i = 0;
  // making sure that the masses are correct depending on whether we asked the IsotopeDistribution solver
  // to return rounded masses
  for (it1 = iso3.begin(), it2 = iso5.begin(); it1 != iso3.end(); ++it1, ++it2, ++i)
  {
    TEST_REAL_SIMILAR(it1->getMZ(), result_mass[i])
    TEST_EQUAL(it2->getMZ(), result_rounded_mass[i])
  }
  solver->setRoundMasses(false);
}
END_SECTION

delete solver;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop
