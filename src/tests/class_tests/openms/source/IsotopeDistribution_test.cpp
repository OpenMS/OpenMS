// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

///////////////////////////

// This one is going to be tested.
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

START_TEST(IsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

IsotopeDistribution* nullPointer = nullptr;

START_SECTION(CoarseIsotopeDistribution())
	CoarseIsotopeDistribution* ptr = nullptr;
	ptr = new CoarseIsotopeDistribution();
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 0)
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION

START_SECTION(CoarseIsotopeDistribution(Size max_isotope))
	CoarseIsotopeDistribution* ptr = new CoarseIsotopeDistribution(117);
	Size max_isotope = ptr->getMaxIsotope();
  TEST_EQUAL(max_isotope, 117)
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION

CoarseIsotopeDistribution* iso = new CoarseIsotopeDistribution();

START_SECTION(CoarseIsotopeDistribution(const CoarseIsotopeDistribution& isotope_distribution))
	CoarseIsotopeDistribution copy;
	copy = *iso;
  for (Size i = 0; i != copy.getContainer().size(); ++i)
  {
    TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
    TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
  }
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
	TEST_EQUAL(copy.getMaxIsotope(), iso->getMaxIsotope())
END_SECTION

START_SECTION(~CoarseIsotopeDistribution())
	CoarseIsotopeDistribution* ptr = new CoarseIsotopeDistribution(117);
	delete ptr;
END_SECTION

START_SECTION(CoarseIsotopeDistribution& operator = (const CoarseIsotopeDistribution& isotope_distribution))
	CoarseIsotopeDistribution copy;
	copy = *iso;
	for (Size i = 0; i != copy.getContainer().size(); ++i)
	{
		TEST_EQUAL(copy.getContainer()[i].getMZ(), iso->getContainer()[i].getMZ())
		TEST_EQUAL(copy.getContainer()[i].getIntensity(), iso->getContainer()[i].getIntensity())
	}
	TEST_EQUAL(copy.getMin(), iso->getMin())
	TEST_EQUAL(copy.getMax(), iso->getMax())
	TEST_EQUAL(copy.size(), iso->size())
	TEST_EQUAL(copy.getMaxIsotope(), iso->getMaxIsotope())
END_SECTION

START_SECTION(void setMaxIsotope(Size max_isotope))
	CoarseIsotopeDistribution iso2;
	iso2.estimateFromPeptideWeight(1234.2);
	TEST_EQUAL(iso2.getMaxIsotope(), 0)
	TEST_EQUAL(iso2.getContainer().size(), 317)
	iso2.setMaxIsotope(117);
	TEST_EQUAL(iso2.getMaxIsotope(), 117)
END_SECTION

START_SECTION(Size getMaxIsotope() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(CoarseIsotopeDistribution operator + (const CoarseIsotopeDistribution& isotope_distribution) const)
	CoarseIsotopeDistribution iso1(1), iso2(1);
	CoarseIsotopeDistribution result = iso1 + iso2;
	TEST_EQUAL(result.size(), 1)
	CoarseIsotopeDistribution::ContainerType container = result.getContainer();
	TEST_EQUAL(container[0].getMZ(), 0)
	TEST_EQUAL(container[0].getIntensity(), 1)
END_SECTION

START_SECTION(CoarseIsotopeDistribution& operator *= (Size factor))
	EmpiricalFormula ef("C222N190O110");
	CoarseIsotopeDistribution id = ef.getIsotopeDistribution(new CoarseIsotopeDistribution(11));
	CoarseIsotopeDistribution::ContainerType container;
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7084, 0.0349429));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7085, 0.109888));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7086, 0.180185));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7087, 0.204395));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7088, 0.179765));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7089, 0.130358));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7090, 0.0809864));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7091, 0.0442441));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7092, 0.0216593));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7093, 0.00963707));
	container.push_back(CoarseIsotopeDistribution::MassAbundance(7094, 0.0039406));

	for (Size i = 0; i != id.size(); ++i)
	{
		TEST_EQUAL(round(id.getContainer()[i].getMZ()), container[i].getMZ())
		TEST_REAL_SIMILAR(id.getContainer()[i].getIntensity(), container[i].getIntensity())
	}

  // test gapped isotope distributions, e.g. bromide 79,81 (missing 80)
  {
    EmpiricalFormula ef("Br2");
    CoarseIsotopeDistribution id = ef.getIsotopeDistribution(new CoarseIsotopeDistribution(5));
    container.clear();
    // the expected results as pairs of
    // [nominal mass, probability]
    // derived via convolution of elemental probabilities; the sum of all probabilities is 1
    // For Br2, this is simply the product of Bromine x Bromine, which
    // has a light isotope (79 Da, ~50% probability) and a heavy isotope (81 Da, ~50% probability)
    container.push_back(CoarseIsotopeDistribution::MassAbundance(158, 0.2569476));  // 79+79, ~ 0.5 * 0.5
    container.push_back(CoarseIsotopeDistribution::MassAbundance(159, 0.0));        // this mass cannot be explained by two Br atoms
    container.push_back(CoarseIsotopeDistribution::MassAbundance(160, 0.49990478)); // 79+81 (or 81+79), ~ 0.5 * 0.5 + 0.5 * 0.5
    container.push_back(CoarseIsotopeDistribution::MassAbundance(161, 0.0));        // same as mass 159
    container.push_back(CoarseIsotopeDistribution::MassAbundance(162, 0.24314761)); // 81+81, ~ 0.5 * 0.5
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
    CoarseIsotopeDistribution id = ef.getIsotopeDistribution(new CoarseIsotopeDistribution(7));
    container.clear();
    container.push_back(CoarseIsotopeDistribution::MassAbundance(170, 0.254198270573));
    container.push_back(CoarseIsotopeDistribution::MassAbundance(171, 0.002749339427));
    container.push_back(CoarseIsotopeDistribution::MassAbundance(172, 0.494555798854));
    container.push_back(CoarseIsotopeDistribution::MassAbundance(173, 0.005348981146));
    container.push_back(CoarseIsotopeDistribution::MassAbundance(174, 0.240545930573));
    container.push_back(CoarseIsotopeDistribution::MassAbundance(175, 0.002601679427));
    for (Size i = 0; i != id.size(); ++i)
    {
      TEST_EQUAL(round(id.getContainer()[i].getMZ()), container[i].getMZ())
      TEST_REAL_SIMILAR(id.getContainer()[i].getIntensity(), container[i].getIntensity())
    }
  }

END_SECTION

START_SECTION(bool operator==(const CoarseIsotopeDistribution &isotope_distribution) const)
	CoarseIsotopeDistribution iso1(1);
	CoarseIsotopeDistribution iso2(2);
	TEST_EQUAL(iso1 == iso2, false)
	iso2.setMaxIsotope(1);
	TEST_EQUAL(iso1 == iso2, true)
	CoarseIsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11))),
    iso4(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso3 == iso4, true)
END_SECTION

START_SECTION(void set(const ContainerType &distribution))
	CoarseIsotopeDistribution iso1(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11))), iso2;
	TEST_EQUAL(iso1 == iso2, false)
	IsotopeDistribution::ContainerType container = iso1.getContainer();
	iso2.set(container);
	TEST_EQUAL(iso1.getContainer() == iso2.getContainer(), true)
	iso2.setMaxIsotope(iso1.getMaxIsotope());
	TEST_EQUAL(iso1 == iso2, true)
END_SECTION

START_SECTION(const ContainerType& getContainer() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(Size getMax() const)
	CoarseIsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso.getMax(), 6)
END_SECTION

START_SECTION(Size getMin() const)
	CoarseIsotopeDistribution iso(EmpiricalFormula("H2").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso.getMin(), 2)
	CoarseIsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso2.getMin(), 48)
END_SECTION

START_SECTION(Size size() const)
	CoarseIsotopeDistribution iso1, iso2(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso1.size(), 1)
	TEST_EQUAL(iso2.size(), 5)
END_SECTION

START_SECTION(void clear())
	CoarseIsotopeDistribution iso2(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso2.size(), 5)
	iso2.clear();
	TEST_EQUAL(iso2.size(), 0)
END_SECTION

START_SECTION(void estimateFromPeptideWeight(double average_weight))
	// hard to test as this is an rough estimate
	CoarseIsotopeDistribution iso(3);
	iso.estimateFromPeptideWeight(100.0);
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.949735)

	iso.estimateFromPeptideWeight(1000.0);
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.586906)

	iso.estimateFromPeptideWeight(10000.0);
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.046495)
END_SECTION

START_SECTION(void CoarseIsotopeDistribution::estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, UInt S_precursor, double average_weight_fragment, UInt S_fragment, const std::vector<UInt>& precursor_isotopes))
	CoarseIsotopeDistribution iso;
	CoarseIsotopeDistribution iso2;
	std::set<UInt> precursor_isotopes;
	// We're isolating the M+2 precursor isotopes
	precursor_isotopes.insert(2);
	// These are regression tests, but the results also follow an expected pattern.

	// With 0 sulfurs, it should be somewhat unlikely for the fragment to be M+2.
	iso.estimateForFragmentFromPeptideWeightAndS(200.0, 0, 100.0, 0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.355445559123552)

	// At such a small size, the regular averagine method should also result in 0 sulfurs.
	// The approximate EmpiricalFormulas should be the same, and therefore so should
	// their isotopic distributions.
	iso2.estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
	iso2.renormalize();

	CoarseIsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso2.begin());
	for (; it1 != iso.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

	// With the only sulfur being in the fragment, it's much more likely that the fragment
	// is M+2.
	iso.estimateForFragmentFromPeptideWeightAndS(200.0, 1, 100.0, 1, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.900804974056174)
	// Both sulfurs are in the fragment, so it's even more likely for the fragment to be M+2
	iso.estimateForFragmentFromPeptideWeightAndS(200.0, 2, 100.0, 2, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.947862830751023)
	// All 3 sulfurs are in the fragment
	iso.estimateForFragmentFromPeptideWeightAndS(200.0, 3, 100.0, 3, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.969454586761089)
	// Any more sulfurs while keeping the masses constant would violate the function preconditions

END_SECTION

START_SECTION(void CoarseIsotopeDistribution::estimateFromPeptideWeightAndS(double average_weight_precursor, UInt S))
	CoarseIsotopeDistribution iso(3);
	CoarseIsotopeDistribution iso2(3);
	// These are regression tests, but the results also follow an expected pattern.

	// With 0 sulfurs, it should be very unlikely for this tiny peptide to be M+2.
	iso.estimateFromPeptideWeightAndS(100.0, 0);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.00290370998965918)

	// At such a small size, the regular averagine method should also result in 0 sulfurs.
	// The approximate EmpiricalFormulas should be the same, and therefore so should
	// their isotopic distributions.
	iso2.estimateFromPeptideWeightAndS(100.0, 0);
	iso2.renormalize();

	CoarseIsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso2.begin());
	for (; it1 != iso.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ());
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

	// With one sulfur, it's more likely that the precursor is M+2 compared to having 0 sulfurs.
	iso.estimateFromPeptideWeightAndS(100.0, 1);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.0439547771832361)
	// With two sulfurs, the M+2 isotope is more likely
	iso.estimateFromPeptideWeightAndS(100.0, 2);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.0804989104418586)
	// With three sulfurs, the M+2 isotope is even more likely
	iso.estimateFromPeptideWeightAndS(100.0, 3);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.rbegin()->getIntensity(), 0.117023432503842)
	// Any more sulfurs while keeping the masses constant would violate the function preconditions

END_SECTION

START_SECTION(void estimateFromRNAWeight(double average_weight))
    // hard to test as this is an rough estimate
    CoarseIsotopeDistribution iso(3);
    iso.estimateFromRNAWeight(100.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.958166)

    iso.estimateFromRNAWeight(1000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.668538)

    iso.estimateFromRNAWeight(10000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.080505)
END_SECTION


START_SECTION(void estimateFromDNAWeight(double average_weight))
    // hard to test as this is an rough estimate
    CoarseIsotopeDistribution iso(3);
    iso.estimateFromDNAWeight(100.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.958166)

    iso.estimateFromDNAWeight(1000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.657083)

    iso.estimateFromDNAWeight(10000.0);
    TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.075138)
END_SECTION

START_SECTION(void estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
	CoarseIsotopeDistribution iso;
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
	iso.estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.954654801320083)

	// This peptide is large enough that the M0 and M+1 precursors are similar in abundance.
	// However, since the fragment is only 1/20th the mass of the precursor, it's
	// much more likely for the extra neutron to be on the complementary fragment.
	// Therefore, it's still much more likely this fragment will be monoisotopic.
	iso.estimateForFragmentFromPeptideWeight(2000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.975984866212216)

	// Same explanation as the previous example.
	iso.estimateForFragmentFromPeptideWeight(20000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995783521351781)

	// Like the first example, the fragment is half the mass of the precursor so
	// the fragments from the M+1 precursor will be roughly 50/50 monoisotopic/M+1.
	// However, this time the peptide is larger than the first example, so the M0
	// and M+1 precursors are also roughly 50/50 in abundance. Therefore, the
	// probability of the M0 fragment should be in the 75% range.
	//                  i.e. (100% * 50%) + (50% * 50%) = 75%
	// M0 frags due to M0 precursor^             ^M0 frags due to M+1 precursor
	iso.estimateForFragmentFromPeptideWeight(2000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.741290977639283)

	// Same explanation as the second example, except now the M+1 precursor is
	// more abundant than the M0 precursor. But, the fragment is so small that
	// it's still more likely for the fragment to be monoisotopic.
	iso.estimateForFragmentFromPeptideWeight(20000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.95467154987681)

	// Same explanation as above.
	iso.estimateForFragmentFromPeptideWeight(20000.0, 10000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.542260764523188)

	// If the fragment is identical to the precursor, then the distribution
	// should be the same as if it was just a precursor that wasn't isolated.
	iso.estimateForFragmentFromPeptideWeight(200.0, 200.0, precursor_isotopes);
	CoarseIsotopeDistribution iso_precursor(2);
	iso_precursor.estimateFromPeptideWeight(200.0);
	CoarseIsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

	for (; it1 != iso.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

END_SECTION

START_SECTION(void estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
	CoarseIsotopeDistribution iso;
	std::set<UInt> precursor_isotopes;
	// We're isolating the M0 and M+1 precursor isotopes
	precursor_isotopes.insert(0);
	precursor_isotopes.insert(1);

	// These are regression tests, but the results also follow an expected pattern.
	// See the comments in estimateForFragmentFromPeptideWeight for an explanation.
	iso.estimateForFragmentFromDNAWeight(200.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.963845242419331)

	iso.estimateForFragmentFromDNAWeight(2000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.978300783455351)

	iso.estimateForFragmentFromDNAWeight(20000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995652529512413)

	iso.estimateForFragmentFromDNAWeight(2000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.776727852910751)

	iso.estimateForFragmentFromDNAWeight(20000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.95504592203456)

	iso.estimateForFragmentFromDNAWeight(20000.0, 10000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.555730613643729)

	iso.estimateForFragmentFromDNAWeight(200.0, 200.0, precursor_isotopes);
	CoarseIsotopeDistribution iso_precursor(2);
	iso_precursor.estimateFromDNAWeight(200.0);
	CoarseIsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

	for (; it1 != iso.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

END_SECTION

START_SECTION(void estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes))
	CoarseIsotopeDistribution iso;
	std::set<UInt> precursor_isotopes;
	// We're isolating the M0 and M+1 precursor isotopes
	precursor_isotopes.insert(0);
	precursor_isotopes.insert(1);

	// These are regression tests, but the results also follow an expected pattern.
	// See the comments in estimateForFragmentFromPeptideWeight for an explanation.
	iso.estimateForFragmentFromRNAWeight(200.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.963845242419331)

	iso.estimateForFragmentFromRNAWeight(2000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.977854088814216)

	iso.estimateForFragmentFromRNAWeight(20000.0, 100.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.995465661923629)

	iso.estimateForFragmentFromRNAWeight(2000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.784037437107401)

	iso.estimateForFragmentFromRNAWeight(20000.0, 1000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.955768644474843)

	iso.estimateForFragmentFromRNAWeight(20000.0, 10000.0, precursor_isotopes);
	iso.renormalize();
	TEST_REAL_SIMILAR(iso.begin()->getIntensity(), 0.558201381343203)

	iso.estimateForFragmentFromRNAWeight(200.0, 200.0, precursor_isotopes);
	CoarseIsotopeDistribution iso_precursor(2);
	iso_precursor.estimateFromRNAWeight(200.0);
	CoarseIsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

	for (; it1 != iso.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

END_SECTION

START_SECTION(void estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes, double C, double H, double N, double O, double S, double P))
	// We are testing that the parameterized version matches the hardcoded version.
	CoarseIsotopeDistribution iso(3);
	CoarseIsotopeDistribution iso2(3);
	std::set<UInt> precursor_isotopes;
	precursor_isotopes.insert(0);
	precursor_isotopes.insert(1);

	iso.estimateForFragmentFromWeightAndComp(2000.0, 1000.0, precursor_isotopes, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0.0);
	iso2.estimateForFragmentFromPeptideWeight(2000.0, 1000.0, precursor_isotopes);
	TEST_EQUAL(iso.begin()->getIntensity(), iso2.begin()->getIntensity());
END_SECTION


START_SECTION(void estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P))
    // We are testing that the parameterized version matches the hardcoded version.
    CoarseIsotopeDistribution iso(3);
    CoarseIsotopeDistribution iso2(3);
    iso.estimateFromWeightAndComp(1000.0, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0.0);
    iso2.estimateFromPeptideWeight(1000.0);
    TEST_EQUAL(iso.begin()->getIntensity(),iso2.begin()->getIntensity());
END_SECTION

START_SECTION(void trimRight(double cutoff))
	CoarseIsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(new CoarseIsotopeDistribution(10)));
	TEST_NOT_EQUAL(iso.size(),3)
	iso.trimRight(0.2);
	TEST_EQUAL(iso.size(),3)
END_SECTION

START_SECTION(void trimLeft(double cutoff))
	CoarseIsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(new CoarseIsotopeDistribution(10)));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	TEST_EQUAL(iso.size(),2)
END_SECTION

START_SECTION(void renormalize())
	CoarseIsotopeDistribution iso(EmpiricalFormula("C160").getIsotopeDistribution(new CoarseIsotopeDistribution(10)));
	iso.trimRight(0.2);
	iso.trimLeft(0.2);
	iso.renormalize();
	double sum = 0;
	for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it)
	{
		sum += it->getIntensity();
	}

	TEST_REAL_SIMILAR(sum, 1.0)
END_SECTION

START_SECTION(IsotopeDistribution& operator+=(const CoarseIsotopeDistribution &isotope_distribution))
	CoarseIsotopeDistribution iso1(EmpiricalFormula("H1").getIsotopeDistribution(new CoarseIsotopeDistribution(11))),
											iso2(EmpiricalFormula("H2").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso1 == iso2, false)
	iso1 += CoarseIsotopeDistribution(EmpiricalFormula("H1").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso1.size() == iso2.size(), true)
	IsotopeDistribution::ConstIterator it1(iso1.begin()), it2(iso2.begin());

	for (; it1 != iso1.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it2->getIntensity(), it2->getIntensity())
	}

END_SECTION

START_SECTION(CoarseIsotopeDistribution operator *(Size factor) const)
	CoarseIsotopeDistribution iso1(EmpiricalFormula("H1").getIsotopeDistribution(new CoarseIsotopeDistribution(11))),
											iso2(EmpiricalFormula("H5").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
	TEST_EQUAL(iso1 == iso2, false)
	CoarseIsotopeDistribution iso3 = iso1 * 5;
	iso3.renormalize();
	iso2.renormalize();

	TEST_EQUAL(iso2.size(), iso3.size())
	IsotopeDistribution::ConstIterator it1(iso2.begin()), it2(iso3.begin());

	for (; it1 != iso2.end(); ++it1, ++it2)
	{
		TEST_EQUAL(it1->getMZ(), it2->getMZ())
		TEST_REAL_SIMILAR(it1->getIntensity(), it2->getIntensity())
	}
END_SECTION

START_SECTION(bool operator!=(const CoarseIsotopeDistribution &isotope_distribution) const)
  CoarseIsotopeDistribution iso1(1);
  CoarseIsotopeDistribution iso2(2);
  TEST_EQUAL(iso1 != iso2, true)
  iso2.setMaxIsotope(1);
  TEST_EQUAL(iso1 != iso2, false)
  CoarseIsotopeDistribution iso3(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11))),
                      iso4(EmpiricalFormula("C4").getIsotopeDistribution(new CoarseIsotopeDistribution(11)));
  TEST_EQUAL(iso3 != iso4, false)
END_SECTION

START_SECTION(CoarseIsotopeDistribution calcFragmentIsotopeDist(const CoarseIsotopeDistribution & comp_fragment_isotope_distribution, const std::set<UInt>& precursor_isotopes))
  CoarseIsotopeDistribution iso1(EmpiricalFormula("C1").getIsotopeDistribution(new CoarseIsotopeDistribution(11))); // fragment
  CoarseIsotopeDistribution iso2(EmpiricalFormula("C2").getIsotopeDistribution(new CoarseIsotopeDistribution(11))); // complementary fragment

  std::set<UInt> precursor_isotopes;
  precursor_isotopes.insert(0);
  precursor_isotopes.insert(1);

  precursor_isotopes.insert(2);
  CoarseIsotopeDistribution iso3;
  iso3.calcFragmentIsotopeDist(iso1,iso2,precursor_isotopes);
  iso3.renormalize();

  IsotopeDistribution::ConstIterator it1(iso1.begin()), it2(iso3.begin());
  // By isolating all the precursor isotopes, the fragment isotopic distribution of a fragment molecule
  // should be the same as if it was the precursor. The probabilities can be slightly different due to
  // numerical issues.
  for (; it1 != iso1.end(); ++it1, ++it2)
  {
	TEST_EQUAL(it1->getMZ(), it2->getMZ())
	TEST_REAL_SIMILAR(it1->getIntensity(), it2->getIntensity())
  }

  precursor_isotopes.erase(precursor_isotopes.find(2));
  CoarseIsotopeDistribution iso4;
  iso4.calcFragmentIsotopeDist(iso1,iso2,precursor_isotopes);
  iso4.renormalize();


  TEST_EQUAL(iso1.getContainer()[0].getMZ(), iso4.getContainer()[0].getMZ())
  TEST_EQUAL(iso1.getContainer()[1].getMZ(), iso4.getContainer()[1].getMZ())
  // Now that we're not isolating every precursor isotope, the probabilities should NOT be similar.
  // Since there's no TEST_REAL_NOT_SIMILAR, we test their similarity to the values they should be
  TEST_REAL_SIMILAR(iso1.getContainer()[0].getIntensity(), 0.989300)
  TEST_REAL_SIMILAR(iso1.getContainer()[1].getIntensity(), 0.010700)

  TEST_REAL_SIMILAR(iso4.getContainer()[0].getIntensity(), 0.989524)
  TEST_REAL_SIMILAR(iso4.getContainer()[1].getIntensity(), 0.010479)
END_SECTION

START_SECTION(Iterator begin())
	NOT_TESTABLE
END_SECTION

START_SECTION(Iterator end())
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator begin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator end() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rbegin())
	NOT_TESTABLE
END_SECTION

START_SECTION(ReverseIterator rend())
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rbegin() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(ConstReverseIterator rend() const)
	NOT_TESTABLE
END_SECTION

delete iso;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

