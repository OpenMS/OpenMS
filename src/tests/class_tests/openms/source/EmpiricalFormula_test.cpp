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
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ElementDB, "$Id$")

/////////////////////////////////////////////////////////////

EmpiricalFormula* e_ptr = nullptr;
EmpiricalFormula* e_nullPointer = nullptr;
const ElementDB * db = ElementDB::getInstance();

START_SECTION(EmpiricalFormula())
  e_ptr = new EmpiricalFormula;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~EmpiricalFormula())
  delete e_ptr;
END_SECTION

START_SECTION(EmpiricalFormula(const String& rhs))
  e_ptr = new EmpiricalFormula("C4");
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
        EmpiricalFormula e0("C5(13)C4H2");
        EmpiricalFormula e1("C5(13)C4");
        EmpiricalFormula e2("(12)C5(13)C4");
        EmpiricalFormula e3("C9");
  TEST_REAL_SIMILAR(e1.getMonoWeight(), e2.getMonoWeight())
  TEST_REAL_SIMILAR(e1.getMonoWeight(), 112.013419)
  TEST_REAL_SIMILAR(e2.getMonoWeight(), 112.013419)
END_SECTION

START_SECTION(EmpiricalFormula(const EmpiricalFormula& rhs))
  EmpiricalFormula ef(*e_ptr);
  TEST_EQUAL(ef == *e_ptr, true)
END_SECTION

START_SECTION((EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge=0)))
  EmpiricalFormula ef(4, db->getElement("C"));
  TEST_EQUAL(ef == *e_ptr, true)
  TEST_EQUAL(ef.getCharge(), 0)
END_SECTION

START_SECTION(const Element* getElement(UInt atomic_number) const)
  const Element* e = db->getElement(6);
  TEST_EQUAL(e->getSymbol(), "C")
END_SECTION

START_SECTION(const Element* getElement(const String& name) const)
  const Element* e = db->getElement("C");
  TEST_EQUAL(e->getSymbol(), "C")
END_SECTION

START_SECTION(SignedSize getNumberOf(const Element* element) const)
  Size num1 = e_ptr->getNumberOf(db->getElement(6));
  TEST_EQUAL(num1, 4);

  Size num2 = e_ptr->getNumberOf(db->getElement("C"));
  TEST_EQUAL(num2, 4);
END_SECTION

START_SECTION(SignedSize getNumberOfAtoms() const)
  Size num4 = e_ptr->getNumberOfAtoms();
  TEST_EQUAL(num4, 4);
END_SECTION

START_SECTION(EmpiricalFormula& operator = (const EmpiricalFormula& rhs))
  EmpiricalFormula ef;
  ef = *e_ptr;
  TEST_EQUAL(*e_ptr == ef, true)
END_SECTION

START_SECTION(EmpiricalFormula operator * (const SignedSize& times) const)
  EmpiricalFormula ef("C3H8");
  ef = ef * 3;
  TEST_EQUAL(ef, "C9H24")
END_SECTION

START_SECTION(EmpiricalFormula& operator += (const EmpiricalFormula& rhs))
  EmpiricalFormula ef("C3");
  ef += ef;
  TEST_EQUAL(ef, "C6")
  EmpiricalFormula ef2("C-6H2");
  ef += ef2;
  TEST_EQUAL(ef, "H2");

  ef = EmpiricalFormula("C");
  TEST_EQUAL(ef, "C")
  ef += EmpiricalFormula("C5");
  TEST_EQUAL(ef, "C6")
  ef += EmpiricalFormula("C-5");
  TEST_EQUAL(ef, "C")
  ef += EmpiricalFormula("C-1H2");
  TEST_EQUAL(ef, "H2")
END_SECTION

START_SECTION(EmpiricalFormula operator + (const EmpiricalFormula& rhs) const)
  EmpiricalFormula ef("C2");
  EmpiricalFormula ef2;
  ef2 = ef + ef;
  TEST_EQUAL(ef2, "C4")
  ef2 = ef2 + EmpiricalFormula("C-4H2");
  TEST_EQUAL(ef2, "H2")
END_SECTION

START_SECTION(EmpiricalFormula& operator -= (const EmpiricalFormula& rhs))
  EmpiricalFormula ef1("C5H12"), ef2("CH12");
  ef1 -= ef2;
  TEST_EQUAL(*e_ptr == ef1, true)
  ef1 -= EmpiricalFormula("C4H-2");
  TEST_EQUAL(ef1, "H2");
END_SECTION

START_SECTION(EmpiricalFormula operator - (const EmpiricalFormula& rhs) const)
  EmpiricalFormula ef1("C5H12"), ef2("CH12");
  EmpiricalFormula ef3, ef4;
  ef3 = ef1 - ef2;
  cerr << *e_ptr << " " << ef3 << endl;
  TEST_EQUAL(*e_ptr == ef3, true)
  ef3 = ef3 - EmpiricalFormula("C4H-2");
  TEST_EQUAL(ef3, "H2");
END_SECTION

START_SECTION(bool isEmpty() const)
  EmpiricalFormula ef;
  TEST_EQUAL(ef.isEmpty(), true)
  TEST_EQUAL(e_ptr->isEmpty(), false)
END_SECTION

START_SECTION(bool hasElement(const Element* element) const)
  const Element* e = db->getElement(6);
  TEST_EQUAL(e_ptr->hasElement(e), true)
  e = db->getElement(1);
  TEST_EQUAL(e_ptr->hasElement(e), false)
END_SECTION

START_SECTION(bool contains(const EmpiricalFormula& ef))

  EmpiricalFormula metabolite("C12H36N2");

  TEST_EQUAL(metabolite.contains(metabolite), true) // contains itself?
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("C-12H36N2")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("C11H36N2")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("N2")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("H36")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("H3")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("P-1")), true)
  TEST_EQUAL(metabolite.contains(EmpiricalFormula()), true)

  TEST_EQUAL(metabolite.contains(EmpiricalFormula("P1")), false)

  // the 'adduct' test
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("KH-2") * -1), true) // make sure we can loose 2H (i.e. we have 2H in the metabolite); K is adducted, so is does not need to be intrinsic
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("K-1H2")), true) // same as above
  TEST_EQUAL(metabolite.contains(EmpiricalFormula("KH-2") * 1), false) // cannot loose K, since we don't have it
END_SECTION

START_SECTION(void setCharge(SignedSize charge))
  e_ptr->setCharge(1);
  NOT_TESTABLE // will be tested in next check
END_SECTION

START_SECTION(SignedSize getCharge() const)
  TEST_EQUAL(e_ptr->getCharge(), 1)
  EmpiricalFormula ef1("C2+");
  TEST_EQUAL(ef1.getCharge(), 1)
  EmpiricalFormula ef2("C2+3");
  TEST_EQUAL(ef2.getCharge(), 3)
END_SECTION

START_SECTION(bool isCharged() const)
  TEST_EQUAL(e_ptr->isCharged(), true)
  e_ptr->setCharge(0);
  TEST_EQUAL(e_ptr->isCharged(), false)
END_SECTION

START_SECTION(double getAverageWeight() const)
  EmpiricalFormula ef("C2");
  const Element* e = db->getElement("C");
  TEST_REAL_SIMILAR(ef.getAverageWeight(), e->getAverageWeight() * 2)
END_SECTION

START_SECTION(bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P))
    // Same stoichiometry as the averagine model
    EmpiricalFormula ef("C494H776N136O148S4");
    EmpiricalFormula ef_approx;
    Int return_flag;
    return_flag = ef_approx.estimateFromWeightAndComp(ef.getAverageWeight(), 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
    // average mass should be the same when using the same stoichiometry
    TOLERANCE_ABSOLUTE(1);
    TEST_REAL_SIMILAR(ef.getAverageWeight(), ef_approx.getAverageWeight());
    // # of elements should be the same when using the same stoichiometry
    for (EmpiricalFormula::ConstIterator itr = ef.begin(); itr != ef.end(); ++itr)
    {
      TEST_EQUAL(itr->second, ef_approx.getNumberOf(itr->first));
    }
    TEST_EQUAL(return_flag, true);

    // Very different stoichiometry than the averagine model
    EmpiricalFormula ef2("C100H100N100O100S100P100");
    return_flag = ef_approx.estimateFromWeightAndComp(ef2.getAverageWeight(), 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
    // average mass should be the same when using a different stoichiometry
    TEST_REAL_SIMILAR(ef2.getAverageWeight(), ef_approx.getAverageWeight());
    // # of elements should be different when using a very different stoichiometry
    for (EmpiricalFormula::ConstIterator itr = ef2.begin(); itr != ef2.end(); ++itr)
    {
      TEST_NOT_EQUAL(itr->second, ef_approx.getNumberOf(itr->first));
    }
    TEST_EQUAL(return_flag, true);

    // Small mass that the model can't fit without using a negative # of hydrogens
    return_flag = ef_approx.estimateFromWeightAndComp(50, 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
    // The same mass can't be achieved because the # hydrogens needed to compensate is negative
    TEST_EQUAL(ef_approx.getAverageWeight() - 50 > 1, true);
    // Don't allow the EmpiricalFormula to have a negative # of hydrogens
    TEST_EQUAL(ef_approx.getNumberOf(db->getElement("H")) >= 0, true);
    // The return flag should now indicate that the estimated formula did not succeed without requesting a negative # of hydrogens
    TEST_EQUAL(return_flag, false);

END_SECTION

START_SECTION(bool estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P))
    EmpiricalFormula ef("C494H776N136O148S4");
    EmpiricalFormula ef_approx;
    EmpiricalFormula ef_approx_S;
    bool return_flag;
    // Using averagine stoichiometry, excluding sulfur.
    return_flag = ef_approx_S.estimateFromWeightAndCompAndS(ef.getAverageWeight(), 4, 4.9384, 7.7583, 1.3577, 1.4773, 0);
    TEST_EQUAL(4, ef_approx_S.getNumberOf(db->getElement("S")));

    // Formula of methionine.
    EmpiricalFormula ef2("C5H9N1O1S1");
    // Using averagine stoichiometry, excluding sulfur.
    return_flag = ef_approx_S.estimateFromWeightAndCompAndS(ef2.getAverageWeight(), 1, 4.9384, 7.7583, 1.3577, 1.4773, 0);
    // Shouldn't need negative hydrogens for this approximation.
    TEST_EQUAL(return_flag, true);
    ef_approx.estimateFromWeightAndComp(ef2.getAverageWeight(), 4.9384, 7.7583, 1.3577, 1.4773, 0.0417, 0);
    // The averagine approximation should result in 0 sulfurs.
    TEST_EQUAL(0, ef_approx.getNumberOf(db->getElement("S")));
    // But with the sulfur-specified averagine version, we forced it be 1
    TEST_EQUAL(1, ef_approx_S.getNumberOf(db->getElement("S")));
    TOLERANCE_ABSOLUTE(1);
    TEST_REAL_SIMILAR(ef_approx.getAverageWeight(), ef_approx_S.getAverageWeight());

END_SECTION

START_SECTION(double getMonoWeight() const)
  EmpiricalFormula ef("C2");
  const Element* e = db->getElement("C");
  TEST_REAL_SIMILAR(ef.getMonoWeight(), e->getMonoWeight() * 2)
  TEST_REAL_SIMILAR(EmpiricalFormula("OH").getMonoWeight(), EmpiricalFormula("HO").getMonoWeight());
  TEST_REAL_SIMILAR(EmpiricalFormula("").getMonoWeight(), 0.0)
END_SECTION

START_SECTION(String toString() const)
  EmpiricalFormula ef("C2H5");
  String str = ef.toString();
  TEST_EQUAL(String(str).hasSubstring("H5"), true)
  TEST_EQUAL(String(str).hasSubstring("C2"), true)
END_SECTION

START_SECTION([EXTRA](friend std::ostream& operator << (std::ostream&, const EmpiricalFormula&)))
  stringstream ss;
  EmpiricalFormula ef("C2H5");
  ss << ef;
  TEST_EQUAL(String(ss.str()).hasSubstring("H5"), true);
  TEST_EQUAL(String(ss.str()).hasSubstring("C2"), true);
END_SECTION

START_SECTION(bool operator != (const EmpiricalFormula& rhs) const)
  EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
  TEST_EQUAL(ef1 != ef2, true)
  TEST_EQUAL(ef1 != ef1, false)
  ef2.setCharge(1);
  TEST_EQUAL(ef2 != *e_ptr, true)
END_SECTION

START_SECTION(bool operator == (const EmpiricalFormula& rhs) const)
  EmpiricalFormula ef1("C2H5"), ef2(*e_ptr);
  TEST_EQUAL(ef1 == ef2, false)
  TEST_EQUAL(ef1 == ef1, true)
  ef2.setCharge(1);
  TEST_EQUAL(ef2 == *e_ptr, false)
END_SECTION

START_SECTION(ConstIterator begin() const)
  EmpiricalFormula ef("C6H12O6");
  Map<String, SignedSize> formula;
  formula["C"] = 6;
  formula["H"] = 12;
  formula["O"] = 6;
  for (EmpiricalFormula::ConstIterator it = ef.begin(); it != ef.end(); ++it)
  {
    TEST_EQUAL(it->second, formula[it->first->getSymbol()])
  }

END_SECTION

START_SECTION(ConstIterator end() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(IsotopeDistribution getIsotopeDistribution(UInt max_depth) const)
  EmpiricalFormula ef("C");
  IsotopeDistribution iso = ef.getIsotopeDistribution(new CoarseIsotopeDistribution(20));
  double result[] = { 0.9893, 0.0107};
  Size i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result[i])
  }
END_SECTION

START_SECTION(IsotopeDistribution getConditionalFragmentIsotopeDist(const EmpiricalFormula& precursor, const std::set<UInt>& precursor_isotopes) const)
  EmpiricalFormula precursor("C2");
  EmpiricalFormula fragment("C");
  std::set<UInt> precursor_isotopes;

  precursor_isotopes.insert(0);
  // isolated precursor isotope is M0
  IsotopeDistribution iso = fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  double result[] = { 1.0 };
  Size i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result[i])
  }

  precursor_isotopes.clear();
  precursor_isotopes.insert(1);
  // isolated precursor isotope is M+1
  iso = fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  double result2[] = { 0.5, 0.5};
  i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result2[i])
  }

  precursor_isotopes.insert(0);
  // isolated precursor isotopes are M0 and M+1
  iso = fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  double result3[] = { 0.98941, 0.01059};
  i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result3[i])
  }

  precursor_isotopes.insert(2);
  // isolated precursor isotopes are M0, M+1, and M+2
  // This is the example found in the comments of the getConditionalFragmentIsotopeDist function.
  // Since we're isolating all the possible precursor isotopes, the fragment isotope distribution
  // should be equivalent to the natural isotope abundances.
  iso = fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  double result4[] = { 0.9893, 0.0107};
  i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result4[i])
  }

  precursor_isotopes.insert(3);
  // isolated precursor isotopes are M0, M+1, M+2, and M+3
  // It's impossible for precursor C2 to have 3 extra neutrons (assuming only natural stable isotopes)
  // Invalid precursor isotopes are ignored and should give the answer as if they were not there
  iso = fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  i = 0;
  for (IsotopeDistribution::ConstIterator it = iso.begin(); it != iso.end(); ++it, ++i)
  {
    TEST_REAL_SIMILAR(it->getIntensity(), result4[i])
  }

  precursor = EmpiricalFormula("C10H10N10O10S2");
  EmpiricalFormula big_fragment = EmpiricalFormula("C9H9N9O9S1");
  EmpiricalFormula small_fragment = EmpiricalFormula("C1H1N1O1S1");

  precursor_isotopes.clear();
  precursor_isotopes.insert(1);
  // isolated precursor isotope is M+1
  IsotopeDistribution big_iso = big_fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);
  IsotopeDistribution small_iso = small_fragment.getConditionalFragmentIsotopeDist(precursor,precursor_isotopes);

  // When we isolate only the M+1 precursor isotope, the big_fragment is more likely to exist as M+1 than M0.
  TEST_EQUAL(big_iso.getContainer()[0].getIntensity() < 0.2, true)
  TEST_EQUAL(big_iso.getContainer()[1].getIntensity() > 0.8, true)

  // The small_fragment, however, is more likely to exist as M0 than M+1.
  TEST_EQUAL(small_iso.getContainer()[0].getIntensity() > 0.8, true)
  TEST_EQUAL(small_iso.getContainer()[1].getIntensity() < 0.2, true)

  // Since the two fragments also happen to be complementary, their probabilities are perfectly reversed.
  IsotopeDistribution::ConstIterator big_it = big_iso.begin();
  IsotopeDistribution::ConstReverseIterator small_it = small_iso.rbegin();
  for (; big_it != big_iso.end(); ++big_it, ++small_it)
  {
    TEST_REAL_SIMILAR(big_it->getIntensity(), small_it->getIntensity())
  }

END_SECTION

START_SECTION(([EXTRA] Check correct charge semantics))
  EmpiricalFormula ef1("H4C+"); // CH4 +1 charge
  const Element * H = db->getElement("H");
  const Element * C = db->getElement("C");

  TEST_EQUAL(ef1.getNumberOf(H), 4)
  TEST_EQUAL(ef1.getNumberOf(C), 1)
  TEST_EQUAL(ef1.getCharge(), 1)
  EmpiricalFormula ef2("H4C1+"); // ""
  TEST_EQUAL(ef2.getNumberOf(H), 4)
  TEST_EQUAL(ef2.getNumberOf(C), 1)
  TEST_EQUAL(ef2.getCharge(), 1)
  EmpiricalFormula ef3("H4C-1+"); // C-1 H4 +1 charge
  TEST_EQUAL(ef3.getNumberOf(H), 4)
  TEST_EQUAL(ef3.getNumberOf(C), -1)
  TEST_EQUAL(ef3.getCharge(), 1)
  EmpiricalFormula ef4("H4C-1"); // C-1 H4 0 charge
  TEST_EQUAL(ef4.getNumberOf(H), 4)
  TEST_EQUAL(ef4.getNumberOf(C), -1)
  TEST_EQUAL(ef4.getCharge(), 0)
  EmpiricalFormula ef5("H4C1-1"); // C1 H4 -1 charge
  TEST_EQUAL(ef5.getNumberOf(H), 4)
  TEST_EQUAL(ef5.getNumberOf(C), 1)
  TEST_EQUAL(ef5.getCharge(), -1)
  EmpiricalFormula ef6("H4C-1-1"); // C-1 H4 -1 charge
  TEST_EQUAL(ef6.getNumberOf(H), 4)
  TEST_EQUAL(ef6.getNumberOf(C), -1)
  TEST_EQUAL(ef6.getCharge(), -1)
  EmpiricalFormula ef7("H4C-1-"); // C-1 H4 -1 charge
  TEST_EQUAL(ef7.getNumberOf(H), 4)
  TEST_EQUAL(ef7.getNumberOf(C), -1)
  TEST_EQUAL(ef7.getCharge(), -1)
  EmpiricalFormula ef8("-"); // -1 Charge
  TEST_EQUAL(ef8.getNumberOf(H), 0)
  TEST_EQUAL(ef8.getNumberOf(C), 0)
  TEST_EQUAL(ef8.getCharge(), -1)
  EmpiricalFormula ef9("+"); // +1 Charge
  TEST_EQUAL(ef9.getNumberOf(H), 0)
  TEST_EQUAL(ef9.getNumberOf(C), 0)
  TEST_EQUAL(ef9.getCharge(), 1)
  EmpiricalFormula ef10("-3"); // -3 Charge
  TEST_EQUAL(ef10.getNumberOf(H), 0)
  TEST_EQUAL(ef10.getNumberOf(C), 0)
  TEST_EQUAL(ef10.getCharge(), -3)
  EmpiricalFormula ef11("+3"); // +3 Charge
  TEST_EQUAL(ef11.getNumberOf(H), 0)
  TEST_EQUAL(ef11.getNumberOf(C), 0)
  TEST_EQUAL(ef11.getCharge(), 3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

