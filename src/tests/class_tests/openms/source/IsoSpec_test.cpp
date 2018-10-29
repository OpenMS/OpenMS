// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsoSpecWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//typedef std::pair<double,double> isopair;
typedef Peak1D isopair;
std::vector<isopair> expected;
expected.push_back(Peak1D( 180.06338810843999453936704 , 0.92211923150129093684768122    ));
expected.push_back(Peak1D( 181.06674294364000843415852 , 0.060338187213731186986365174   ));
expected.push_back(Peak1D( 181.06760524584001359471586 , 0.0021130960055518325027557047  ));
expected.push_back(Peak1D( 181.06966485435998492903309 , 0.0012805273467946340793660598  ));
expected.push_back(Peak1D( 182.06763310194000382580271 , 0.01137744132753025667892377    ));
expected.push_back(Peak1D( 182.07009777883999390724057 , 0.0016450768656347837786552146  ));
expected.push_back(Peak1D( 182.07096008103999906779791 , 0.00013826886808985885154825446 ));
expected.push_back(Peak1D( 182.07301968955999882382457 , 8.3790356109809577852924611e-05 ));
expected.push_back(Peak1D( 183.07098793713998929888476 , 0.00074447442519563354536987765 ));
expected.push_back(Peak1D( 183.0718502393399944594421  , 2.1726787058638343864448716e-05 ));
expected.push_back(Peak1D( 183.07345261403997938032262 , 2.3920974041512910987549237e-05 ));
expected.push_back(Peak1D( 183.07390984786002263717819 , 1.5799610569594269005293946e-05 ));
expected.push_back(Peak1D( 184.07187809544001311223838 , 5.8491247994869661650744336e-05 ));
expected.push_back(Peak1D( 184.07434277234000319367624 , 2.0297554674751325817869813e-05 ));

std::vector<isopair> expected_oms;
expected_oms.push_back(Peak1D( 180.06339038280000863778696 , 0.92263317941561073798339976    ));
expected_oms.push_back(Peak1D( 181.06674538280000774648215 , 0.059873700450778437331944559   ));
expected_oms.push_back(Peak1D( 181.06760738279999145561305 , 0.0021087279716237856790062022  ));
expected_oms.push_back(Peak1D( 181.06966713090000098418386 , 0.001273380225742650438680581   ));
expected_oms.push_back(Peak1D( 182.06764438280001172643097 , 0.011376032168236337518973933   ));
expected_oms.push_back(Peak1D( 182.07010038280000685517734 , 0.0016189442373783591386932068  ));
expected_oms.push_back(Peak1D( 182.07096238279999056430825 , 0.00013684457672024180033450158 ));
expected_oms.push_back(Peak1D( 182.07302213090000009287905 , 8.2635209633747614224076605e-05 ));
expected_oms.push_back(Peak1D( 183.07099938280001083512616 , 0.00073824045954083467209472236 ));
expected_oms.push_back(Peak1D( 183.07186138279999454425706 , 2.1667113372227351532611078e-05 ));
expected_oms.push_back(Peak1D( 183.07345538280000596387254 , 2.3346748674918047107952959e-05 ));
expected_oms.push_back(Peak1D( 183.07392113090000407282787 , 1.5700729969000005987024918e-05 ));
expected_oms.push_back(Peak1D( 184.07189838280001481507497 , 5.8444185791655326584186775e-05 ));
expected_oms.push_back(Peak1D( 184.07435438280000994382135 , 1.9961521148266482778097647e-05 ));

IsoSpecWrapper* ptr = nullptr;
IsoSpecWrapper* nullPointer = nullptr;
START_SECTION((IsoSpecThresholdWrapper(const EmpiricalFormula&, double, bool)))
  ptr = new IsoSpecThresholdWrapper(EmpiricalFormula("C10"), 0.5, false);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IsoSpecThresholdWrapper()))
  delete ptr;
END_SECTION

START_SECTION(( void run() ))
{
  double threshold = 1e-5;
  bool absolute = false;
  IsoSpecThresholdWrapper iso(EmpiricalFormula("C6H12O6"), threshold, absolute);
  std::vector<Peak1D> iso_result(iso.run());

  TEST_EQUAL(iso_result.size(), 14)

  // std::cout.precision(26);
  std::sort(iso_result.begin(), iso_result.end(),  [](isopair a, isopair b) {return a.getPos() < b.getPos();});

  for (Size i = 0; i != expected.size(); ++i)
  {
    TEST_REAL_SIMILAR(iso_result[i].getPos(), expected[i].getPos());
    TEST_REAL_SIMILAR(iso_result[i].getIntensity(), expected[i].getIntensity());
  }

  // human insulin
  std::vector<Peak1D> iso_result2 = IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), threshold, absolute).run();
  TEST_EQUAL(iso_result2.size(), 5402)

  std::vector<Peak1D> iso_result3 = IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), 0.01, false).run();
  TEST_EQUAL(iso_result3.size(), 269)
}
END_SECTION

START_SECTION(( [EXTRA] void run(const std::string&) ))
{
  double threshold = 1e-5;
  bool absolute = true;
  std::vector<Peak1D> iso_result(IsoSpecThresholdWrapper(EmpiricalFormula("C6H12O6"), threshold, absolute).run());

  TEST_EQUAL(iso_result.size(), 14)

  std::sort(iso_result.begin(), iso_result.end(),  [](isopair a, isopair b) {return a.getPos() < b.getPos();});

  for (Size i = 0; i != expected.size(); ++i)
  {
    TEST_REAL_SIMILAR(iso_result[i].getPos(), expected[i].getPos());
    TEST_REAL_SIMILAR(iso_result[i].getIntensity(), expected[i].getIntensity());
  }

  // human insulin
  std::vector<Peak1D> iso_result2(IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), threshold, absolute).run());
  TEST_EQUAL(iso_result2.size(), 1720)

  std::vector<Peak1D> iso_result3(IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), 0.01, true).run());
  TEST_EQUAL(iso_result3.size(), 21)
}
END_SECTION

START_SECTION(( 
    void run(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities) ))
{

  EmpiricalFormula ef ("C6H12O6");

  std::vector<int> isotopeNumbers;
  std::vector<int> atomCounts;
  std::vector<std::vector<double> > isotopeMasses;
  std::vector<std::vector<double> > isotopeProbabilities;

  for (auto elem : ef)
  {
    atomCounts.push_back( elem.second );

    std::vector<double> masses;
    std::vector<double> probs;
    for (auto iso : elem.first->getIsotopeDistribution())
    {
      if (iso.getIntensity() <= 0.0) continue; // Note: there will be a segfault if one of the intensities is zero!
      masses.push_back(iso.getMZ());
      probs.push_back(iso.getIntensity());
    }
    isotopeNumbers.push_back( masses.size() );
    isotopeMasses.push_back(masses);
    isotopeProbabilities.push_back(probs);
  }

  // ----------------------------------- 
  // Start
  // ----------------------------------- 
  {
    double threshold = 1e-5;
    bool absolute = false;
    std::vector<Peak1D> iso_results(IsoSpecThresholdWrapper(isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities, threshold, absolute).run());

    TEST_EQUAL(iso_results.size(), 14)

    std::sort(iso_results.begin(), iso_results.end(),  [](isopair a, isopair b) {return a.getPos() < b.getPos();});

    for (Size i = 0; i != expected_oms.size(); ++i)
    {
      TEST_REAL_SIMILAR(iso_results[i].getPos(), expected_oms[i].getPos());
      TEST_REAL_SIMILAR(iso_results[i].getIntensity(), expected_oms[i].getIntensity());
    }
  }

  // TEST exception:
  // We cannot have zero values as input data
  double threshold = 1e-5;
  isotopeNumbers[0] += 1;
  isotopeMasses[0].push_back(3.0160492699999998933435563);
  isotopeProbabilities[0].push_back(0.0);
  TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecThresholdWrapper(isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities, threshold, false).run());

}
END_SECTION

#if 0
START_SECTION(( [STRESSTEST] void run(const std::string&) ))
{
  // Do some stress testing of the library...
  // this is close to the performance of IsoSpec by itself
  int sum = 0;
  for (Size k = 0; k < 2e5; k++)
  {
    double threshold = 1e-2;
    bool absolute = false;
    IsoSpecThresholdWrapper iso("C520H817N139O147", threshold, absolute);
    auto res = iso.run();
    sum += res.size();
  }
  TEST_EQUAL(sum, 140*2*1e5)
}
END_SECTION
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
