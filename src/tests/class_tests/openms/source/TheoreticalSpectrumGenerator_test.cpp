// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(TheoreticalSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TheoreticalSpectrumGenerator* ptr = 0;
TheoreticalSpectrumGenerator* nullPointer = 0;

START_SECTION(TheoreticalSpectrumGenerator())
  ptr = new TheoreticalSpectrumGenerator();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source))
  TheoreticalSpectrumGenerator copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~TheoreticalSpectrumGenerator())
  delete ptr;
END_SECTION

ptr = new TheoreticalSpectrumGenerator();
AASequence peptide = AASequence::fromString("IFSQVGK");

START_SECTION(TheoreticalSpectrumGenerator& operator = (const TheoreticalSpectrumGenerator& tsg))
  TheoreticalSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void addPeaks(RichPeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, Int charge = 1))
  RichPeakSpectrum y_spec, b_spec, a_spec;
  ptr->addPeaks(y_spec, peptide, Residue::YIon, 1);
  ptr->addPeaks(b_spec, peptide, Residue::BIon, 1);
  ptr->addPeaks(a_spec, peptide, Residue::AIon, 1);
  TOLERANCE_ABSOLUTE(0.001)
  double y_result[] = {147.113, 204.135, 303.203, 431.262, 518.294, 665.362};
  for (Size i = 0; i != y_spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(y_spec[i].getPosition()[0], y_result[i])
  }
  double b_result[] = {/*115.1,*/ 261.16, 348.192, 476.251, 575.319, 632.341};
  for (Size i = 0; i != b_spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(b_spec[i].getPosition()[0], b_result[i])
  }

  double a_result[] = {/*87.1048,*/ 233.165, 320.197, 448.256, 547.324, 604.346};
  for (Size i = 0; i != a_spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(a_spec[i].getPosition()[0], a_result[i])
  }

  RichPeakSpectrum y_spec2;
  ptr->addPeaks(y_spec2, peptide, Residue::YIon, 2);
  TOLERANCE_ABSOLUTE(0.01)
  for (Size i = 0; i != y_spec2.size(); ++i)
  {
    TEST_REAL_SIMILAR(y_spec2[i].getPosition()[0], (y_result[i]+1.0)/2.0)
  }
END_SECTION

START_SECTION(void addAbundantImmoniumIons(RichPeakSpectrum& spec))
  RichPeakSpectrum spec;
  ptr->addAbundantImmoniumIons(spec);
  TEST_EQUAL(spec.size(), 5)
  TEST_REAL_SIMILAR(spec[0].getPosition()[0], 86.09698)
  TEST_REAL_SIMILAR(spec[1].getPosition()[0], 110.0718)
  TEST_REAL_SIMILAR(spec[2].getPosition()[0], 120.0813)
  TEST_REAL_SIMILAR(spec[3].getPosition()[0], 136.0762)
  TEST_REAL_SIMILAR(spec[4].getPosition()[0], 159.0922)
END_SECTION


START_SECTION(void addPrecursorPeaks(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
  RichPeakSpectrum spec;
  ptr->addPrecursorPeaks(spec, peptide, 1);
  double result[] = {760.4352, 761.4192, 778.4457};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  RichPeakSpectrum spec2;
  ptr->addPrecursorPeaks(spec2, peptide, 2);
  double result2[] = {380.7212, 381.2132, 389.7265};
  for (Size i = 0; i != spec2.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec2[i].getPosition()[0], result2[i])
  }

END_SECTION

START_SECTION(void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
  RichPeakSpectrum spec;
  ptr->getSpectrum(spec, peptide, 1);
  TEST_EQUAL(spec.size(), 11)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {/*114.091,*/ 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getSpectrum(spec, peptide, 2);
  TEST_EQUAL(spec.size(), 22)

  spec.clear(true);
  Param param(ptr->getParameters());
  param.setValue("add_first_prefix_ion", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, peptide, 1);
  TEST_EQUAL(spec.size(), 12)

  double result2[] = {114.091, 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result2[i])
  }


  AASequence new_peptide = AASequence::fromString("DFPLANGER");
  /**  From http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
   Seq    #       A            B            C            X            Y            Z         # (+1) 
    D     1     88.03990    116.03481    133.06136       -        1018.49583   1001.46928    9 
    F     2    235.10831    263.10323    280.12978    929.44815    903.46888    886.44233    8 
    P     3    332.16108    360.15599    377.18254    782.37973    756.40047    739.37392    7 
    I     4    445.24514    473.24005    490.26660    685.32697    659.34771    642.32116    6 
    A     5    516.28225    544.27717    561.30372    572.24291    546.26364    529.23709    5 
    N     6    630.32518    658.32009    675.34664    501.20579    475.22653    458.19998    4 
    G     7    687.34664    715.34156    732.36811    387.16287    361.18360    344.15705    3 
    E     8    816.38924    844.38415    861.41070    330.14140    304.16214    287.13559    2 
    R     9    972.49035   1000.48526       -         201.09881    175.11955    158.09300    1 
  **/

  double result_all[Size(52)] = { 88.03990,116.03481,133.06136,1018.49583,1001.46928,235.10831,263.10323,280.12978,929.44815,903.46888,886.44233,
    332.16108,360.15599,377.18254,782.37973,756.40047,739.37392, 445.24514,473.24005,490.26660,685.32697,659.34771,642.32116,
    516.28225,544.27717,561.30372,572.24291,546.26364,529.23709, 630.32518,658.32009,675.34664,501.20579,475.22653,458.19998,
    687.34664,715.34156,732.36811,387.16287,361.18360,344.15705, 816.38924,844.38415,861.41070,330.14140,304.16214,287.13559,
    /*972.49035 basically neutral loss of prec,*/1000.48526,201.09881,175.11955,158.09300};
  std::sort(result_all,result_all+Size(51));
  spec.clear(true);

  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_precursor_peaks", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, new_peptide, 1);
  TEST_EQUAL(spec.size(), 51)

  vector<double> generated;
  for (Size i = 0; i != spec.size(); ++i)
  {
    generated.push_back(spec[i].getPosition()[0]);
  }

  std::sort(generated.begin(),generated.end());
  for (Size i = 0; i != generated.size(); ++i)
  {
    TEST_REAL_SIMILAR(generated[i], result_all[i])
  }
  


END_SECTION

START_SECTION(([EXTRA] bugfix test where losses lead to formulae with negative element frequencies))
{
  AASequence tmp_aa = AASequence::fromString("RDAGGPALKK");
  RichPeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param params;

  params.setValue("add_isotopes", "true");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  t_gen.setParameters(params);

  t_gen.getSpectrum(tmp, tmp_aa,1);
  t_gen.addPeaks(tmp, tmp_aa, Residue::AIon);
  TEST_EQUAL(tmp.size(), 212)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
