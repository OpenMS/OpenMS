// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumAnnotator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAnnotator* ptr = nullptr;
SpectrumAnnotator* null_ptr = nullptr;
START_SECTION(SpectrumAnnotator())
{
	ptr = new SpectrumAnnotator();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~SpectrumAnnotator())
{
	delete ptr;
}
END_SECTION

TheoreticalSpectrumGenerator tg = TheoreticalSpectrumGenerator();
Param tgp(tg.getDefaults());
tgp.setValue("add_metainfo", "true");
tgp.setValue("add_y_ions", "true");
tgp.setValue("add_b_ions", "true");
tg.setParameters(tgp);
SpectrumAlignment sa;
Param sap = sa.getDefaults();
sap.setValue("tolerance", 0.1, "...");
sa.setParameters(sap);
SpectrumAnnotator annot;
AASequence peptide = AASequence::fromString("IFSQVGK");
PeptideHit hit;
hit.setSequence(peptide);
hit.setCharge(2);
PeakSpectrum spec;
spec.setMSLevel(2);
Peak1D p;
double peaklist[] = {147.113, 204.135, 303.203, 431.262, 518.294, 665.362, 261.16, 348.192, 476.251, 575.319, 632.341};
size_t pls = 11; // peaklist size

for (Size i = 0; i != pls; ++i)
{
  p.setIntensity(1.1f);
  p.setMZ(peaklist[i]);
  spec.push_back(p);
}
PeptideIdentification pi;
pi.setHits(std::vector<PeptideHit>(1,hit));

START_SECTION((void SpectrumAnnotator::annotateMatches(PeakSpectrum& spec, const PeptideHit& ph, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const))

  annot.annotateMatches(spec, hit, tg, sa);
  string annotlist[] = {"y1+", "y2+", "b2+", "y3+", "b3+", "y4+", "b4+", "y5+", "b5+", "b6+", "y6+"};

  PeakSpectrum::StringDataArray types = spec.getStringDataArrays().front();

  ABORT_IF(spec.size() != types.size() || types.size() != pls)
  for (size_t i = 0; i < spec.size(); ++i)
  {
    TEST_STRING_EQUAL(types[i],annotlist[i])
  }
  TEST_REAL_SIMILAR(spec.getMetaValue("fragment_mass_tolerance"),0.1)
END_SECTION

START_SECTION((void SpectrumAnnotator::addIonMatchStatistics(PeptideIdentification& pi, const MSSpectrum& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa) const))

  annot.addIonMatchStatistics(pi,spec,tg,sa);
  for (size_t i = 0; i < pi.getHits().size(); ++i)
  {
    TEST_EQUAL(pi.getHits()[i].getMetaValue("peak_number"),11)
    TEST_EQUAL(pi.getHits()[i].getMetaValue("sum_intensity"),12.1)
    TEST_EQUAL(pi.getHits()[i].getMetaValue("matched_ion_number"),11)
    TEST_EQUAL(pi.getHits()[i].getMetaValue("matched_intensity"),12.1)
    TEST_STRING_EQUAL(pi.getHits()[i].getMetaValue("matched_ions"),"y1+,y2+,b2+,y3+,b3+,y4+,b4+,y5+,b5+,b6+,y6+")
    TEST_STRING_EQUAL(pi.getHits()[i].getMetaValue("max_series_type"),"y")
    TEST_EQUAL(pi.getHits()[i].getMetaValue("max_series_size"),6)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("sn_by_matched_intensity"),0)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("sn_by_median_intensity"),0)
    TEST_EQUAL(pi.getHits()[i].getMetaValue("precursor_in_ms2"), false)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("topN_meanfragmenterror"),0.00051117)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("topN_MSEfragmenterror"),0)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("topN_stddevfragmenterror"),0.0002534)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("median_fragment_error"),0.0003167)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("IQR_fragment_error"),0.000486)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("NTermIonCurrentRatio"),0.454545)
    TEST_REAL_SIMILAR(pi.getHits()[i].getMetaValue("CTermIonCurrentRatio"),0.545454)
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



