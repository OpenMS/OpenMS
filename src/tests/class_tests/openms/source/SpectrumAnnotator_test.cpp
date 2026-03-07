// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <algorithm>
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
    TEST_STRING_EQUAL(types[i], annotlist[i])
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

START_SECTION((void SpectrumAnnotator::addPeakAnnotationsToPeptideHit(PeptideHit& ph, const PeakSpectrum& spec, const TheoreticalSpectrumGenerator& tg, const SpectrumAlignment& sa, bool include_unmatched_peaks) const))
  // Create a fresh PeptideHit for testing
  PeptideHit test_hit;
  test_hit.setSequence(peptide);
  test_hit.setCharge(2);

  // Create a fresh spectrum (copy of the test data)
  PeakSpectrum test_spec;
  test_spec.setMSLevel(2);
  Peak1D test_p;
  for (Size i = 0; i != pls; ++i)
  {
    test_p.setIntensity(1.1f);
    test_p.setMZ(peaklist[i]);
    test_spec.push_back(test_p);
  }

  // Test with include_unmatched_peaks = false (default behavior)
  annot.addPeakAnnotationsToPeptideHit(test_hit, test_spec, tg, sa, false);

  // Verify the PeakAnnotations were set correctly
  const std::vector<PeptideHit::PeakAnnotation>& pas = test_hit.getPeakAnnotations();
  
  // We expect 11 annotations (same as the number of matched peaks)
  TEST_EQUAL(pas.size(), 11)
  
  // Verify the annotations contain the expected ion names
  StringList expected_ions = ListUtils::create<String>("y1+,y2+,b2+,y3+,b3+,y4+,b4+,y5+,b5+,b6+,y6+");
  StringList found_ions;
  for (const auto& pa : pas)
  {
    found_ions.push_back(pa.annotation);
    // Verify that mz and intensity are set
    TEST_NOT_EQUAL(pa.mz, -1.0)
    TEST_REAL_SIMILAR(pa.intensity, 1.1f)
  }
  
  // Sort both lists for comparison (since order may vary based on spectrum alignment)
  std::sort(expected_ions.begin(), expected_ions.end());
  std::sort(found_ions.begin(), found_ions.end());
  TEST_EQUAL(found_ions.size(), expected_ions.size())
  for (size_t i = 0; i < found_ions.size(); ++i)
  {
    TEST_STRING_EQUAL(found_ions[i], expected_ions[i])
  }

  // Test with include_unmatched_peaks = true
  // Add some unmatched peaks to the spectrum
  PeakSpectrum test_spec_with_unmatched;
  test_spec_with_unmatched.setMSLevel(2);
  // Add some peaks that will NOT match
  Peak1D unmatched_p;
  unmatched_p.setIntensity(0.5f);
  unmatched_p.setMZ(100.0); // This m/z won't match any theoretical ion
  test_spec_with_unmatched.push_back(unmatched_p);
  unmatched_p.setMZ(1000.0); // This m/z also won't match
  test_spec_with_unmatched.push_back(unmatched_p);
  // Add the original matched peaks
  for (Size i = 0; i != pls; ++i)
  {
    test_p.setIntensity(1.1f);
    test_p.setMZ(peaklist[i]);
    test_spec_with_unmatched.push_back(test_p);
  }

  PeptideHit test_hit2;
  test_hit2.setSequence(peptide);
  test_hit2.setCharge(2);
  
  annot.addPeakAnnotationsToPeptideHit(test_hit2, test_spec_with_unmatched, tg, sa, true);
  
  const std::vector<PeptideHit::PeakAnnotation>& pas2 = test_hit2.getPeakAnnotations();
  
  // We expect 13 annotations (11 matched + 2 unmatched)
  TEST_EQUAL(pas2.size(), 13)
  
  // Count annotated and unannotated peaks
  size_t annotated_count = 0;
  size_t unannotated_count = 0;
  for (const auto& pa : pas2)
  {
    if (pa.annotation.empty())
    {
      unannotated_count++;
    }
    else
    {
      annotated_count++;
    }
    // All peaks should have mz and intensity set
    TEST_NOT_EQUAL(pa.mz, -1.0)
  }
  
  TEST_EQUAL(annotated_count, 11)  // 11 matched peaks
  TEST_EQUAL(unannotated_count, 2) // 2 unmatched peaks
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



