// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Mohammed Alhigaylan$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
///////////////////////////

#include <cmath>

using namespace OpenMS;
using namespace std;

namespace
{
  /// A mass trace with a Gaussian elution profile, sampled once per second and
  /// prepared the way ElutionPeakDetection hands traces to FeatureFindingMetabo.
  MassTrace makeGaussianTrace(double mz, double apex_rt, double height, double sigma, const std::string& label)
  {
    std::vector<Peak2D> peaks;
    std::vector<double> intensities;
    for (int i = -10; i <= 10; ++i)
    {
      Peak2D p;
      p.setRT(apex_rt + i);
      p.setMZ(mz);
      p.setIntensity(height * std::exp(-0.5 * i * i / (sigma * sigma)));
      peaks.push_back(p);
      intensities.push_back(p.getIntensity());
    }
    MassTrace mt(peaks);
    mt.setLabel(label);
    mt.updateWeightedMeanRT();
    mt.updateWeightedMeanMZ();
    mt.setCentroidSD(0.001);
    mt.setSmoothedIntensities(intensities);
    mt.estimateFWHM(true);
    return mt;
  }

  /// Assemble the given traces (FeatureFindingMetabo sorts its input, hence the copy).
  FeatureMap assemble(std::vector<MassTrace> traces, const Param& changed_params)
  {
    FeatureFindingMetabo ffm;
    Param p = ffm.getParameters();
    for (Param::ParamIterator it = changed_params.begin(); it != changed_params.end(); ++it)
    {
      p.setValue(it.getName(), it->value);
    }
    ffm.setParameters(p);
    FeatureMap fm;
    std::vector<std::vector<MSChromatogram>> chromatograms;
    ffm.run(traces, fm, chromatograms);
    return fm;
  }

  /// Run the mass trace pipeline with the settings of the FeatureFinderMetabo TOPP tool.
  FeatureMap runPipeline(const std::string& mzml_file, const Param& changed_ffm_params)
  {
    PeakMap exp;
    MzMLFile().load(mzml_file, exp);
    exp.sortSpectra(true);

    MassTraceDetection mtd;
    Param p = mtd.getParameters();
    p.setValue("noise_threshold_int", 10.0);
    p.setValue("chrom_peak_snr", 3.0);
    mtd.setParameters(p);
    std::vector<MassTrace> traces, split_traces;
    mtd.run(exp, traces);

    ElutionPeakDetection epd;
    p = epd.getParameters();
    p.setValue("chrom_peak_snr", 3.0);
    p.setValue("chrom_fwhm", 5.0);
    p.setValue("width_filtering", "fixed");
    epd.setParameters(p);
    epd.detectPeaks(traces, split_traces);

    return assemble(split_traces, changed_ffm_params);
  }

  /// The feature whose monoisotopic trace has the given m/z (or nullptr).
  const Feature* findFeature(const FeatureMap& fm, double mz, double mz_tolerance = 0.01)
  {
    for (const Feature& f : fm)
    {
      if (std::fabs(f.getMZ() - mz) <= mz_tolerance)
      {
        return &f;
      }
    }
    return nullptr;
  }
}

START_TEST(FeatureFindingMetabo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFindingMetabo* ptr = nullptr;
FeatureFindingMetabo* null_ptr = nullptr;
START_SECTION(FeatureFindingMetabo())
{
  ptr = new FeatureFindingMetabo();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFindingMetabo())
{
  delete ptr;
}
END_SECTION

// load a mzML file for testing the algorithm
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_input1.mzML"), input);

FeatureMap test_fm;

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

std::vector<std::vector< OpenMS::MSChromatogram > > chromatograms;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
test_epd.detectPeaks(output_mt, splitted_mt);

FuzzyStringComparator fsc;
fsc.setAcceptableRelative(1.001);
fsc.setAcceptableAbsolute(1);
StringList sl;
sl.push_back("xml-stylesheet");
sl.push_back("<featureMap");
sl.push_back("<feature id");
fsc.setWhitelist(sl);

//std::cout << "\n\n" << fsc.compareStrings("529090", "529091") << "\n\n\n";

START_SECTION((void run(std::vector< MassTrace > &, FeatureMap &, chromatograms &)))
{
  FeatureFindingMetabo test_ffm;
  // run with non-default setting (C13 isotope distance)
  Param p = test_ffm.getParameters();
  p.setValue("mz_scoring_13C", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 79);

  // run with default settings (from paper using charge+isotope# dependent distances)
  p.setValue("report_convex_hulls", "true");
  p.setValue("mz_scoring_13C", "false");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 76);
  // --> this gives less features, i.e. more isotope clusters (but the input data is simulated and highly weird -- should be replaced at some point)

  //todo the new isotope m/z scoring should produce similar results, but still has to be tested.
  p.setValue("mz_scoring_by_elements", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 75);
  // --> this gives less features, i.e. more isotope clusters (but the input data is simulated and highly weird -- should be replaced at some point)
}
END_SECTION

START_SECTION(([EXTRA] run with historic settings))
{
  // The RT overlap measured against the longer trace and missing isotope traces
  // counted as zero intensity (the behavior before issue #1942) assemble fewer
  // traces: the features of each run below are 4 to 5 more than with the defaults.
  FeatureFindingMetabo test_ffm;
  Param p = test_ffm.getParameters();
  p.setValue("isotope_rt_overlap_reference", "longer");
  p.setValue("isotope_filtering_missing_traces", "zero");

  // C13 isotope distance
  p.setValue("mz_scoring_13C", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  // one feature less than before the MassTrace::estimateFWHM right-flank fix (#10052): the RT
  // co-elution score divides the FWHM-border overlap by the trace FWHM, so a corrected FWHM moves
  // hypotheses across 'min_isotope_rt_overlap' and one more trace pair is assembled as isotopes
  TEST_EQUAL(test_fm.size(), 83);

  // charge+isotope# dependent distances (from paper)
  p.setValue("report_convex_hulls", "true");
  p.setValue("mz_scoring_13C", "false");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 81);

  // test annotation of input
  std::string tmp_file;
  NEW_TMP_FILE(tmp_file);
  FeatureXMLFile().store(tmp_file, test_fm);
  TEST_EQUAL(fsc.compareFiles(tmp_file, OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_output1.featureXML")), true);

  p.setValue("mz_scoring_by_elements", "true");
  test_ffm.setParameters(p);
  test_ffm.run(splitted_mt, test_fm, chromatograms);
  TEST_EQUAL(test_fm.size(), 80);
}
END_SECTION

// Test with ion mobility data
START_SECTION(([EXTRA] run with ion mobility data))
{
  // Load LC-IMS-MS data
  PeakMap im_input;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("../../../topp/FeatureFinderMetabo_6_input.mzML"), im_input);

  // Run MassTraceDetection with ion mobility tolerance
  MassTraceDetection mtd_im;
  Param mtd_p = mtd_im.getParameters();
  mtd_p.setValue("mass_error_ppm", 10.0);
  mtd_p.setValue("noise_threshold_int", 10.0);
  mtd_p.setValue("ion_mobility_tolerance", 0.01);
  mtd_im.setParameters(mtd_p);
  std::vector<MassTrace> im_traces;
  mtd_im.run(im_input, im_traces);
  TEST_EQUAL(im_traces.empty(), false);

  // Verify mass traces contain IM data
  bool has_im = false;
  for (const auto& mt : im_traces)
  {
    if (mt.containsIMData())
    {
      has_im = true;
      break;
    }
  }
  TEST_EQUAL(has_im, true);

  // Run ElutionPeakDetection
  ElutionPeakDetection epd_im;
  std::vector<MassTrace> im_splitted;
  epd_im.detectPeaks(im_traces, im_splitted);
  TEST_EQUAL(im_splitted.empty(), false);

  // Verify IM data survives elution peak detection splitting
  has_im = false;
  for (const auto& mt : im_splitted)
  {
    if (mt.containsIMData())
    {
      has_im = true;
      break;
    }
  }
  TEST_EQUAL(has_im, true);

  // Run FeatureFindingMetabo with IM-aware settings
  FeatureFindingMetabo ffm_im;
  Param ffm_p = ffm_im.getParameters();
  ffm_p.setValue("mz_scoring_13C", "true");
  ffm_im.setParameters(ffm_p);
  FeatureMap im_fm;
  std::vector<std::vector<OpenMS::MSChromatogram>> im_chromatograms;
  ffm_im.run(im_splitted, im_fm, im_chromatograms);
  TEST_EQUAL(im_fm.empty(), false);

  // Verify that output features contain ion mobility metadata
  bool has_im_meta = false;
  for (const auto& f : im_fm)
  {
    if (f.metaValueExists("masstrace_centroid_im"))
    {
      has_im_meta = true;
      break;
    }
  }
  TEST_EQUAL(has_im_meta, true);
}
END_SECTION

START_SECTION(([EXTRA] defaults))
{
  Param p = FeatureFindingMetabo().getParameters();
  TEST_EQUAL(p.getValue("isotope_rt_overlap_reference"), "shorter")
  TEST_EQUAL(p.getValue("isotope_filtering_missing_traces"), "unknown")
}
END_SECTION

START_SECTION(([EXTRA] isotope_rt_overlap_reference))
{
  // An isotope trace of low intensity is shorter than the monoisotopic trace, because
  // only its most intense part exceeds the noise: the two overlap in 4 s here, which
  // is all of the isotope trace (FWHM 3.5 s), but less than half of the monoisotopic
  // trace (FWHM 9.4 s).
  std::vector<MassTrace> traces;
  traces.push_back(makeGaussianTrace(300.1, 100.0, 1.0e5, 4.0, "mono"));
  traces.push_back(makeGaussianTrace(301.102, 100.0, 2.0e4, 1.5, "iso1"));

  FeatureMap fm = assemble(traces, Param());
  TEST_EQUAL(fm.size(), 1)
  TEST_EQUAL(fm[0].getCharge(), 1)
  TEST_EQUAL(fm[0].getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES), 2)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 300.1)

  Param longer;
  longer.setValue("isotope_rt_overlap_reference", "longer");
  fm = assemble(traces, longer);
  TEST_EQUAL(fm.size(), 2)
  TEST_EQUAL(fm[0].getCharge(), 0)
  TEST_EQUAL(fm[1].getCharge(), 0)

  // the apex of the longer trace has to lie within the shorter trace
  traces[1] = makeGaussianTrace(301.102, 105.0, 2.0e4, 1.5, "iso1");
  fm = assemble(traces, Param());
  TEST_EQUAL(fm.size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] isotope_filtering_missing_traces))
{
  // A molecule of 850 Da whose second and third isotopic trace are below the
  // detection limit. The first isotopic trace has 70% of the monoisotopic intensity,
  // which is common at this mass (about 60 carbon atoms).
  std::vector<MassTrace> traces;
  traces.push_back(makeGaussianTrace(850.5, 100.0, 1.0e5, 3.0, "mono"));
  traces.push_back(makeGaussianTrace(851.502, 100.0, 0.7e5, 3.0, "iso1"));

  // unknown intensity of the missing traces: accepted
  FeatureMap fm = assemble(traces, Param());
  TEST_EQUAL(fm.size(), 1)
  TEST_EQUAL(fm[0].getCharge(), 1)
  TEST_EQUAL(fm[0].getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES), 2)
  TEST_EQUAL(fm[0].getMetaValue("legal_isotope_pattern"), 1)

  // missing traces counted as zero intensity: the model rejects the pattern
  Param zero;
  zero.setValue("isotope_filtering_missing_traces", "zero");
  fm = assemble(traces, zero);
  TEST_EQUAL(fm.size(), 2)
  TEST_EQUAL(fm[0].getCharge(), 0)
  TEST_EQUAL(fm[1].getCharge(), 0)

  // a first isotopic trace that fits the mass on its own is accepted either way
  traces[1] = makeGaussianTrace(851.502, 100.0, 0.45e5, 3.0, "iso1");
  TEST_EQUAL(assemble(traces, Param()).size(), 1)
  TEST_EQUAL(assemble(traces, zero).size(), 1)

  // a first isotopic trace that is more intense than the monoisotopic trace is
  // impossible at this mass, whatever the intensity of the missing traces is
  traces[1] = makeGaussianTrace(851.502, 100.0, 1.3e5, 3.0, "iso1");
  TEST_EQUAL(assemble(traces, Param()).size(), 2)
  TEST_EQUAL(assemble(traces, zero).size(), 2)

  // the setting does not apply without a metabolite isotope model
  Param none = zero;
  none.setValue("isotope_filtering_model", "none");
  TEST_EQUAL(assemble(traces, none).size(), 1)
}
END_SECTION

START_SECTION(([EXTRA] isotope traces of low intensity (issue #1942)))
{
  // Excerpts (m/z 518-529, RT 160-240 s) of two runs of the same study, taken from
  // issue #1942: a doubly charged ion at m/z 523.286 elutes at 200 s with isotope
  // traces at m/z 523.786 and 524.288. Its traces are 3 to 5 scans wide.
  const std::string low = OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_issue1942_low.mzML");
  const std::string high = OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_issue1942_high.mzML");

  Param longer, zero, historic;
  longer.setValue("isotope_rt_overlap_reference", "longer");
  zero.setValue("isotope_filtering_missing_traces", "zero");
  historic.setValue("isotope_rt_overlap_reference", "longer");
  historic.setValue("isotope_filtering_missing_traces", "zero");

  // the ion is found in both runs
  FeatureMap fm = runPipeline(low, Param());
  const Feature* f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 2)
  TEST_EQUAL(f->getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES), 3)

  fm = runPipeline(high, Param());
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 2)
  // the trace at m/z 524.288 also holds a co-eluting ion that starts there
  TEST_EQUAL(f->getMetaValue(Constants::UserParam::NUM_OF_MASSTRACES), 2)

  // the historic settings leave the monoisotopic trace on its own in both runs
  fm = runPipeline(low, historic);
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 0)

  fm = runPipeline(high, historic);
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 0)

  // run of low intensity: measuring the RT overlap against the shorter trace suffices
  fm = runPipeline(low, zero);
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 2)

  // run of high intensity: needs both settings
  fm = runPipeline(high, zero);
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 0)

  fm = runPipeline(high, longer);
  f = findFeature(fm, 523.286);
  TEST_NOT_EQUAL(f, nullptr)
  ABORT_IF(f == nullptr)
  TEST_EQUAL(f->getCharge(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST



