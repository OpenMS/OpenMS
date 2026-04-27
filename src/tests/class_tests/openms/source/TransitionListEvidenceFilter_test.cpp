// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace OpenMS;
using namespace OpenSwath;
using namespace std;

namespace
{
  ChromExtractParams makeExtractParams(double mz_window)
  {
    ChromExtractParams params;
    params.min_upper_edge_dist = 0.0;
    params.mz_extraction_window = mz_window;
    params.im_extraction_window = -1.0;
    params.ppm = false;
    params.extraction_function = "tophat";
    params.rt_extraction_window = -1.0;
    params.extra_rt_extract = 0.0;
    return params;
  }

  MSSpectrum makeSpectrum(double rt, const vector<pair<double, double>>& peaks)
  {
    MSSpectrum spectrum;
    spectrum.setRT(rt);
    for (const auto& peak_data : peaks)
    {
      Peak1D peak;
      peak.setMZ(peak_data.first);
      peak.setIntensity(peak_data.second);
      spectrum.push_back(peak);
    }
    spectrum.sortByPosition();
    return spectrum;
  }

  void appendProfilePeak(vector<pair<double, double>>& peaks, double center_mz, double apex_intensity)
  {
    const vector<pair<double, double>> shape{
      {-0.020, 0.005},
      {-0.015, 0.030},
      {-0.010, 0.120},
      {-0.005, 0.450},
      { 0.000, 1.000},
      { 0.005, 0.450},
      { 0.010, 0.120},
      { 0.015, 0.030},
      { 0.020, 0.005}
    };

    for (const auto& point : shape)
    {
      peaks.emplace_back(center_mz + point.first, apex_intensity * point.second);
    }
  }

  MSSpectrum makeIMSpectrum(double rt, const vector<tuple<double, double, double>>& peaks)
  {
    MSSpectrum spectrum;
    spectrum.setRT(rt);
    MSSpectrum::FloatDataArray im_array;
    im_array.setName("Ion Mobility");
    for (const auto& peak_data : peaks)
    {
      Peak1D peak;
      peak.setMZ(get<0>(peak_data));
      peak.setIntensity(get<1>(peak_data));
      spectrum.push_back(peak);
      im_array.push_back(get<2>(peak_data));
    }
    spectrum.getFloatDataArrays().push_back(im_array);
    return spectrum;
  }

  SwathMap makeSwathMap(bool ms1,
                        double lower,
                        double upper,
                        const vector<MSSpectrum>& spectra,
                        double im_lower = -1.0,
                        double im_upper = -1.0)
  {
    PeakMap* peak_map = new PeakMap;
    for (const auto& spectrum : spectra)
    {
      peak_map->addSpectrum(spectrum);
    }
    shared_ptr<PeakMap> exp(peak_map);

    SwathMap map = im_lower >= 0.0 ?
      SwathMap(lower, upper, (lower + upper) / 2.0, im_lower, im_upper, ms1) :
      SwathMap(lower, upper, (lower + upper) / 2.0, ms1);
    map.sptr = SpectrumAccessPtr(new SpectrumAccessOpenMS(exp));
    return map;
  }

  void addCompound(LightTargetedExperiment& exp,
                   const string& id,
                   double precursor_mz,
                   const vector<double>& fragment_mzs,
                   bool decoy = false,
                   double precursor_im = -1.0)
  {
    LightCompound compound;
    compound.id = id;
    compound.sequence = id;
    compound.rt = 100.0;
    compound.drift_time = precursor_im;
    compound.charge = 2;
    compound.protein_refs.push_back("protein");
    exp.compounds.push_back(compound);

    for (Size i = 0; i < fragment_mzs.size(); ++i)
    {
      LightTransition transition;
      transition.transition_name = id + "_tr" + std::to_string(i);
      transition.peptide_ref = id;
      transition.precursor_mz = precursor_mz;
      transition.precursor_im = precursor_im;
      transition.product_mz = fragment_mzs[i];
      transition.library_intensity = 1000.0 - static_cast<double>(i);
      transition.setDetectingTransition(true);
      transition.setQuantifyingTransition(true);
      transition.setDecoy(decoy);
      exp.transitions.push_back(transition);
    }
  }

  LightTargetedExperiment makeLibrary()
  {
    LightTargetedExperiment exp;
    LightProtein protein;
    protein.id = "protein";
    exp.proteins.push_back(protein);

    addCompound(exp, "PEP_A", 500.0, {100.0, 110.0, 120.0});
    addCompound(exp, "PEP_B", 600.0, {200.0, 210.0, 220.0});
    addCompound(exp, "DECOY_PEP", 700.0, {300.0, 310.0, 320.0}, true);
    return exp;
  }

  TransitionListEvidenceFilter makeFilter(const String& evidence_sources, Size min_supported = 1)
  {
    TransitionListEvidenceFilter filter;
    Param params = filter.getParameters();
    params.setValue("enabled", "true");
    params.setValue("evidence_sources", evidence_sources);
    params.setValue("min_supported_precursors", static_cast<Int>(min_supported));
    params.setValue("ms2_min_fragment_hits", 2);
    params.setValue("ms2_top_transitions_per_precursor", 3);
    params.setValue("ms1_top_peaks_per_spectrum", 10);
    params.setValue("ms2_top_peaks_per_spectrum", 10);
    filter.setParameters(params);
    return filter;
  }
}

START_TEST(TransitionListEvidenceFilter, "$Id$")

START_SECTION(TransitionListEvidenceFilter())
{
  TransitionListEvidenceFilter filter;
  TEST_EQUAL(filter.getParameters().getValue("enabled").toString(), "false")
  TEST_EQUAL(filter.getParameters().getValue("evidence_sources").toString(), "hybrid")
}
END_SECTION

START_SECTION((filter() - MS1 support))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{500.004, 1000.0}, {700.0, 5000.0}})}));

  TransitionListEvidenceFilter filter = makeFilter("ms1");
  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 1);

  TEST_EQUAL(result.total_target_precursors, 2)
  TEST_EQUAL(result.supported_precursors, 1)
  TEST_EQUAL(result.ms1_supported, 1)
  TEST_EQUAL(result.ms2_supported, 0)
  TEST_EQUAL(result.filtered_targets.compounds.size(), 1)
  TEST_EQUAL(result.filtered_targets.compounds[0].id, "PEP_A")
  TEST_EQUAL(result.filtered_targets.transitions.size(), 3)
}
END_SECTION

START_SECTION((filter() - MS2 support))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0, {makeSpectrum(12.0, {{100.002, 1200.0}, {110.003, 900.0}, {200.0, 500.0}})}));

  TransitionListEvidenceFilter filter = makeFilter("ms2");
  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 1);

  TEST_EQUAL(result.total_target_precursors, 2)
  TEST_EQUAL(result.supported_precursors, 1)
  TEST_EQUAL(result.ms1_supported, 0)
  TEST_EQUAL(result.ms2_supported, 1)
  TEST_EQUAL(result.filtered_targets.compounds.size(), 1)
  TEST_EQUAL(result.filtered_targets.compounds[0].id, "PEP_A")
  TEST_EQUAL(result.evidence[0].ms2_best_fragment_hits, 2)
}
END_SECTION

START_SECTION((filter() - hybrid support and decoy exclusion))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{500.001, 1000.0}, {700.0, 9000.0}})}));
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0, {makeSpectrum(12.0, {{200.002, 1100.0}, {210.002, 1000.0}, {300.0, 8000.0}, {310.0, 8000.0}})}));

  TransitionListEvidenceFilter filter = makeFilter("hybrid", 2);
  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 1);

  TEST_EQUAL(result.total_target_precursors, 2)
  TEST_EQUAL(result.supported_precursors, 2)
  TEST_EQUAL(result.ms1_supported, 1)
  TEST_EQUAL(result.ms2_supported, 1)
  TEST_EQUAL(result.filtered_targets.compounds.size(), 2)
  TEST_EQUAL(result.filtered_targets.transitions.size(), 6)
  TEST_EQUAL(transition_exp.compounds.size(), 3)
  TEST_EQUAL(transition_exp.transitions.size(), 9)
  for (const auto& transition : result.filtered_targets.transitions)
  {
    TEST_EQUAL(transition.getDecoy(), false)
    TEST_NOT_EQUAL(transition.getPeptideRef(), "DECOY_PEP")
  }
}
END_SECTION

START_SECTION((filter() - peak picking path))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  vector<pair<double, double>> profile_peaks;
  appendProfilePeak(profile_peaks, 100.0, 2000.0);
  appendProfilePeak(profile_peaks, 110.0, 1500.0);
  appendProfilePeak(profile_peaks, 120.0, 1000.0);
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0,
                                    {makeSpectrum(12.0, profile_peaks)}));

  TransitionListEvidenceFilter filter = makeFilter("ms2");
  Param params = filter.getParameters();
  params.setValue("peak_picking:enabled", "true");
  params.setValue("peak_picking:PeakPickerHiRes:signal_to_noise", 0.0);
  filter.setParameters(params);

  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.08), false, 1);
  TEST_EQUAL(result.supported_precursors, 1)
  TEST_EQUAL(result.filtered_targets.compounds[0].id, "PEP_A")
}
END_SECTION

START_SECTION((filter() - PASEF precursor ion mobility auto scale))
{
  LightTargetedExperiment transition_exp;
  LightProtein protein;
  protein.id = "protein";
  transition_exp.proteins.push_back(protein);
  addCompound(transition_exp, "PEP_IM", 500.0, {100.0, 110.0, 120.0}, false, 0.0007);

  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0,
                                    {makeIMSpectrum(12.0, {{100.002, 1200.0, 0.70}, {110.002, 900.0, 0.70}})},
                                    0.60, 0.80));

  TransitionListEvidenceFilter filter = makeFilter("ms2");
  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), true, 1);

  TEST_EQUAL(result.supported_precursors, 1)
  TEST_REAL_SIMILAR(result.precursor_im_scale, 1000.0)
  TEST_REAL_SIMILAR(result.evidence[0].precursor_im, 0.7)
  TEST_REAL_SIMILAR(result.filtered_targets.compounds[0].getDriftTime(), 0.7)
  TEST_REAL_SIMILAR(result.filtered_targets.transitions[0].getPrecursorIM(), 0.7)
}
END_SECTION

START_SECTION((filter() - PASEF precursor ion mobility charge transform))
{
  LightTargetedExperiment transition_exp;
  LightProtein protein;
  protein.id = "protein";
  transition_exp.proteins.push_back(protein);
  addCompound(transition_exp, "PEP_IM_CHARGE", 500.0, {100.0, 110.0, 120.0}, false, 0.35);

  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0,
                                    {makeIMSpectrum(12.0, {{100.002, 1200.0, 0.70}, {110.002, 900.0, 0.70}})},
                                    0.60, 0.80));

  TransitionListEvidenceFilter filter = makeFilter("ms2");
  TransitionListEvidenceFilter::Result result = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), true, 1);

  TEST_EQUAL(result.supported_precursors, 1)
  TEST_REAL_SIMILAR(result.precursor_im_scale, 1.0)
  TEST_EQUAL(result.precursor_im_scaled_by_charge, true)
  TEST_REAL_SIMILAR(result.evidence[0].precursor_im, 0.7)
  TEST_REAL_SIMILAR(result.filtered_targets.compounds[0].getDriftTime(), 0.7)
  TEST_REAL_SIMILAR(result.filtered_targets.transitions[0].getPrecursorIM(), 0.7)
}
END_SECTION

START_SECTION((filter() - too few supported precursors))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0, {makeSpectrum(10.0, {{500.001, 1000.0}})}));

  TransitionListEvidenceFilter filter = makeFilter("ms1", 2);
  TEST_EXCEPTION(Exception::IllegalArgument, filter.filter(swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 1))
}
END_SECTION

START_SECTION((filter() - threaded consistency))
{
  LightTargetedExperiment transition_exp = makeLibrary();
  vector<SwathMap> swath_maps;
  swath_maps.push_back(makeSwathMap(true, 0.0, 0.0,
                                    {makeSpectrum(10.0, {{500.001, 1000.0}}),
                                     makeSpectrum(11.0, {{600.001, 1200.0}})}));
  swath_maps.push_back(makeSwathMap(false, 400.0, 650.0,
                                    {makeSpectrum(12.0, {{100.002, 1100.0}, {110.002, 1000.0}}),
                                     makeSpectrum(13.0, {{200.002, 1400.0}, {210.002, 1300.0}})}));

  TransitionListEvidenceFilter filter = makeFilter("hybrid", 2);
  TransitionListEvidenceFilter::Result single_thread = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 1);
  TransitionListEvidenceFilter::Result two_threads = filter.filter(
    swath_maps, transition_exp, makeExtractParams(0.02), makeExtractParams(0.02), false, 2);

  TEST_EQUAL(single_thread.supported_precursors, two_threads.supported_precursors)
  TEST_EQUAL(single_thread.ms1_supported, two_threads.ms1_supported)
  TEST_EQUAL(single_thread.ms2_supported, two_threads.ms2_supported)
  TEST_EQUAL(single_thread.filtered_targets.transitions.size(), two_threads.filtered_targets.transitions.size())
}
END_SECTION

END_TEST
