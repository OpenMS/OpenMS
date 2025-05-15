// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/QC/PSMExplainedIonCurrent.h>
#include <cmath>
#include <random>
//////////////////////////

using namespace OpenMS;

// Functions to create input data

std::mt19937 gen(0);

void addRandomPeaks(MSSpectrum& spec, double total_intensity = 1.0, Int number_of_peaks = 1)
{
  double peak_intensity = total_intensity / number_of_peaks;
  spec.sortByPosition();
  std::uniform_real_distribution<> distr((*spec.begin()).getMZ(), (*(spec.end() - 1)).getMZ());
  for (Int i = 0; i < number_of_peaks; ++i)
  {
    spec.emplace_back(distr(gen), peak_intensity);
  }
}

// create a MSSpectrum with Precursor, MSLevel and RT without Peaks
const MSSpectrum createMSSpectrum(UInt ms_level, double rt, const String& id, Precursor::ActivationMethod precursor_method = Precursor::ActivationMethod::CID)
{
  Precursor precursor;
  std::set<Precursor::ActivationMethod> am;
  am.insert(precursor_method);
  precursor.setActivationMethods(am);

  MSSpectrum ms_spec;
  ms_spec.setRT(rt);
  ms_spec.setMSLevel(ms_level);
  ms_spec.setPrecursors({precursor});
  ms_spec.setNativeID(id);

  return ms_spec;
}

// create a MSSpectrum with Precursor, MSLevel, RT and fill it with Peaks
const MSSpectrum createMSSpectrum(UInt ms_level, double rt, const String& id, const AASequence& seq, Int charge, const Param& theo_gen_params,
                                  Precursor::ActivationMethod precursor_method = Precursor::ActivationMethod::CID)
{
  MSSpectrum ms_spec;
  ms_spec.setRT(rt);
  ms_spec.setMSLevel(ms_level);
  ms_spec.setNativeID(id);

  TheoreticalSpectrumGenerator t;
  if (!theo_gen_params.empty())
    t.setParameters(theo_gen_params);

  t.getSpectrum(ms_spec, seq, 1, charge <= 2 ? 1 : 2);
  std::set<Precursor::ActivationMethod> am;
  am.insert(precursor_method);
  ms_spec.getPrecursors()[0].setActivationMethods(am);

  return ms_spec;
}

// create a PeptideIdentifiaction with a PeptideHit (sequence, charge), rt and mz
// default values for sequence PEPTIDE
const PeptideIdentification createPeptideIdentification(const String& id, const String& sequence = String("PEPTIDE"), Int charge = 3, double mz = 266)
{
  PeptideHit peptide_hit;
  peptide_hit.setSequence(AASequence::fromString(sequence));
  peptide_hit.setCharge(charge);

  PeptideIdentification peptide_id;
  peptide_id.setSpectrumReference( id);
  peptide_id.setMZ(mz);
  peptide_id.setHits({peptide_hit});

  return peptide_id;
}

START_TEST(PSMExplainedIonCurrent, "$Id$")

AASequence::fromString("").empty();

// Generate test data

// MSExperiment
Param p;

// create b- and y-ion spectrum of peptide sequence HIMALAYA with charge 1
PeakSpectrum ms_spec_2_himalaya = createMSSpectrum(2, 3.7, "XTandem::1", AASequence::fromString("HIMALAYA"), 1, p);
addRandomPeaks(ms_spec_2_himalaya, 7.0); // add 7 to 13 -> correctness should be 13/20

// create c- and z-ion spectrum of peptide sequence ALABAMA with charge 2
TheoreticalSpectrumGenerator theo_gen_al;
p = theo_gen_al.getParameters();
p.setValue("add_c_ions", "true");
p.setValue("add_zp1_ions", "true");
p.setValue("add_b_ions", "false");
p.setValue("add_y_ions", "false");
PeakSpectrum ms_spec_2_alabama = createMSSpectrum(2, 2, "XTandem::2", AASequence::fromString("ALABAMA"), 2, p, Precursor::ActivationMethod::ECD);
addRandomPeaks(ms_spec_2_alabama, 5.0); // add 5 to 10 -> correctness should be 2/3

MSSpectrum empty_spec;

MSExperiment exp;
exp.setSpectra({empty_spec, ms_spec_2_alabama, ms_spec_2_himalaya});

// MSExperiment with no given fragmentation method (falls back to CID)
MSExperiment exp_no_pc(exp);
exp_no_pc[0].setPrecursors({});

// MSExperiment with MS1 Spectrum
MSExperiment exp_ms1(exp);
exp_ms1.setSpectra({createMSSpectrum(1, 5, "XTandem::3")});

// MSExperiment with Sori activation
MSExperiment exp_sori(exp);
exp_sori.setSpectra({createMSSpectrum(2, 7, "XTandem::5", Precursor::ActivationMethod::SORI)});

// MSExperiment with himalaya, spectrum with peaks with intensity 0 & empty spectrum
MSExperiment failing_exp(exp);
PeakSpectrum zero_peaks(createMSSpectrum(2, 4, "XTandem::6"));
zero_peaks.emplace_back(10, 0);
zero_peaks.emplace_back(20, 0);
failing_exp.setSpectra({ms_spec_2_himalaya, zero_peaks, empty_spec});

// map the MSExperiment
QCBase::SpectraMap spectra_map;

// PeptideIdentifications
PeptideIdentification empty_id;
empty_id.setRT(6);
PeptideIdentification himalaya = createPeptideIdentification("XTandem::1", "HIMALAYA", 1, 888);
PeptideIdentification alabama = createPeptideIdentification("XTandem::2", "ALABAMA", 2, 264);
PeptideIdentification no_hit_id(himalaya);
no_hit_id.setHits({});

std::vector<PeptideIdentification> pep_ids({himalaya, alabama, empty_id});

// ProteinIdentifications
ProteinIdentification protId;
ProteinIdentification::SearchParameters param;
param.fragment_mass_tolerance_ppm = false;
param.fragment_mass_tolerance = 0.3;
protId.setSearchParameters(param);

// FeatureMap
FeatureMap fmap;

Feature empty_feat;

fmap.setUnassignedPeptideIdentifications(pep_ids);
fmap.push_back(empty_feat);
fmap.setProteinIdentifications({protId});

PSMExplainedIonCurrent* ptr = nullptr;
PSMExplainedIonCurrent* nulpt = nullptr;
START_SECTION(PSMExplainedIonCurrent())
{
  ptr = new PSMExplainedIonCurrent();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION

START_SECTION(~PSMExplainedIonCurrent())
{
  delete ptr;
}
END_SECTION

PSMExplainedIonCurrent psm_corr;

START_SECTION(void compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20))
{
  spectra_map.calculateMap(exp);
  //--------------------------------------------------------------------
  // test with valid input - default parameter
  //--------------------------------------------------------------------
  psm_corr.compute(fmap, exp, spectra_map);
  std::vector<PSMExplainedIonCurrent::Statistics> result = psm_corr.getResults();

  TEST_REAL_SIMILAR(result[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with valid input - ToleranceUnit PPM
  //--------------------------------------------------------------------

  PSMExplainedIonCurrent psm_corr_ppm;
  psm_corr_ppm.compute(fmap, exp, spectra_map, QCBase::ToleranceUnit::PPM, 6);
  std::vector<PSMExplainedIonCurrent::Statistics> result_ppm = psm_corr_ppm.getResults();

  TEST_REAL_SIMILAR(result_ppm[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result_ppm[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with valid input and flags
  //--------------------------------------------------------------------
  PSMExplainedIonCurrent psm_corr_flag_da;
  psm_corr_flag_da.compute(fmap, exp, spectra_map, QCBase::ToleranceUnit::DA, 1);
  std::vector<PSMExplainedIonCurrent::Statistics> result_flag_da = psm_corr_flag_da.getResults();

  TEST_REAL_SIMILAR(result_flag_da[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result_flag_da[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with missing toleranceUnit and toleranceValue in featureMap
  //--------------------------------------------------------------------

  // featureMap with missing ProteinIdentifications
  {
    FeatureMap fmap_auto = fmap;
    fmap_auto.getProteinIdentifications().clear();
    TEST_EXCEPTION(Exception::MissingInformation, psm_corr.compute(fmap_auto, exp, spectra_map, QCBase::ToleranceUnit::AUTO))
  }

  //--------------------------------------------------------------------
  // test with no given fragmentation method
  //--------------------------------------------------------------------
  spectra_map.calculateMap(exp_no_pc);
  psm_corr.compute(fmap, exp_no_pc, spectra_map);
  TEST_REAL_SIMILAR(psm_corr.getResults()[1].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(psm_corr.getResults()[1].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with matching ms1 spectrum
  //--------------------------------------------------------------------
  {
    spectra_map.calculateMap(exp_ms1);
    // fmap with PeptideIdentification with spec_ref matching to a MS1 Spectrum
    FeatureMap fmap_ms1(fmap);
    fmap_ms1.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::3")});
    TEST_EXCEPTION(Exception::IllegalArgument, psm_corr.compute(fmap_ms1, exp_ms1, spectra_map))
  }

  //--------------------------------------------------------------------
  // test with fragmentation method SORI, which is not supported
  //--------------------------------------------------------------------
  {
    // put PeptideIdentification with spec_ref matching to MSSpectrum with fragmentation method SORI to fmap
    FeatureMap fmap_sori;
    fmap_sori.setProteinIdentifications({protId});
    fmap_sori.setUnassignedPeptideIdentifications({createPeptideIdentification("XTandem::5")});

    spectra_map.calculateMap(exp_sori);
    TEST_EXCEPTION(Exception::InvalidParameter, psm_corr.compute(fmap_sori, exp_sori, spectra_map))
  }

  //--------------------------------------------------------------------
  // Only failing inputs
  //--------------------------------------------------------------------
  {
    // put PeptideIdentification without hits & with spec_ref matching to spectrum with no peaks, empty spectrum to fmap
    FeatureMap failing_fmap;
    failing_fmap.setProteinIdentifications({protId});
    failing_fmap.setUnassignedPeptideIdentifications({no_hit_id, createPeptideIdentification("XTandem::6")});

    spectra_map.calculateMap(failing_exp);

    PSMExplainedIonCurrent psm_corr_failing;
    TEST_EXCEPTION(Exception::MissingInformation, psm_corr_failing.compute(failing_fmap, failing_exp, spectra_map))
  }
}
END_SECTION

START_SECTION(compute(std::vector<PeptideIdentification>& pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum,
                      ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20))
{
  spectra_map.calculateMap(exp);
  //--------------------------------------------------------------------
  // test with valid input - default parameter
  //--------------------------------------------------------------------
  psm_corr.compute(pep_ids, param, exp, spectra_map);
  std::vector<PSMExplainedIonCurrent::Statistics> result = psm_corr.getResults();

  TEST_REAL_SIMILAR(result[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with valid input - ToleranceUnit PPM
  //--------------------------------------------------------------------

  PSMExplainedIonCurrent psm_corr_ppm;
  psm_corr_ppm.compute(pep_ids, param, exp, spectra_map, QCBase::ToleranceUnit::PPM, 6);
  std::vector<PSMExplainedIonCurrent::Statistics> result_ppm = psm_corr_ppm.getResults();

  TEST_REAL_SIMILAR(result_ppm[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result_ppm[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with valid input and flags
  //--------------------------------------------------------------------
  PSMExplainedIonCurrent psm_corr_flag_da;
  psm_corr_flag_da.compute(pep_ids, param, exp, spectra_map, QCBase::ToleranceUnit::DA, 1);
  std::vector<PSMExplainedIonCurrent::Statistics> result_flag_da = psm_corr_flag_da.getResults();

  TEST_REAL_SIMILAR(result_flag_da[0].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(result_flag_da[0].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with missing toleranceUnit and toleranceValue from params
  //--------------------------------------------------------------------

  // featureMap with missing ProteinIdentifications
  {
    ProteinIdentification::SearchParameters no_params;
    TEST_EXCEPTION(Exception::MissingInformation, psm_corr.compute(pep_ids, no_params, exp, spectra_map, QCBase::ToleranceUnit::AUTO))
  }

  //--------------------------------------------------------------------
  // test with no given fragmentation method
  //--------------------------------------------------------------------
  spectra_map.calculateMap(exp_no_pc);
  psm_corr.compute(pep_ids, param, exp_no_pc, spectra_map);
  TEST_REAL_SIMILAR(psm_corr.getResults()[1].average_correctness, (13. / 20 + 10. / 15) / 2.)
  TEST_REAL_SIMILAR(psm_corr.getResults()[1].variance_correctness, 0.000138)

  //--------------------------------------------------------------------
  // test with matching ms1 spectrum
  //--------------------------------------------------------------------
  {
    spectra_map.calculateMap(exp_ms1);
    // PeptideIdentification with spec_ref matching to a MS1 Spectrum
    std::vector<PeptideIdentification> ms1_pep({createPeptideIdentification("XTandem::3")});
    TEST_EXCEPTION(Exception::IllegalArgument, psm_corr.compute(ms1_pep, param, exp_ms1, spectra_map))
  }

  //--------------------------------------------------------------------
  // test with fragmentation method SORI, which is not supported
  //--------------------------------------------------------------------
  {
    spectra_map.calculateMap(exp_sori);
    // PeptideIdentification with spec_ref matching to MSSpectrum with fragmentation method SORI
    std::vector<PeptideIdentification> sori_id({createPeptideIdentification("XTandem::5")});

    TEST_EXCEPTION(Exception::InvalidParameter, psm_corr.compute(sori_id, param, exp_sori, spectra_map))
  }

  //--------------------------------------------------------------------
  // Only failing inputs
  //--------------------------------------------------------------------
  {
    // put PeptideIdentification with spec_ref matching to MSSpectrum with no peaks to fmap
    std::vector<PeptideIdentification> failing_ids({no_hit_id, createPeptideIdentification("XTandem::6")});

    spectra_map.calculateMap(failing_exp);

    PSMExplainedIonCurrent psm_corr_failing;
    TEST_EXCEPTION(Exception::MissingInformation, psm_corr_failing.compute(failing_ids, param, failing_exp, spectra_map))
  }
}
END_SECTION

START_SECTION(const String& getName() const override)
{
  TEST_EQUAL(psm_corr.getName(), "PSMExplainedIonCurrent");
}
END_SECTION

START_SECTION(const std::vector<Statistics>& getResults() const)
{
  // tested in compute tests above
  NOT_TESTABLE;
}
END_SECTION

START_SECTION(QCBase::Status requirements() const override)
{
  QCBase::Status stat = QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  TEST_EQUAL(psm_corr.requirements() == stat, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
