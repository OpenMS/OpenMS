// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <cmath>
#include <sstream>

using namespace OpenMS;
using namespace std;

namespace
{
  class TestableFeatureFinderIdentificationAlgorithm : public FeatureFinderIdentificationAlgorithm
  {
  public:
    void setPeptideCount(Size count) { n_peptides_ = count; }
    Size getPeptideCount() const { return n_peptides_; }
    void writeStatistics(const FeatureMap& features) const { statistics_(features); }
  };
}

START_TEST(FeatureFinderIdentificationAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderIdentificationAlgorithm* ptr = 0;
FeatureFinderIdentificationAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderIdentificationAlgorithm())
{
  ptr = new FeatureFinderIdentificationAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFinderIdentificationAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] summary statistics cannot wrap a missing-peptide count)
{
  // Exercise the reporting boundary with deliberately inconsistent counters.
  // This is the shape that previously made an unsigned subtraction print a
  // huge "peptides without features" value instead of a non-negative count.
  FeatureMap features;
  IdentificationData::Run& run = features.getIdentificationData().addRun("run");
  const auto source = run.addSource(IdentificationData::SourceFile());
  const auto query = run.addIdentification(source, IdentificationData::Observation());
  IdentificationData::MatchData match;
  match.representation = "PEPTIDEK";
  match.charge = 2;
  const auto match_id = run.addMatch(query, match);

  Feature feature;
  feature.setIntensity(1.0);
  feature.addIDMatch({run.getUuid(), match_id});
  features.push_back(feature);

  TestableFeatureFinderIdentificationAlgorithm algo;
  algo.setPeptideCount(0);

  ostringstream captured_info;
  OPENMS_LOG_INFO.insert(captured_info);
  algo.writeStatistics(features);
  OPENMS_LOG_INFO.remove(captured_info);

  TEST_TRUE(captured_info.str().find("0 peptides without features") != string::npos)
}
END_SECTION


START_SECTION([EXTRA] run() rejects MS data without spectra)
{
  // Chromatogram extraction indexes the experiment unconditionally, so running on an
  // empty experiment used to read out of bounds instead of reporting the problem
  // (https://github.com/OpenMS/OpenMS/issues/9980).
  FeatureFinderIdentificationAlgorithm algo;

  ProteinIdentification protein;
  protein.setIdentifier("run0");
  vector<ProteinIdentification> proteins{protein};

  PeptideIdentification peptide;
  peptide.setIdentifier("run0");
  peptide.setRT(100.0);
  peptide.setMZ(500.0);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  peptide.setHits({hit});
  PeptideIdentificationList peptides{peptide};

  FeatureMap features;
  FeatureMap seeds;

  TEST_EQUAL(algo.getMSData().empty(), true)
  TEST_EXCEPTION(Exception::IllegalArgument,
                 algo.run(peptides, proteins, features, seeds, ""))
  TEST_EXCEPTION(Exception::IllegalArgument,
                 algo.run(IdentificationData(), features, seeds, ""))
}
END_SECTION

START_SECTION([EXTRA] run() counts the distinct identified peptides for the summary)
{
  // The summary statistics report the peptides of the run, not a stale or uninitialised count: three distinct
  // peptides, one of them identified twice, eluting in a small MS1-only experiment.
  const std::vector<std::pair<std::string, double>> peptides_rt = {{"PEPTIDEK", 100.0}, {"PEPTIDEK", 101.0}, {"SAMPLER", 130.0},
                                                                   {"ELVISLIVESK", 160.0}};
  PeakMap ms1;
  for (int scan = 0; scan < 111; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(80.0 + scan);
    spectrum.setNativeID("scan=" + std::to_string(scan + 1));
    for (const auto& [sequence, apex] : std::vector<std::pair<std::string, double>>{{"PEPTIDEK", 100.0}, {"SAMPLER", 130.0}, {"ELVISLIVESK", 160.0}})
    {
      const double mz = AASequence::fromString(sequence).getMZ(2);
      const double height = 1e6 * std::exp(-0.5 * std::pow((spectrum.getRT() - apex) / 3.0, 2));
      for (int isotope = 0; isotope < 3; ++isotope)
      {
        spectrum.push_back(Peak1D(mz + isotope * Constants::C13C12_MASSDIFF_U / 2, height / (isotope + 1)));
      }
    }
    spectrum.sortByPosition();
    ms1.addSpectrum(spectrum);
  }

  ProteinIdentification protein;
  protein.setIdentifier("run0");
  vector<ProteinIdentification> proteins{protein};
  PeptideIdentificationList peptides;
  for (const auto& [sequence, rt] : peptides_rt)
  {
    PeptideIdentification peptide;
    peptide.setIdentifier("run0");
    peptide.setRT(rt);
    peptide.setMZ(AASequence::fromString(sequence).getMZ(2));
    peptide.setScoreType("q-value"); // identification data needs scored matches
    peptide.setHigherScoreBetter(false);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setCharge(2);
    hit.setScore(0.01);
    peptide.setHits({hit});
    peptides.push_back(peptide);
  }

  TestableFeatureFinderIdentificationAlgorithm algo;
  Param param = algo.getParameters();
  param.setValue("model:type", "none"); // the synthetic peaks need no elution model
  algo.setParameters(param);
  algo.setMSData(std::move(ms1));
  FeatureMap features;
  algo.run(peptides, proteins, features, FeatureMap(), "");
  TEST_EQUAL(algo.getPeptideCount(), 3)

  ostringstream captured_info;
  OPENMS_LOG_INFO.insert(captured_info);
  algo.writeStatistics(features);
  OPENMS_LOG_INFO.remove(captured_info);
  TEST_TRUE(captured_info.str().find("3 peptides identified") != string::npos)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


