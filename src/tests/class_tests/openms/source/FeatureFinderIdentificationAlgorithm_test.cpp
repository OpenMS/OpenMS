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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <cmath>
#include <set>
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

  using PeptidesRT = std::vector<std::pair<std::string, double>>;

  PeakMap makeSyntheticMS1(const PeptidesRT& apexes)
  {
    PeakMap ms1;
    for (int scan = 0; scan < 111; ++scan)
    {
      MSSpectrum spectrum;
      spectrum.setMSLevel(1);
      spectrum.setRT(80.0 + scan);
      spectrum.setNativeID("scan=" + std::to_string(scan + 1));
      for (const auto& [sequence, apex] : apexes)
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
    return ms1;
  }

  PeptideIdentificationList makePeptideIDs(const PeptidesRT& peptides_rt)
  {
    PeptideIdentificationList peptides;
    for (const auto& [sequence, rt] : peptides_rt)
    {
      PeptideIdentification peptide;
      peptide.setIdentifier("run0");
      peptide.setRT(rt);
      peptide.setMZ(AASequence::fromString(sequence).getMZ(2));
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(sequence));
      hit.setCharge(2);
      peptide.setHits({hit});
      peptides.push_back(peptide);
    }
    return peptides;
  }
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
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  PeptideIdentification peptide;
  peptide.setHits({hit});

  Feature feature;
  feature.setIntensity(1.0);
  feature.setPeptideIdentifications({peptide});
  FeatureMap features;
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
}
END_SECTION

START_SECTION([EXTRA] run() counts the distinct identified peptides for the summary)
{
  // The summary statistics report the peptides of the run, not a stale or uninitialised count: three distinct
  // peptides, one of them identified twice, eluting in a small MS1-only experiment.
  PeakMap ms1 = makeSyntheticMS1({{"PEPTIDEK", 100.0}, {"SAMPLER", 130.0}, {"ELVISLIVESK", 160.0}});

  ProteinIdentification protein;
  protein.setIdentifier("run0");
  vector<ProteinIdentification> proteins{protein};
  PeptideIdentificationList peptides = makePeptideIDs({{"PEPTIDEK", 100.0}, {"PEPTIDEK", 101.0}, {"SAMPLER", 130.0}, {"ELVISLIVESK", 160.0}});

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

START_SECTION([EXTRA] getChromatograms() returns the chromatograms of every extraction batch)
{
  const PeptidesRT peptides_rt = {{"PEPTIDEK", 100.0}, {"SAMPLER", 130.0}, {"ELVISLIVESK", 160.0}};
  PeakMap ms1 = makeSyntheticMS1(peptides_rt);

  ProteinIdentification protein;
  protein.setIdentifier("run0");
  vector<ProteinIdentification> proteins{protein};
  PeptideIdentificationList peptides = makePeptideIDs(peptides_rt);

  FeatureFinderIdentificationAlgorithm algo;
  Param param = algo.getParameters();
  param.setValue("model:type", "none");
  param.setValue("extract:batch_size", 1);
  Param param_keep = param;
  param_keep.setValue("keep_chromatograms", "true");
  algo.setParameters(param_keep);
  algo.setMSData(PeakMap(ms1));
  FeatureMap features;
  algo.run(peptides, proteins, features, FeatureMap(), "");

  std::set<std::string> transition_ids, chrom_ids;
  for (const auto& transition : algo.getLibrary().getTransitions())
  {
    transition_ids.insert(transition.getNativeID());
  }
  for (const auto& chrom : algo.getChromatograms().getChromatograms())
  {
    chrom_ids.insert(chrom.getNativeID());
  }
  TEST_NOT_EQUAL(transition_ids.size(), 0)
  TEST_EQUAL(algo.getChromatograms().getNrChromatograms(), transition_ids.size())
  TEST_TRUE(chrom_ids == transition_ids)
  TEST_NOT_EQUAL(features.size(), 0)

  FeatureFinderIdentificationAlgorithm algo_default;
  algo_default.setParameters(param);
  algo_default.setMSData(std::move(ms1));
  FeatureMap features_default;
  algo_default.run(peptides, proteins, features_default, FeatureMap(), "");
  TEST_EQUAL(algo_default.getChromatograms().getNrChromatograms(), 0)
  TEST_EQUAL(features_default.size(), features.size())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


