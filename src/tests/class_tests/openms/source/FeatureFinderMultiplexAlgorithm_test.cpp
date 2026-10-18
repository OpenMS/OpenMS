// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <cmath>
#include <map>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderMultiplexAlgorithm, "$Id$")

FeatureFinderMultiplexAlgorithm* ptr = 0;
FeatureFinderMultiplexAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderMultiplexAlgorithm())
{
  ptr = new FeatureFinderMultiplexAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~FeatureFinderMultiplexAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual void run()))
{
  MzMLFile mzml_file;
  MSExperiment exp;
  ConsensusMap result;
  
  mzml_file.getOptions().addMSLevel(1);
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_input.mzML"), exp);
  exp.updateRanges();
  
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_parameters.ini"), param);
  param = param.copy("FeatureFinderMultiplex:1:",true);
  param.remove("in");
  param.remove("out");
  param.remove("out_multiplets");
  param.remove("log");
  param.remove("debug");
  param.remove("threads");
  param.remove("no_progress");
  param.remove("force");
  param.remove("test");
  
  FeatureFinderMultiplexAlgorithm algorithm;
  algorithm.setParameters(param);
  algorithm.run(exp, true);
  result = algorithm.getConsensusMap();
  
  TEST_EQUAL(result.size(), 2);
  
  double L = result[0].getFeatures().begin()->getIntensity();
  double H = (++(result[0].getFeatures().begin()))->getIntensity();

  // Check that the HEAVY:LIGHT ratio is close to the expected 3:1 ratio
  TOLERANCE_ABSOLUTE(0.2);
  TEST_REAL_SIMILAR(H/L, 3.0);
  TOLERANCE_ABSOLUTE(1e-5);

  // every handle refers to the feature it was built from (its elution region), by unique id
  const FeatureMap& features = algorithm.getFeatureMap();
  TEST_EQUAL(features.size(), 4);
  std::map<UInt64, const Feature*> feature_of;
  for (const Feature& f : features) { feature_of[f.getUniqueId()] = &f; }
  for (const ConsensusFeature& cf : result)
  {
    for (const FeatureHandle& handle : cf.getFeatures())
    {
      TEST_EQUAL(feature_of.count(handle.getUniqueId()), 1)
      if (feature_of.count(handle.getUniqueId()) == 1)
      {
        const Feature& f = *feature_of[handle.getUniqueId()];
        TEST_REAL_SIMILAR(f.getMZ(), handle.getMZ())
        TEST_REAL_SIMILAR(f.getRT(), handle.getRT())
        TEST_EQUAL(f.getConvexHulls().empty(), false)
      }
    }
  }

  // the peaks of the detected patterns and the peaks the detection worked on are available
  TEST_EQUAL(algorithm.getBlacklist().empty(), false)
  TEST_EQUAL(algorithm.getCentroidedExperiment().empty(), false)
}
END_SECTION

// A triple SILAC run (centroided) with one peptide whose light form is absent: medium and heavy only.
// With knock-outs the two have to be found as one doublet, i.e. the doublet pattern must not
// require a peak at the absent light position.
START_SECTION((knock-out detection of a multiplet without its lightest channel))
{
  const double mz_medium = 600.0, charge = 2.0;
  const double lys4 = 4.0251069836, lys8 = 8.0141988132;
  const double mz_heavy = mz_medium + (lys8 - lys4) / charge;
  const std::vector<double> isotopes = {1.0, 0.6, 0.25, 0.08};

  MSExperiment exp;
  for (int i = 0; i < 60; ++i)
  {
    const double rt = 100.0 + i * 1.0;
    const double elution = std::exp(-0.5 * std::pow((rt - 130.0) / 5.0, 2));
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(rt);
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (const double mz : {mz_medium, mz_heavy})
    {
      for (Size k = 0; k < isotopes.size(); ++k)
      {
        if (1e6 * elution * isotopes[k] > 100.0)
        {
          spectrum.push_back(Peak1D(mz + k * Constants::C13C12_MASSDIFF_U / charge, 1e6 * elution * isotopes[k]));
        }
      }
    }
    spectrum.sortByPosition();
    exp.addSpectrum(spectrum);
  }
  exp.updateRanges();

  FeatureFinderMultiplexAlgorithm algorithm;
  Param param = algorithm.getParameters();
  param.setValue("algorithm:labels", "[][Lys4,Arg6][Lys8,Arg10]");
  param.setValue("algorithm:charge", "2:2");
  param.setValue("algorithm:isotopes_per_peptide", "3:4");
  param.setValue("algorithm:knock_out", "true");
  algorithm.setParameters(param);
  MSExperiment input = exp;
  algorithm.run(input, false);

  const ConsensusMap& multiplets = algorithm.getConsensusMap();
  TEST_EQUAL(multiplets.size(), 1)
  if (multiplets.size() == 1)
  {
    TEST_EQUAL(multiplets[0].size(), 2)
    TOLERANCE_ABSOLUTE(0.01);
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().begin()->getMZ(), mz_medium)
    TEST_REAL_SIMILAR(multiplets[0].getFeatures().rbegin()->getMZ(), mz_heavy)
  }

  // centroided input has a blacklist, too (the resolver tells absent from not quantifiable channels with it)
  Size blacklisted = 0;
  for (const MSSpectrum& spectrum : algorithm.getBlacklist()) { blacklisted += spectrum.size(); }
  TEST_NOT_EQUAL(blacklisted, 0)
}
END_SECTION

END_TEST
