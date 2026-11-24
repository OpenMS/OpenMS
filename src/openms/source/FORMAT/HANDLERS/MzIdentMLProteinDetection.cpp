// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLProteinDetection.h>

#include <OpenMS/ANALYSIS/ID/Scores.h>

namespace OpenMS::Internal
{
  namespace
  {
    // OpenMS inference engine names (ProteinIdentification::getInferenceEngine()) and their PSI-MS software terms
    struct EngineTerm
    {
      const char* engine;
      const char* term;
    };
    constexpr EngineTerm engine_terms[] = {
      {"Epifany", "EPIFANY"},                          // MS:1003118
      {"TOPPProteinInference", "TOPP ProteinInference"}, // MS:1002203
      {"Percolator", "Percolator"}                     // MS:1001490
    };
  }

  const std::string MzIdentMLProteinDetection::CLUSTER_SCORE_PARAM = "cluster score";
  const std::string MzIdentMLProteinDetection::SCORE_TYPE_PARAM = "protein score type";
  const std::string MzIdentMLProteinDetection::HIGHER_SCORE_BETTER_PARAM = "higher protein score better";

  std::string MzIdentMLProteinDetection::proteinScoreAccession(const std::string& score_type, const std::string& inference_engine)
  {
    if (Scores::isScoreType(score_type, Scores::IDType::PP))
    {
      // the engine-specific term says where the probability comes from; both are read back as "Posterior Probability"
      return inference_engine == "Epifany" ? "MS:1003119" /* EPIFANY:Protein posterior probability */ : "MS:1002367" /* probability for proteins */;
    }
    if (Scores::isScoreType(score_type, Scores::IDType::QVAL))
    {
      return "MS:1001869"; // protein-level q-value
    }
    return "";
  }

  std::string MzIdentMLProteinDetection::groupScoreAccession(const std::string& score_type)
  {
    if (Scores::isScoreType(score_type, Scores::IDType::PP))
    {
      return "MS:1002376"; // protein group-level probability
    }
    if (Scores::isScoreType(score_type, Scores::IDType::QVAL))
    {
      return "MS:1002373"; // protein group-level q-value
    }
    return "";
  }

  std::string MzIdentMLProteinDetection::groupThresholdAccession(const std::string& score_type)
  {
    if (Scores::isScoreType(score_type, Scores::IDType::QVAL))
    {
      return "MS:1002369"; // protein group-level global FDR: groups with a q-value up to the threshold were kept
    }
    return "MS:1002486"; // protein group-level statistical threshold
  }

  std::string MzIdentMLProteinDetection::scoreTypeOfAccession(const std::string& accession)
  {
    if (accession == "MS:1003119" || accession == "MS:1002367" || accession == "MS:1002376")
    {
      return "Posterior Probability";
    }
    if (accession == "MS:1001869" || accession == "MS:1003117" || accession == "MS:1002373")
    {
      return "q-value";
    }
    return "";
  }

  std::string MzIdentMLProteinDetection::softwareTermName(const std::string& inference_engine)
  {
    for (const EngineTerm& et : engine_terms)
    {
      if (inference_engine == et.engine) { return et.term; }
    }
    return "";
  }

  std::string MzIdentMLProteinDetection::inferenceEngineOfSoftware(const std::string& software_name)
  {
    for (const EngineTerm& et : engine_terms)
    {
      if (software_name == et.term) { return et.engine; }
    }
    return software_name;
  }
} // namespace OpenMS::Internal
