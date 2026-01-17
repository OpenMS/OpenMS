// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/Scores.h>

namespace OpenMS
{

  // Static member initialization
  std::map<Scores::IDType, std::set<String>> Scores::id_type_to_str_;
  std::map<Scores::IDType, bool> Scores::id_type_to_better_;
  bool Scores::maps_initialized_ = false;

  void Scores::initializeMaps_()
  {
    if (maps_initialized_) return;

    id_type_to_str_ = {
      //TODO introduce real meaningful score names for XTandem, Mascot etc. (e.g., hyperscore)
      {IDType::RAW, {"svm", "MS:1001492", "XTandem", "OMSSA", "SEQUEST:xcorr", "Mascot", "mvh", "hyperscore", "ln(hyperscore)"}},
      //TODO find out reasonable raw scores for SES that provide E-Values as main score or see below
      //TODO there is no test for spectraST idXML, so I don't know its score
      //TODO check if we should combine RAW and RAW_EVAL:
      // What if a SE does not have an e-value score (spectrast, OMSSA, crux/sequest, myrimatch),
      // then you need additional if's/try's
      {IDType::RAW_EVAL, {"expect", "SpecEValue", "E-Value", "evalue", "MS:1002053", "MS:1002257"}},
      {IDType::PP, {"Posterior Probability"}},
      {IDType::PEP, {"Posterior Error Probability", "pep", "PEP", "posterior_error_probability", "MS:1001493"}}, // TODO add CV terms
      {IDType::FDR, {"FDR", "fdr", "false discovery rate"}},
      {IDType::QVAL, {"q-value", "qvalue", "MS:1001491", "q-Value", "qval"}}
    };

    id_type_to_better_ = {
      {IDType::RAW, true}, //TODO this might actually not always be true
      {IDType::RAW_EVAL, false},
      {IDType::PP, true},
      {IDType::PEP, false},
      {IDType::FDR, false},
      {IDType::QVAL, false}
    };

    maps_initialized_ = true;
  }

  bool Scores::isScoreType(const String& score_name, IDType type)
  {
    initializeMaps_();

    String chopped = score_name;
    if (chopped.hasSuffix("_score"))
    {
      chopped = chopped.chop(6);
    }
    const std::set<String>& possible_types = id_type_to_str_.at(type);
    return possible_types.find(chopped) != possible_types.end();
  }

  Scores::IDType Scores::parseIDType(const String& score_type_in)
  {
    String score_type = score_type_in;
    if (score_type.hasSuffix("_score"))
    {
      score_type = score_type.chop(6);
    }
    score_type.toLower();
    score_type.erase(std::remove_if(score_type.begin(), score_type.end(),
              [](unsigned char c) { return c == '-' || c == '_' || c == ' '; }),
              score_type.end());

    const std::map<String, IDType> s_to_type = {
      {"raw", IDType::RAW},
      {"rawevalue", IDType::RAW_EVAL},
      {"qvalue", IDType::QVAL},
      {"fdr", IDType::FDR},
      {"falsediscoveryrate", IDType::FDR},
      {"pep", IDType::PEP},
      {"posteriorerrorprobability", IDType::PEP},
      {"posteriorprobability", IDType::PP},
      {"pp", IDType::PP}
    };

    if (auto it = s_to_type.find(score_type); it != s_to_type.end())
    {
      return it->second;
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, String("Unknown score type '") + score_type_in + "'.");
    }
  }

  bool Scores::isHigherBetter(IDType type)
  {
    initializeMaps_();
    return id_type_to_better_.at(type);
  }

  std::vector<String> Scores::getAllIDScoreNames()
  {
    initializeMaps_();
    std::vector<String> names;
    for (const auto& [type, name_set] : id_type_to_str_)
    {
      for (const auto& name : name_set)
      {
        names.push_back(name);
      }
    }
    return names;
  }

  const std::set<String>& Scores::getIDNamesForType(IDType type)
  {
    initializeMaps_();
    return id_type_to_str_.at(type);
  }

  bool Scores::findIDTypeByName(const String& name, IDType& type)
  {
    initializeMaps_();
    for (const auto& [scoretype, names] : id_type_to_str_)
    {
      if (names.find(name) != names.end())
      {
        type = scoretype;
        return true;
      }
    }
    return false;
  }

} // namespace OpenMS
