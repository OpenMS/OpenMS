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

  const Scores::Maps_& Scores::getMaps_()
  {
    static const Maps_ maps = [] {
      Maps_ m;
      m.type_to_str = {
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

      m.type_to_better = {
        {IDType::RAW, true}, //TODO this might actually not always be true
        {IDType::RAW_EVAL, false},
        {IDType::PP, true},
        {IDType::PEP, false},
        {IDType::FDR, false},
        {IDType::QVAL, false}
      };

      return m;
    }();
    return maps;
  }

  bool Scores::isScoreType(const std::string& score_name, IDType type)
  {
    const auto& maps = getMaps_();

    std::string chopped = score_name;
    if (StringUtils::hasSuffix(chopped, "_score"))
    {
      chopped = StringUtils::chop(chopped, 6);
    }
    const std::set<std::string>& possible_types = maps.type_to_str.at(type);
    return possible_types.contains(chopped);
  }

  Scores::IDType Scores::parseIDType(const std::string& score_type_in)
  {
    std::string score_type = score_type_in;
    if (StringUtils::hasSuffix(score_type, "_score"))
    {
      score_type = StringUtils::chop(score_type, 6);
    }
    StringUtils::toLower(score_type);
    score_type.erase(std::remove_if(score_type.begin(), score_type.end(),
              [](unsigned char c) { return c == '-' || c == '_' || c == ' '; }),
              score_type.end());

    const std::map<std::string, IDType> s_to_type = {
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
                                          OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Unknown score type '") + score_type_in + "'.");
    }
  }

  bool Scores::isHigherBetter(IDType type)
  {
    const auto& maps = getMaps_();
    return maps.type_to_better.at(type);
  }

  std::vector<std::string> Scores::getAllIDScoreNames()
  {
    const auto& maps = getMaps_();
    std::vector<std::string> names;
    for (const auto& [type, name_set] : maps.type_to_str)
    {
      for (const auto& name : name_set)
      {
        names.push_back(name);
      }
    }
    return names;
  }

  const std::set<std::string>& Scores::getIDNamesForType(IDType type)
  {
    const auto& maps = getMaps_();
    return maps.type_to_str.at(type);
  }

  bool Scores::findIDTypeByName(const std::string& name, IDType& type)
  {
    const auto& maps = getMaps_();
    for (const auto& [scoretype, names] : maps.type_to_str)
    {
      if (names.contains(name))
      {
        type = scoretype;
        return true;
      }
    }
    return false;
  }

  std::string Scores::normalizeScoreName(const std::string& score_name)
  {
    std::string result = score_name;
    if (StringUtils::hasSuffix(result, "_score"))
    {
      result = StringUtils::chop(result, 6);
    }
    return result;
  }

  bool Scores::isKnownScoreType(const std::string& score_name)
  {
    const auto& maps = getMaps_();
    std::string normalized = normalizeScoreName(score_name);

    for (const auto& [type, names] : maps.type_to_str)
    {
      if (names.contains(normalized))
      {
        return true;
      }
    }
    return false;
  }

} // namespace OpenMS
