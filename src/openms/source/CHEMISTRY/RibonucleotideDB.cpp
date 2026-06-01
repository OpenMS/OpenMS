// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;
namespace OpenMS
{
  RibonucleotideDB::RibonucleotideDB() : max_code_length_(0)
  {
    std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers;
    providers.push_back(std::make_unique<ModomicsJSONDataProvider>("CHEMISTRY/Modomics.json"));
    providers.push_back(std::make_unique<RibonucleotideTSVDataProvider>("CHEMISTRY/Custom_RNA_modifications.tsv"));
    loadFromProviders_(providers);
  }

  RibonucleotideDB::RibonucleotideDB(std::vector<std::unique_ptr<RibonucleotideDataProvider>> providers)
    : max_code_length_(0)
  {
    loadFromProviders_(providers);
  }

  RibonucleotideDB* RibonucleotideDB::getInstance()
  {
    static RibonucleotideDB* db_ = new RibonucleotideDB(); // Meyers' singleton -> thread safe
    return db_;
  }

  void RibonucleotideDB::loadFromProviders_(std::vector<std::unique_ptr<RibonucleotideDataProvider>>& providers)
  {
    // Collect deferred ambiguity entries (need all ribonucleotides loaded first)
    std::vector<std::tuple<std::string, String, String>> deferred_ambiguities;

    for (auto& provider : providers)
    {
      auto entries = provider->loadRibonucleotides();
      for (auto& entry : entries)
      {
        if (entry.isAmbiguous())
        {
          deferred_ambiguities.emplace_back(
            entry.ribo->getCode(),
            entry.alternative_code_1,
            entry.alternative_code_2);
        }
        code_map_[entry.ribo->getCode()] = ribonucleotides_.size();
        max_code_length_ = std::max(max_code_length_, entry.ribo->getCode().size());
        ribonucleotides_.push_back(std::move(entry.ribo));
      }
    }

    // Now resolve ambiguities (all codes are indexed)
    for (const auto& [code, alt1, alt2] : deferred_ambiguities)
    {
      try
      {
        ambiguity_map_[code] = std::make_pair(getRibonucleotide(alt1), getRibonucleotide(alt2));
      }
      catch (Exception::ElementNotFound& e)
      {
        OPENMS_LOG_ERROR << "Error resolving ambiguity for " << code << ": " << e.what() << std::endl;
      }
    }
  }

  RibonucleotideDB::ConstRibonucleotidePtr RibonucleotideDB::getRibonucleotide(const std::string& code)
  {
    std::unordered_map<std::string, Size>::const_iterator pos = code_map_.find(code);
    if (pos == code_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code);
    }
    return ribonucleotides_[pos->second].get();
  }


  RibonucleotideDB::ConstRibonucleotidePtr RibonucleotideDB::getRibonucleotidePrefix(const std::string& seq)
  {
    std::string prefix = seq.substr(0, max_code_length_);
    while (!prefix.empty())
    {
      auto pos = code_map_.find(prefix);
      if (pos != code_map_.end())
      {
        return ribonucleotides_[pos->second].get();
      }
      prefix = prefix.substr(0, prefix.size() - 1);
    }
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, seq);
  }


  pair<RibonucleotideDB::ConstRibonucleotidePtr, RibonucleotideDB::ConstRibonucleotidePtr> RibonucleotideDB::getRibonucleotideAlternatives(const std::string& code)
  {
    auto pos = ambiguity_map_.find(code);
    if (pos == ambiguity_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code);
    }
    return pos->second;
  }
} // namespace OpenMS
