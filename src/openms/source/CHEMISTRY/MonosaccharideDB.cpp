// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/MonosaccharideDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

namespace OpenMS
{
  MonosaccharideDB::MonosaccharideDB()
  {
    loadFromJSON_();
  }

  MonosaccharideDB* MonosaccharideDB::getInstance()
  {
    // Meyers' singleton - thread safe in C++11 and later
    static MonosaccharideDB* instance_ = new MonosaccharideDB();
    return instance_;
  }

  void MonosaccharideDB::loadFromJSON_()
  {
    using json = nlohmann::json;

    String path = "CHEMISTRY/monosaccharides.json";
    String full_path = File::find(path);

    OPENMS_LOG_DEBUG << "Loading monosaccharide data from " << full_path << "\n";

    // Use std::filesystem::path to support non-ASCII paths on Windows (wide-string open)
    std::ifstream file{std::filesystem::path{std::string(full_path)}};
    if (!file.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    json data;
    try
    {
      std::ostringstream content;
      content << file.rdbuf();
      data = json::parse(content.str());
    }
    catch (const json::parse_error& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        full_path, String("JSON parse error: ") + e.what());
    }

    // Check for required structure
    if (!data.contains("monosaccharides"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        full_path, "JSON file missing 'monosaccharides' key");
    }

    const auto& monos = data["monosaccharides"];

    for (auto it = monos.begin(); it != monos.end(); ++it)
    {
      const String symbol = it.key();
      const auto& entry = it.value();

      Monosaccharide mono;
      mono.symbol = symbol;

      // Required fields
      if (!entry.contains("mass"))
      {
        OPENMS_LOG_WARN << "Monosaccharide '" << symbol << "' missing 'mass' field, skipping\n";
        continue;
      }
      mono.mass = entry["mass"].get<double>();

      // Optional fields with defaults
      if (entry.contains("name"))
      {
        mono.name = entry["name"].get<std::string>();
      }
      else
      {
        mono.name = symbol;
      }

      if (entry.contains("formula"))
      {
        mono.formula = entry["formula"].get<std::string>();
      }

      if (entry.contains("synonyms") && entry["synonyms"].is_array())
      {
        for (const auto& syn : entry["synonyms"])
        {
          String synonym = syn.get<std::string>();
          mono.synonyms.push_back(synonym);
          // Map synonym to primary symbol for lookup
          synonym_to_symbol_[synonym] = symbol;
        }
      }

      monosaccharides_[symbol] = mono;
      // Also map primary symbol to itself for uniform lookup
      synonym_to_symbol_[symbol] = symbol;
    }

    OPENMS_LOG_DEBUG << "Loaded " << monosaccharides_.size() << " monosaccharides\n";
  }

  bool MonosaccharideDB::hasSymbol(const String& symbol) const
  {
    return synonym_to_symbol_.find(symbol) != synonym_to_symbol_.end();
  }

  const MonosaccharideDB::Monosaccharide* MonosaccharideDB::getMonosaccharide(const String& symbol) const
  {
    // First check if it's a synonym
    auto syn_it = synonym_to_symbol_.find(symbol);
    if (syn_it == synonym_to_symbol_.end())
    {
      return nullptr;
    }

    // Look up by primary symbol
    auto mono_it = monosaccharides_.find(syn_it->second);
    if (mono_it == monosaccharides_.end())
    {
      return nullptr;
    }

    return &(mono_it->second);
  }

  const MonosaccharideDB::Monosaccharide& MonosaccharideDB::getMonosaccharideOrThrow(const String& symbol) const
  {
    const Monosaccharide* mono = getMonosaccharide(symbol);
    if (mono == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, symbol);
    }
    return *mono;
  }

  std::vector<String> MonosaccharideDB::getAllSymbols() const
  {
    std::vector<String> symbols;
    symbols.reserve(monosaccharides_.size());
    for (const auto& pair : monosaccharides_)
    {
      symbols.push_back(pair.first);
    }
    return symbols;
  }

  Size MonosaccharideDB::getNumberOfMonosaccharides() const
  {
    return monosaccharides_.size();
  }

} // namespace OpenMS
