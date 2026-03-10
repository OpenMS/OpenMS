// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>

using namespace std;

/// @brief Specialize nlohmann::adl_serializer for OpenMS::EmpiricalFormula so that nlohmann::json
// knows how to (de)serialize it: from_json constructs an EmpiricalFormula from a JSON string,to_json
// converts it back to its string form, all explicitly.
namespace nlohmann {
//
    template <>
    struct adl_serializer<OpenMS::EmpiricalFormula>
    {
        static void from_json(const json& j, OpenMS::EmpiricalFormula& ef)
        {
            std::string formula_string = j.get<std::string>();
            ef = OpenMS::EmpiricalFormula(formula_string);
        }

        static void to_json(json& j, const OpenMS::EmpiricalFormula& ef)
        {
            j = ef.toString();
        }
    };
}

namespace OpenMS
{

  namespace // anonymous namespace for file-local helpers
  {

    // A structure for storing a pointer to a ribo in the database, as well as the possible alternatives if it is ambiguous (eg a methyl group that for which we can't determine the localization)
    struct ParsedEntry_
    {
      unique_ptr<Ribonucleotide> ribo;
      String alternative_1;
      String alternative_2;
      bool isAmbiguous () { return !alternative_1.empty(); }
    };

    // All valid JSON ribonucleotides must at minimum have elements defining name, short_name, reference_moiety, and formula
    // @throw Exception::MissingInformation if some of the required info for the entry is missing
    void entryIsWellFormed_(const nlohmann::json::value_type& entry)
    {
      if (entry.find("name") == entry.cend())
      {
        String msg = "\"name\" entry missing for ribonucleotide";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg);
      }
      if (entry.find("short_name") == entry.cend())
      {
        String msg = "\"short_name\" entry missing for ribonucleotide";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg);
      }
      if (entry.find("reference_moiety") == entry.cend())
      {
        String msg = "\"reference_moiety\" entry missing for ribonucleotide";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg);
      }
      if (entry.find("formula") == entry.cend())
      {
        String msg = "\"formula\" entry missing for ribonucleotide";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    // Return the Empirical formula for the ribo with a base-loss. Ideally we store these in the JSON, otherwise its guessed from the code.
    EmpiricalFormula getBaseLossFormula_(const nlohmann::json::value_type& entry)
    {
      String code = entry.at("short_name").get<std::string>();
      // If we have an explicitly defined baseloss_formula
      if (auto e = entry.find("baseloss_formula"); e != entry.cend() && !e->is_null())
      {
        return EmpiricalFormula(e->get<std::string>());
      }
      //TODO: Calculate base loss formula from SMILES
      else // If we don't have a defined baseloss_formula calculate it from our shortCode
      {
        if (code.hasPrefix('d')) // handle deoxyribose, possibly with methyl mod
        {
          return EmpiricalFormula("C5H10O4");
        }
        else if (code.hasSuffix('m')) // mod. attached to the ribose, not base
        {
          return EmpiricalFormula("C6H12O5");
        }
        else if (code.hasSuffix("m*")) // check if we have both a sulfer and a 2'-O methyl
        {
          return EmpiricalFormula("C6H12O5");
        }
        else if (code.hasSuffix("Ar(p)") ||  code.hasSuffix("Gr(p)"))
        {
          return EmpiricalFormula("C10H19O21P");
        }
        else
        {
          return EmpiricalFormula("C5H10O5");
        }
      }
    }

    // Generate an entry from a JSON object.
    ParsedEntry_ parseEntry_(const nlohmann::json::value_type& entry)
    {
      ParsedEntry_ parsed;
      auto ribo = std::make_unique<Ribonucleotide>();
      ribo->setName(entry.at("name").template get<std::string>());
      String code = entry.at("short_name").get<std::string>();
      ribo->setCode(code);
      // NewCode doesn't exist any more, we use the same shortname for compatibility
      ribo->setNewCode(code);

      // Handle moiety
      if (!entry["reference_moiety"].is_array())
      {
        String msg = "\"reference_moiety\" must be a JSON array";
        throw Exception::InvalidValue(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg, entry["reference_moiety"].dump());
      }
      if (entry["reference_moiety"].size() == 1 && string(entry.at("reference_moiety").at(0)).length() == 1)
      {
        ribo->setOrigin(string(entry.at("reference_moiety").at(0))[0]);
        ribo->setTermSpecificity(Ribonucleotide::ANYWHERE); // due to format changes we get the terminal specificity from the moieties, modomics contains base specific terminals, but they can be represented by the wild-card ones
      }
      else if (entry["reference_moiety"].size() == 4) // if all moieties are possible it might be a terminal
      {
        ribo->setOrigin('X'); // Use X as any unmodified
        if (code.hasSuffix("pN"))
        {
          ribo->setTermSpecificity(Ribonucleotide::FIVE_PRIME);
        }
        else if (code.hasSuffix("p") && code.hasPrefix("N"))
        {
          ribo->setTermSpecificity(Ribonucleotide::THREE_PRIME);
        }
        else
        {
          ribo->setTermSpecificity(Ribonucleotide::ANYWHERE); //other nonspecific mods
        }
      }
      else
      {
        String msg = "we don't support bases with multiple reference moieties or multicharacter moieties.";
        throw Exception::InvalidValue(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, msg, entry["reference_moiety"].dump());
      }

      if (entry.find("abbrev") != entry.cend())
      {
        ribo->setHTMLCode(entry.at("abbrev").get<std::string>()); //This is the single letter unicode representation that only SOME mods have
      }
      ribo->setFormula(entry.at("formula").get<OpenMS::EmpiricalFormula>());

      if (auto e = entry.find("mass_avg"); e != entry.cend() && !e->is_null())
      {
        ribo->setAvgMass(e->get<double>());
      }
      if (std::abs(ribo->getAvgMass() - ribo->getFormula().getAverageWeight()) >= 0.01)
      {
        OPENMS_LOG_DEBUG << "Average mass of " << code << " differs substantially from its formula mass.\n";
      }

      if (auto e = entry.find("mass_monoiso"); e != entry.cend() && !e->is_null())
      {
        ribo->setMonoMass(e->get<double>());
      }
      else
      {
        OPENMS_LOG_DEBUG << "Monoisotopic mass of " << code << " is not defined. Calculating from formula\n";
        ribo->setMonoMass(ribo->getFormula().getMonoWeight());
      }
      if ( std::abs(ribo->getMonoMass() - ribo->getFormula().getMonoWeight()) >= 0.01)
      {
        OPENMS_LOG_DEBUG << "Monoisotopic mass of " << code << " differs substantially from its formula mass.\n";
      }

      // Handle base loss formula
      ribo->setBaselossFormula(getBaseLossFormula_(entry));

      // Handle ambiguities
      if (code.hasSuffix('?') || code.hasSuffix("?*")) // ambiguity code -> fill the map
      {
        if (!entry.contains("alternatives"))
        {
          String msg = "Ambiguous mod without alternative found in " + code;
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code, msg);
        }
        const auto& alts = entry.at("alternatives");
        if (!alts.is_array() || alts.size() < 2)
        {
          String msg = "\"alternatives\" must be an array with at least 2 entries for " + code;
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code, msg);
        }
        parsed.alternative_1 = alts.at(0).get<std::string>();
        parsed.alternative_2 = alts.at(1).get<std::string>();
      }

      parsed.ribo = std::move(ribo);
      return parsed;
    }

  } // anonymous namespace

  ModomicsJSONDataProvider::ModomicsJSONDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<RibonucleotideEntry> ModomicsJSONDataProvider::loadRibonucleotides()
  {
    using json = nlohmann::json;

    String full_path = File::find(filename_);

    // Use std::filesystem::path to support non-ASCII paths on Windows (wide-string open)
    std::ifstream file{std::filesystem::path{std::string(full_path)}};
    if (!file.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    Size line_count = 0;
    json mod_obj;
    try
    {
      std::ostringstream content;
      content << file.rdbuf();
      mod_obj = json::parse(content.str());
    }
    catch (nlohmann::json::parse_error& e)
    {
      OPENMS_LOG_ERROR << "Error: Failed to parse Modomics JSON. Reason:\n" << e.what() << endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   full_path, String("JSON parse error: ") + e.what());
    }

    std::vector<RibonucleotideEntry> result;

    for (auto& element : mod_obj)
    {
      line_count++;
      try
      {
        // Throw an exception if we are straight up missing necessary elements of the JSON
        entryIsWellFormed_(element);

        ParsedEntry_ entry = parseEntry_(element);

        // there are some weird exotic mods in modomics that don't have codes. We ignore them
        if (entry.ribo->getCode() != "")
        {
          RibonucleotideEntry result_entry;
          result_entry.ribo = std::move(entry.ribo);
          if (entry.isAmbiguous())
          {
            result_entry.alternative_code_1 = entry.alternative_1;
            result_entry.alternative_code_2 = entry.alternative_2;
          }
          result.push_back(std::move(result_entry));
        }
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input element " << line_count << ". Reason:\n" << e.getName() << " - " << e.what() << "\nSkipping this element." << endl;
      }
      catch (std::exception& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input element " << line_count << ". Reason:\n" << e.what() << "\nSkipping this element." << endl;
      }
    }

    return result;
  }

} // namespace OpenMS
