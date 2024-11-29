// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <nlohmann/json.hpp>


using namespace std;

namespace OpenMS
{
  // A structure for storing a pointer to a ribo in the database, as well as the possible alternatives if it is ambiguous (eg a methyl group that for which we can't determine the localization)
  struct ParsedEntry_
    {
      unique_ptr<Ribonucleotide> ribo;
      String alternative_1;
      String alternative_2;
      bool isAmbiguous () { return !alternative_1.empty(); }
    };

  RibonucleotideDB::RibonucleotideDB() : max_code_length_(0)
  {
    // Modomics mods were retreived from https://www.genesilico.pl/modomics/api/modifications
    readFromJSON_("CHEMISTRY/Modomics.json");
    OPENMS_LOG_DEBUG << "Loading modomics RNA Modifications from "<< File::find("CHEMISTRY/Modomics.json") <<"\n";
    
    // We still use the old tsv format for custom mods
    readFromFile_("CHEMISTRY/Custom_RNA_modifications.tsv");
    OPENMS_LOG_DEBUG << "Loading custom RNA Modifications from "<< File::find("CHEMISTRY/Custom_RNA_modifications.tsv") <<"\n";
    
    if (File::exists("CHEMISTRY/User_Modifications.tsv"))
    {
      OPENMS_LOG_INFO << "Loading user specified Modifications from TSV\n";
    }
    if (File::exists("CHEMISTRY/User_Modifications.json"))
    {
      OPENMS_LOG_INFO << "Loading user specified Modifications from JSON\n";
    }
  }

  RibonucleotideDB* RibonucleotideDB::getInstance()
  {
    static RibonucleotideDB* db_ = new RibonucleotideDB(); // Meyers' singleton -> thread safe
    return db_;
  }

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
    String code = entry.at("short_name");
    // If we have an explicitly defined baseloss_formula
    if (auto e = entry.find("baseloss_formula"); e != entry.cend() && !e->is_null())
    {
      return EmpiricalFormula(*e);
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
    unique_ptr<Ribonucleotide> ribo (new Ribonucleotide());
    ribo->setName(entry.at("name"));
    String code = entry.at("short_name");
    ribo->setCode(code);
    // NewCode doesn't exist any more, we use the same shortname for compatibility
    ribo->setNewCode(code);

    // Handle moiety
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
                                              OPENMS_PRETTY_FUNCTION, msg, entry["reference_moiety"]);
    }
    
    if (entry.find("abbrev") != entry.cend())
    {
      ribo->setHTMLCode(entry.at("abbrev")); //This is the single letter unicode representation that only SOME mods have
    }
    ribo->setFormula(EmpiricalFormula(entry.at("formula")));
    if ( !(entry.find("mass_avg") == entry.cend()) && !(entry.at("mass_avg").is_null()))
    {
      ribo->setAvgMass(entry.at("mass_avg"));
    }
    if (std::abs(ribo->getAvgMass() - ribo->getFormula().getAverageWeight()) >= 0.01)
    {
      OPENMS_LOG_DEBUG << "Average mass of " << code << " differs substantially from its formula mass.\n";
    }

    if (auto e = entry.find("mass_monoiso"); e != entry.cend() && !e->is_null())
    {
      ribo->setMonoMass(*e);
    }
    else
    {
      OPENMS_LOG_DEBUG << "Monoisotopic mass of " << code << " is not defined. Calculating from formula\n";
      ribo->setMonoMass(ribo->getFormula().getMonoWeight());
    }
    if ( std::abs(ribo->getMonoMass() - ribo->getFormula().getMonoWeight()) >= 0.01)
    {
      OPENMS_LOG_DEBUG << "Average mass of " << code << " differs substantially from its formula mass.\n";
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
      parsed.alternative_1 = string(entry.at("alternatives").at(0)), parsed.alternative_2 = string(entry.at("alternatives").at(1)); // we always have exactly two ambiguities
    }

    parsed.ribo = std::move(ribo);
    return parsed;
  }

   // Read from a JSON file into a RibonucleotideDB
  void RibonucleotideDB::readFromJSON_(const std::string& path)
  {
    using json = nlohmann::json;

    String full_path = File::find(path);

    // the input file is Unicode encoded, so we need Qt to read it:
    QFile file(full_path.toQString());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    QTextStream source(&file);
    source.setCodec("UTF-8");
    Size line_count = 0;
    json mod_obj;
    try
    {
      mod_obj = json::parse(String(source.readAll()));
    }
    catch (Exception::ParseError& e)
    {
      OPENMS_LOG_ERROR << "Error: Failed to parse Modomics JSON. Reason:\n" << e.getName() << " - " << e.what() << endl;
      throw;
    }
    for (auto& element : mod_obj)
    {
      line_count++;
      try
      {
        // Throw an exception if we are straight up missing necessary elements of the JSON
        entryIsWellFormed_(element);

        ParsedEntry_ entry = parseEntry_(element);

        unique_ptr<Ribonucleotide> ribo = std::move(entry.ribo);
        if (entry.isAmbiguous()) // Handle the ambiguity map
        {
          ambiguity_map_[ribo->getCode()] = make_pair(getRibonucleotide(entry.alternative_1), getRibonucleotide(entry.alternative_2));
        }
        // there are some weird exotic mods in modomics that don't have codes. We ignore them
        if (ribo->getCode() != "")
        {
          code_map_[ribo->getCode()] = ribonucleotides_.size();
          max_code_length_ = max(max_code_length_, ribo->getCode().size());
          ribonucleotides_.push_back(std::move(ribo));
        }
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input element " << line_count << ". Reason:\n" << e.getName() << " - " << e.what() << "\nSkipping this line." << endl;
      }
    }
  }
  
  // Read entries from a TSV file
  void RibonucleotideDB::readFromFile_(const std::string& path)
  {
    String full_path = File::find(path);

    String header = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";

    // the input file is Unicode encoded, so we need Qt to read it:
    QFile file(full_path.toQString());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    QTextStream source(&file);
    source.setCodec("UTF-8");
    Size line_count = 1;
    String line = source.readLine();
    while (line[0] == '#') // skip leading comments
    {
      line = source.readLine();
      ++line_count;
    }
    if (!line.hasPrefix(header)) // additional columns are allowed
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, msg);
    }

    QChar prime(0x2032); // Unicode "prime" character
    while (!source.atEnd())
    {
      line_count++;
      QString row = source.readLine();

      // replace all "prime" characters with apostrophes (e.g. in "5'", "3'"):
      row.replace(prime, '\'');
      try
      {
        unique_ptr<Ribonucleotide> ribo = parseRow_(row.toStdString(), line_count);
        code_map_[ribo->getCode()] = ribonucleotides_.size();
        max_code_length_ = max(max_code_length_, ribo->getCode().size());
        ribonucleotides_.push_back(std::move(ribo));
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input line " << line_count << ". Reason:\n" << e.getName() << " - " << e.what() << "\nSkipping this line." << endl;
      }
    }
  }

  //Parse a row in a TSV file
  const unique_ptr<Ribonucleotide> RibonucleotideDB::parseRow_(const std::string& row, Size line_count)
  {
    vector<String> parts;
    String(row).split('\t', parts);
    if (parts.size() < 9)
    {
      String msg = "9 tab-separated fields expected, found " + String(parts.size()) + " in line " + String(line_count);
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, msg);
    }
    unique_ptr<Ribonucleotide> ribo (new Ribonucleotide());
    ribo->setName(parts[0]);
    if (parts[1].hasSuffix("QtRNA")) // use just "Q" instead of "QtRNA"
    {
      ribo->setCode(parts[1].chop(4));
    }
    else
    {
      ribo->setCode(parts[1]);
    }
    ribo->setNewCode(parts[2]);
    if (parts[3] == "preQ0base") // queuosine and its derivatives
    {
      ribo->setOrigin('G'); // queuosine replaces "G" in tRNA-Asp/Asn
    }
    else if (parts[3].size() == 1) // A, C, G, U
    {
      ribo->setOrigin(parts[3][0]);
    }
    // "parts[4]" is the Unicode equivalent to "parts[5]", so we can skip it
    ribo->setHTMLCode(parts[5]);
    if (!parts[6].empty() && (parts[6] != "-"))
    {
      ribo->setFormula(EmpiricalFormula(parts[6]));
    }
    if (!parts[7].empty() && (parts[7] != "None"))
    {
      ribo->setMonoMass(parts[7].toDouble());
      if ((ribo->getMonoMass() == 0.0) && (!ribo->getFormula().isEmpty()))
      {
        ribo->setMonoMass(ribo->getFormula().getMonoWeight());
      }
    }
    if (!parts[8].empty() && (parts[8] != "None"))
    {
      ribo->setAvgMass(parts[8].toDouble());
      if ((ribo->getAvgMass() == 0.0) && (!ribo->getFormula().isEmpty()))
      {
        ribo->setAvgMass(ribo->getFormula().getAverageWeight());
      }
    }
    // Modomics' "new code" contains information on terminal specificity:
    if ((!parts[2].empty()) && parts[2].back() == 'N') // terminal mod., exception: "GN"
    {
      if (parts[2].hasSubstring("55") || (parts[2] == "N"))
      {
        ribo->setTermSpecificity(Ribonucleotide::FIVE_PRIME);
      }
      else if (parts[2].hasSubstring("33"))
      {
        ribo->setTermSpecificity(Ribonucleotide::THREE_PRIME);
      }
    }
    else // default specificity is "ANYWHERE"; now set formula after base loss:
    {
      if (parts[1].front() == 'd') // handle deoxyribose, possibly with methyl mod
      {
        ribo->setBaselossFormula(EmpiricalFormula("C5H10O4"));
      }
      else if (parts[1].back() == 'm') // mod. attached to the ribose, not base
      {
        ribo->setBaselossFormula(EmpiricalFormula("C6H12O5"));
      }
      else if (parts[1].substr(parts[1].size() - 2) == "m*") // check if we have both a sulfer and a 2'-O methyl
      {
        ribo->setBaselossFormula(EmpiricalFormula("C6H12O5"));
      }
      else if (parts[1].back() == '?') // ambiguity code -> fill the map
      {
        if (parts.size() < 10)
        {
          String msg = "10th field expected for ambiguous modification in line " + String(line_count);
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, msg);
        }
        String code1 = parts[9].prefix(' '), code2 = parts[9].suffix(' ');
        ambiguity_map_[parts[1]] = make_pair(getRibonucleotide(code1), getRibonucleotide(code2));
      }
      else if ((parts[1] == "Ar(p)") || (parts[1] == "Gr(p)"))
      {
        ribo->setBaselossFormula(EmpiricalFormula("C10H19O21P"));
      }
    }
    return ribo;
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
