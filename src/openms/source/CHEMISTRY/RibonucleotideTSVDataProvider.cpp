// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <filesystem>
#include <fstream>

using namespace std;

namespace OpenMS
{

  namespace // anonymous namespace for file-local helpers
  {

    // Parse a single row from a TSV file into a RibonucleotideEntry.
    // Extracted from RibonucleotideDB::parseRow_().
    RibonucleotideEntry parseRow_(const std::string& row, Size line_count)
    {
      vector<String> parts;
      String(row).split('\t', parts);
      if (parts.size() < 9)
      {
        String msg = "9 tab-separated fields expected, found " + String(parts.size()) + " in line " + String(line_count);
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, msg);
      }
      RibonucleotideEntry entry;
      auto ribo = std::make_unique<Ribonucleotide>();
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
        else if (parts[1].size() >= 2 && parts[1].substr(parts[1].size() - 2) == "m*") // check if we have both a sulfer and a 2'-O methyl
        {
          ribo->setBaselossFormula(EmpiricalFormula("C6H12O5"));
        }
        else if (parts[1].back() == '?') // ambiguity code -> store alternative codes
        {
          if (parts.size() < 10)
          {
            String msg = "10th field expected for ambiguous modification in line " + String(line_count);
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, msg);
          }
          String code1 = parts[9].prefix(' '), code2 = parts[9].suffix(' ');
          entry.alternative_code_1 = code1;
          entry.alternative_code_2 = code2;
        }
        else if ((parts[1] == "Ar(p)") || (parts[1] == "Gr(p)"))
        {
          ribo->setBaselossFormula(EmpiricalFormula("C10H19O21P"));
        }
      }
      entry.ribo = std::move(ribo);
      return entry;
    }

  } // anonymous namespace

  RibonucleotideTSVDataProvider::RibonucleotideTSVDataProvider(const String& filename)
    : filename_(filename)
  {
  }

  std::vector<RibonucleotideEntry> RibonucleotideTSVDataProvider::loadRibonucleotides()
  {
    String full_path = File::find(filename_);

    String header = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";

    // Use std::filesystem::path to support non-ASCII paths on Windows (wide-string open)
    std::ifstream file{std::filesystem::path{std::string(full_path)}};
    if (!file.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, full_path);
    }

    Size line_count = 1;
    String line;
    std::getline(file, line);
    while (!line.empty() && line[0] == '#') // skip leading comments
    {
      std::getline(file, line);
      ++line_count;
    }
    if (!line.hasPrefix(header)) // additional columns are allowed
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, msg);
    }

    const std::string prime = "\xE2\x80\xB2"; // UTF-8 encoding of Unicode "prime" character U+2032

    std::vector<RibonucleotideEntry> result;

    while (std::getline(file, line))
    {
      line_count++;

      // replace all "prime" characters with apostrophes (e.g. in "5'", "3'"):
      std::string::size_type pos = 0;
      while ((pos = line.find(prime, pos)) != std::string::npos)
      {
        line.replace(pos, prime.size(), "'");
        ++pos;
      }
      try
      {
        RibonucleotideEntry entry = parseRow_(line, line_count);
        result.push_back(std::move(entry));
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to parse input line " << line_count << ". Reason:\n" << e.getName() << " - " << e.what() << "\nSkipping this line." << endl;
      }
    }

    return result;
  }

} // namespace OpenMS
