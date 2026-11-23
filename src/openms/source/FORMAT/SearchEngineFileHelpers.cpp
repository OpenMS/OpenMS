// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include "SearchEngineFileHelpers.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS::Internal::SearchEngineFileHelpers
{

  void handlePTMs(const std::string& modification_line,
                  const std::string& modifications_filename,
                  const bool monoisotopic,
                  std::map<std::string, std::vector<std::string> >& PTMname_residues_mass_type,
                  const char* function)
  {
    PTMname_residues_mass_type.clear();
    // to store the information about modifications from the ptm xml file
    map<std::string, pair<std::string, std::string> > ptm_informations;
    if (!modification_line.empty())       // if modifications are used look whether whether composition and residues (and type and name) is given, the name (type) is used (then the modifications file is needed) or only the mass and residues (and type and name) is given
    {
      vector<std::string> modifications, mod_parts;
      StringUtils::split(modification_line, ':', modifications);       // get the single modifications

      // to get masses from a formula
      EmpiricalFormula add_formula, substract_formula;

      std::string types = "OPT#FIX#";
      std::string name, residues, mass, type;

      // 0 - mass; 1 - composition; 2 - ptm name
      Int mass_or_composition_or_name(-1);

      for (const std::string& mod_i : modifications)
      {
        if (mod_i.empty())
        {
          continue;
        }
        // clear the formulae
        add_formula = substract_formula = EmpiricalFormula();
        name = residues = mass = type = "";

        // get the single parts of the modification string
        StringUtils::split(mod_i, ',', mod_parts);
        mass_or_composition_or_name = -1;

        // check whether the first part is a mass, composition or name

        // check whether it is a mass
        try
        {
          mass = mod_parts.front();
          // to check whether the first part is a mass, it is converted into a float and then back into a string and compared to the given string
          // remove + signs because they don't appear in a float
          if (StringUtils::hasPrefix(mass, "+"))
          {
            mass.erase(0, 1);
          }
          if (StringUtils::hasSuffix(mass, "+"))
          {
            mass.erase(mass.length() - 1, 1);
          }
          if (StringUtils::hasSuffix(mass, "-"))             // a - sign at the end will not be converted
          {
            mass.erase(mass.length() - 1, 1);
            mass.insert(0, "-");
          }
          // if it is a mass
          if (!StringUtils::toStr(StringUtils::toFloat(mass)).empty()) // just check if conversion does not throw, i.e. consumes the whole string
          {
            mass_or_composition_or_name = 0;
          }
        }
        catch (Exception::ConversionError & /*c_e*/)
        {
          mass_or_composition_or_name = -1;
        }

        // check whether it is a name (look it up in the corresponding file)
        if (mass_or_composition_or_name == -1)
        {
          if (ptm_informations.empty())             // if the ptm xml file has not been read yet, read it
          {
            if (!File::exists(modifications_filename))
            {
              throw Exception::FileNotFound(__FILE__, __LINE__, function, modifications_filename);
            }
            if (!File::readable(modifications_filename))
            {
              throw Exception::FileNotReadable(__FILE__, __LINE__, function, modifications_filename);
            }

            // getting all available modifications from a file
            PTMXMLFile().load(modifications_filename, ptm_informations);
          }
          // if the modification cannot be found
          if (ptm_informations.contains(mod_parts.front()))
          {
            mass = ptm_informations[mod_parts.front()].first;             // composition
            residues = ptm_informations[mod_parts.front()].second;             // residues
            name = mod_parts.front();             // name

            mass_or_composition_or_name = 2;
          }
        }

        // check whether it's an empirical formula / if a composition was given, get the mass
        if (mass_or_composition_or_name == -1)
        {
          mass = mod_parts.front();
        }
        if (mass_or_composition_or_name == -1 || mass_or_composition_or_name == 2)
        {
          // check whether there is a positive and a negative formula
          std::string::size_type pos = mass.find("-");
          try
          {
            if (pos != std::string::npos)
            {
              add_formula = EmpiricalFormula(StringUtils::substr(mass, 0, pos));
              substract_formula = EmpiricalFormula(StringUtils::substr(mass, ++pos));
            }
            else
            {
              add_formula = EmpiricalFormula(mass);
            }
            // sum up the masses
            if (monoisotopic)
            {
              mass =StringUtils::toStr(add_formula.getMonoWeight() - substract_formula.getMonoWeight());
            }
            else
            {
              mass =StringUtils::toStr(add_formula.getAverageWeight() - substract_formula.getAverageWeight());
            }
            if (mass_or_composition_or_name == -1)
            {
              mass_or_composition_or_name = 1;
            }
          }
          catch (Exception::ParseError & /*pe*/)
          {
            PTMname_residues_mass_type.clear();
            throw Exception::ParseError(__FILE__, __LINE__, function, mod_i, "There's something wrong with this modification. Aborting!");
          }
        }

        // now get the residues
        mod_parts.erase(mod_parts.begin());
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            PTMname_residues_mass_type.clear();
            throw Exception::ParseError(__FILE__, __LINE__, function, mod_i, "No residues for modification given. Aborting!");
          }

          // get the residues
          residues = mod_parts.front();
          StringUtils::substitute(residues, '*', 'X');
          StringUtils::toUpper(residues);
          mod_parts.erase(mod_parts.begin());
        }

        // get the type
        if (mod_parts.empty())
        {
          type = "OPT";
        }
        else
        {
          type = mod_parts.front();
          StringUtils::toUpper(type);
          if (types.contains(type))
          {
            mod_parts.erase(mod_parts.begin());
          }
          else
          {
            type = "OPT";
          }
        }

        if (mod_parts.size() > 1)
        {
          PTMname_residues_mass_type.clear();
          throw Exception::ParseError(__FILE__, __LINE__, function, mod_i, "There's something wrong with the type of this modification. Aborting!");
        }

        // get the name
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            name = "PTM_" + StringUtils::toStr(PTMname_residues_mass_type.size());
          }
          else
          {
            name = mod_parts.front();
          }
        }

        // insert the modification
        if (!PTMname_residues_mass_type.contains(name))
        {
          PTMname_residues_mass_type[name] = vector<std::string>(3);
          PTMname_residues_mass_type[name][0] = residues;
          // mass must not have more than 5 digits after the . (otherwise the test may fail)
          PTMname_residues_mass_type[name][1] = StringUtils::substr(mass, 0, mass.find(".") + 6);
          PTMname_residues_mass_type[name][2] = type;
        }
        else
        {
          PTMname_residues_mass_type.clear();
          throw Exception::ParseError(__FILE__, __LINE__, function, mod_i, "There's already a modification with this name. Aborting!");
        }
      }
    }
  }

  void getACAndACType(std::string line, std::string& accession, std::string& accession_type)
  {
    std::string swissprot_prefixes = "JLOPQUX";
    /// @todo replace this by general FastA implementation? (Martin)
    accession.clear();
    accession_type.clear();
    // if it's a FASTA line
    if (StringUtils::hasPrefix(line, ">"))
    {
      line.erase(0, 1);
    }
    if (!line.empty() && (line[line.length() - 1] < 33))
    {
      line.resize(line.length() - 1);
    }
    StringUtils::trim(line);

    // if it's a swissprot accession
    if (StringUtils::hasPrefix(line, "tr") || StringUtils::hasPrefix(line, "sp"))
    {
      accession = StringUtils::substr(line, 3, line.find('|', 3) - 3);
      accession_type = "SwissProt";
    }
    else if (StringUtils::hasPrefix(line, "gi"))
    {
      std::string::size_type snd(line.find('|', 3));
      std::string::size_type third(0);
      if (snd != std::string::npos)
      {
        third = line.find('|', ++snd) + 1;

        accession = StringUtils::substr(line, third, line.find('|', third) - third);
        accession_type = StringUtils::substr(line, snd, third - 1 - snd);
      }
      if (accession_type == "gb")
      {
        accession_type = "GenBank";
      }
      else if (accession_type == "emb")
      {
        accession_type = "EMBL";
      }
      else if (accession_type == "dbj")
      {
        accession_type = "DDBJ";
      }
      else if (accession_type == "ref")
      {
        accession_type = "NCBI";
      }
      else if ((accession_type == "sp") || (accession_type == "tr"))
      {
        accession_type = "SwissProt";
      }
      else if (accession_type == "gnl")
      {
        accession_type = accession;
        snd = line.find('|', third);
        third = line.find('|', ++snd);
        if (third != std::string::npos)
        {
          accession = StringUtils::substr(line, snd, third - snd);
        }
        else
        {
          third = line.find(' ', snd);
          if (third != std::string::npos)
          {
            accession = StringUtils::substr(line, snd, third - snd);
          }
          else
          {
            accession = StringUtils::substr(line, snd);
          }
        }
      }
      else
      {
        std::string::size_type pos1(line.find('(', 0));
        std::string::size_type pos2(0);
        if (pos1 != std::string::npos)
        {
          pos2 = line.find(')', ++pos1);
          if (pos2 != std::string::npos)
          {
            accession = StringUtils::substr(line, pos1, pos2 - pos1);
            if ((accession.size() == 6) && (std::string(swissprot_prefixes).find(accession[0], 0) != std::string::npos))
            {
              accession_type = "SwissProt";
            }
            else
            {
              accession.clear();
            }
          }
        }
        if (accession.empty())
        {
          accession_type = "gi";
          if (snd != std::string::npos)
          {
            accession = StringUtils::substr(line, 3, snd - 4);
          }
          else
          {
            snd = line.find(' ', 3);
            if (snd != std::string::npos)
            {
              accession = StringUtils::substr(line, 3, snd - 3);
            }
            else
            {
              accession = StringUtils::substr(line, 3);
            }
          }
        }
      }
    }
    else if (StringUtils::hasPrefix(line, "ref"))
    {
      accession = StringUtils::substr(line, 4, line.find('|', 4) - 4);
      accession_type = "NCBI";
    }
    else if (StringUtils::hasPrefix(line, "gnl"))
    {
      line.erase(0, 3);
      accession_type = StringUtils::substr(line, 0, line.find('|', 0));
      accession = StringUtils::substr(line, accession_type.length() + 1);
    }
    else if (StringUtils::hasPrefix(line, "lcl"))
    {
      line.erase(0, 4);
      accession_type = "lcl";
      accession = line;
    }
    else
    {
      std::string::size_type pos1(line.find('(', 0));
      std::string::size_type pos2(0);
      if (pos1 != std::string::npos)
      {
        pos2 = line.find(')', ++pos1);
        if (pos2 != std::string::npos)
        {
          accession = StringUtils::substr(line, pos1, pos2 - pos1);
          if ((accession.size() == 6) && (std::string(swissprot_prefixes).find(accession[0], 0) != std::string::npos))
          {
            accession_type = "SwissProt";
          }
          else
          {
            accession.clear();
          }
        }
      }
      if (accession.empty())
      {
        pos1 = line.find('|');
        accession = StringUtils::substr(line, 0, pos1);
        if ((accession.size() == 6) && (std::string(swissprot_prefixes).find(accession[0], 0) != std::string::npos))
        {
          accession_type = "SwissProt";
        }
        else
        {
          pos1 = line.find(' ');
          accession = StringUtils::substr(line, 0, pos1);
          if ((accession.size() == 6) && (std::string(swissprot_prefixes).find(accession[0], 0) != std::string::npos))
          {
            accession_type = "SwissProt";
          }
          else
          {
            accession = StringUtils::substr(line, 0, 6);
            if (std::string(swissprot_prefixes).find(accession[0], 0) != std::string::npos)
            {
              accession_type = "SwissProt";
            }
            else
            {
              accession.clear();
            }
          }
        }
      }
    }
    if (accession.empty())
    {
      accession = StringUtils::trim(line);
      accession_type = "unknown";
    }
  }

} // namespace OpenMS::Internal::SearchEngineFileHelpers
