// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  PepNovoOutfile::PepNovoOutfile() = default;

  PepNovoOutfile::PepNovoOutfile(const PepNovoOutfile &) = default;

  PepNovoOutfile::~PepNovoOutfile() = default;

  PepNovoOutfile & PepNovoOutfile::operator=(const PepNovoOutfile &) = default;

  bool PepNovoOutfile::operator==(const PepNovoOutfile &) const
  {
    return true;
  }

  void
  PepNovoOutfile::load(
    const std::string & result_filename,
    PeptideIdentificationList & peptide_identifications,
    ProteinIdentification & protein_identification,
    const double & score_threshold,
    const IndexPosMappingType & index_to_precursor,
    const map<std::string, std::string> & pnovo_modkey_to_mod_id
    )
  {
    // generally used variables
    StringList substrings;
    map<std::string, Int> columns;
    PeptideHit peptide_hit;

    std::string line, identifier, filename, sequence, sequence_with_mods;
    std::string score_type = "PepNovo", version = "unknown";

    DateTime datetime = DateTime::now();     // there's no date given from PepNovo
    protein_identification.setDateTime(datetime);

    peptide_identifications.clear();
    PeptideIdentification peptide_identification;
    protein_identification = ProteinIdentification();

    // open the result
    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      if (!File::exists(result_filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else if (!File::readable(result_filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
    }

    Size line_number(0);     // used to report in which line an error occurred
    Size id_count(0);        // number of IDs seen (not necessarily the ones finally returned)

    getSearchEngineAndVersion(result_filename, protein_identification);
    //if information could not be retrieved from the outfile use defaults
    if (protein_identification.getSearchEngineVersion().empty())
    {
      protein_identification.setSearchEngine("PepNovo");
      protein_identification.setSearchEngineVersion(version);
    }
    identifier = protein_identification.getSearchEngine() + "_" + datetime.getDate();
    protein_identification.setIdentifier(identifier);

    map<std::string, std::string> mod_mask_map;
    const vector<std::string> & mods = protein_identification.getSearchParameters().variable_modifications;
    for (const std::string& mod_it : mods)
    {
      if (mod_it.empty())
        continue;
      //cout<<*mod_it<<endl;
      if (pnovo_modkey_to_mod_id.contains(mod_it))
      {
        //cout<<keys_to_id.find(*mod_it)->second<<endl;
        const ResidueModification* tmp_mod = ModificationsDB::getInstance()->getModification(pnovo_modkey_to_mod_id.find(mod_it)->second);
        if (StringUtils::prefix(mod_it, 1) == "^" || StringUtils::prefix(mod_it, 1) == "$")
        {
          mod_mask_map[mod_it] = "(" + tmp_mod->getId() + ")";
        }
        else
        {
          mod_mask_map[mod_it] =StringUtils::toStr(tmp_mod->getOrigin()) + "(" + tmp_mod->getId() + ")";
        }
      }
      else
      {
        if (StringUtils::prefix(mod_it, 1) != "^" && StringUtils::prefix(mod_it, 1) != "$")
        {
          mod_mask_map[mod_it] = StringUtils::prefix(mod_it, 1) + "[" + StringUtils::substr(mod_it, 1) + "]";
          //cout<<mod_mask_map[*mod_it]<<endl;
        }
        else
        {
          mod_mask_map[mod_it] = "[" + mod_it + "]";
          //cout<<mod_mask_map[*mod_it]<<endl;
        }
      }
    }


    Size index;
    while (getline(result_file, line))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1); // remove weird EOL character
      }
      StringUtils::trim(line);
      ++line_number;
      if (StringUtils::hasPrefix(line, ">> "))         // >> 1 /home/shared/pepnovo/4611_raw_ms2_picked.mzXML.1001.2.dta
      {
        ++id_count;
        if (!peptide_identification.empty() && !peptide_identification.getHits().empty())
        {
          peptide_identifications.push_back(peptide_identification);
        }

        StringUtils::split(line, ' ', substrings);
        //std::string index = File::basename(StringUtils::substr(line, line.find(' ', strlen(">> ")) + 1));
        if (substrings.size() < 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns (spectrum Id) in file in line " + StringUtils::toStr(line_number) + std::string(" (should be 2 or more)!"), result_filename);
        }

        try
        {
          index = StringUtils::toInt32(StringUtils::trimmed(substrings[2]));
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected an index number in line " + StringUtils::toStr(line_number) + " at position 2 (line was: '" + line + "')!", result_filename);
        }

        //cout<<"INDEX: "<<index<<endl;
        peptide_identification = PeptideIdentification();
        bool success = false;
        if (!index_to_precursor.empty())
        {
          if (index_to_precursor.contains(index))
          {
            peptide_identification.setRT(index_to_precursor.find(index)->second.first);
            peptide_identification.setMZ(index_to_precursor.find(index)->second.second);
            success = true;
          }
          else
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Index '" + StringUtils::toStr(index) + "' in line '" + line + "' not found in index table (line was: '" + line + "')!", result_filename);
          }
        }

        if (!success)
        { // try to reconstruct from title entry (usually sensible when MGF is supplied to PepNovo)
          try
          {
            if (substrings.size() >= 4)
            {
              StringList parts = ListUtils::create<std::string>(substrings[3], '_');
              if (parts.size() >= 2)
              {
                peptide_identification.setRT(StringUtils::toDouble(parts[1]));
                peptide_identification.setMZ(StringUtils::toDouble(parts[0]));
                success = true;
              }
            }
          }
          catch (...)
          {

          }
          if (!success)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor could not be reconstructed from title '" + substrings[3] + "' in line '" + line + "' (line was: '" + line + "')!", result_filename);
          }
        }
        peptide_identification.setSignificanceThreshold(score_threshold);
        peptide_identification.setScoreType(score_type);
        peptide_identification.setIdentifier(identifier);
      }
      else if (StringUtils::hasPrefix(line, "#Index"))         // #Index  Prob    Score   N-mass  C-Mass  [M+H]   Charge  Sequence
      {
        if (columns.empty())           // map the column names to their column number
        {
          StringUtils::split(line, '\t', substrings);
          for (vector<std::string>::const_iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i)
          {
            if ((*s_i) == "#Index")
            {
              columns["Index"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "RnkScr")
            {
              columns["RnkScr"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "PnvScr")
            {
              columns["PnvScr"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "N-Gap")
            {
              columns["N-Gap"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "C-Gap")
            {
              columns["C-Gap"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "[M+H]")
            {
              columns["[M+H]"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "Charge")
            {
              columns["Charge"] = s_i - substrings.begin();
            }
            else if ((*s_i) == "Sequence")
            {
              columns["Sequence"] = s_i - substrings.begin();
            }
          }

          if (columns.size() != 8)
          {
            result_file.close();
            result_file.clear();
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns in file in line " + StringUtils::toStr(line_number) + std::string(" (should be 8)!"), result_filename);
          }
        }
        while (getline(result_file, line))
        {
          ++line_number;
          if (!line.empty() && (line[line.length() - 1] < 33))
          {
            line.resize(line.length() - 1);
          }
          StringUtils::trim(line);

          if (line.empty())
          {
            break;
          }
          StringUtils::split(line, '\t', substrings);
          if (!substrings.empty())
          {
            if (substrings.size() != 8)
            {
              result_file.close();
              result_file.clear();
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns in file in line " + StringUtils::toStr(line_number) + std::string(" (should be 8)!"), result_filename);
            }
            if (StringUtils::toFloat(substrings[columns["RnkScr"]]) >= score_threshold)
            {
              peptide_hit = PeptideHit();
              peptide_hit.setCharge(StringUtils::toInt32(substrings[columns["Charge"]]));
              peptide_hit.setRank(StringUtils::toInt32(substrings[columns["Index"]]) + 1);
              peptide_hit.setScore(StringUtils::toFloat(substrings[columns["RnkScr"]]));
              peptide_hit.setMetaValue("PnvScr", StringUtils::toFloat(substrings[columns["PnvScr"]]));
              peptide_hit.setMetaValue("N-Gap", StringUtils::toFloat(substrings[columns["N-Gap"]]));
              peptide_hit.setMetaValue("C-Gap", StringUtils::toFloat(substrings[columns["C-Gap"]]));
              peptide_hit.setMetaValue("MZ", StringUtils::toFloat(substrings[columns["[M+H]"]]));
              sequence = substrings[columns["Sequence"]];


              for (map<std::string, std::string>::iterator mask_it = mod_mask_map.begin(); mask_it != mod_mask_map.end(); ++mask_it)
              {
                if (StringUtils::hasPrefix(mask_it->first, "^") && StringUtils::hasSubstring(sequence, mask_it->first))
                {
                  StringUtils::substitute(sequence, mask_it->first, "");
                  sequence = mask_it->second + sequence;
                }
                //cout<<mask_it->first<<" "<<mask_it->second<<endl;
                StringUtils::substitute(sequence, mask_it->first, mask_it->second);
              }
              peptide_hit.setSequence(AASequence::fromString(sequence));
              peptide_identification.insertHit(peptide_hit);
            }
          }
        }
      }
    }
    if (!peptide_identifications.empty() || !peptide_identification.getHits().empty())
    {
      peptide_identifications.push_back(peptide_identification);
    }

    result_file.close();
    result_file.clear();

    OPENMS_LOG_INFO << "Parsed " << id_count << " ids, retained " << peptide_identifications.size() << "." << std::endl;

  }

  void
  PepNovoOutfile::getSearchEngineAndVersion(
    const std::string & pepnovo_output_without_parameters_filename,
    ProteinIdentification & protein_identification)
  {
    ifstream pepnovo_output_without_parameters(pepnovo_output_without_parameters_filename.c_str());
    if (!pepnovo_output_without_parameters)
    {
      if (!File::exists(pepnovo_output_without_parameters_filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pepnovo_output_without_parameters_filename);
      }
      else if (!File::readable(pepnovo_output_without_parameters_filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pepnovo_output_without_parameters_filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pepnovo_output_without_parameters_filename);
      }
    }

    ProteinIdentification::SearchParameters search_param;
    // searching for something like this: PepNovo v1.03
    std::string line;
    vector<std::string> substrings;
    while (getline(pepnovo_output_without_parameters, line))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      StringUtils::trim(line);
      if (line.empty())
      {
        continue;
      }
      if (StringUtils::hasPrefix(line, "PepNovo"))
      {
        StringUtils::split(line, ', ', substrings);
        if (substrings.size() == 2) //previous version of PepNovo
        {
          protein_identification.setSearchEngine(StringUtils::trim(substrings[0]));
          protein_identification.setSearchEngineVersion(StringUtils::trim(substrings[1])); //else something is strange and we use defaults later
        }
        else
        {
          StringUtils::split(line, ' ', substrings);
          if (substrings.size() == 3)
          {
            protein_identification.setSearchEngine(StringUtils::trim(substrings[0]));
            protein_identification.setSearchEngineVersion(StringUtils::trim(substrings[2])); //else something is strange and we use defaults later
          }
        }
      }
      if (StringUtils::hasPrefix(line, "PM"))
      {
        StringUtils::split(line, ' ', substrings);
        search_param.precursor_mass_tolerance = StringUtils::toFloat(substrings.back());
      }
      if (StringUtils::hasPrefix(line, "Fragment"))
      {
        StringUtils::split(line, ' ', substrings);
        search_param.fragment_mass_tolerance = StringUtils::toFloat(substrings.back());
      }
      if (StringUtils::hasPrefix(line, "PTM"))
      {
        StringUtils::split(line, ':', substrings);
        substrings.erase(substrings.begin());
        for (vector<std::string>::iterator ptm_it = substrings.begin(); ptm_it != substrings.end(); ++ptm_it)
        {
          StringUtils::trim(*ptm_it);
        }
        if (!substrings.empty() && substrings[0] != "None")
        {
          search_param.variable_modifications = substrings;
        }
      }
      if (StringUtils::hasPrefix(line, ">>"))
      {
        break;
      }
    }
    protein_identification.setSearchParameters(search_param);
  }

} //namespace OpenMS
