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
    vector<PeptideIdentification> & peptide_identifications,
    ProteinIdentification & protein_identification,
    const double & score_threshold,
    const IndexPosMappingType & index_to_precursor,
    const map<String, String> & pnovo_modkey_to_mod_id
    )
  {
    // generally used variables
    StringList substrings;
    map<String, Int> columns;
    PeptideHit peptide_hit;

    String
      line,
      score_type = "PepNovo",
      version = "unknown",
      identifier,
      filename,
      sequence,
      sequence_with_mods;

    DateTime datetime = DateTime::now();     // there's no date given from PepNovo
    protein_identification.setDateTime(datetime);

    peptide_identifications.clear();
    PeptideIdentification peptide_identification;
    protein_identification = ProteinIdentification();

    // open the result
    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
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

    map<String, String> mod_mask_map;
    const vector<String> & mods = protein_identification.getSearchParameters().variable_modifications;
    for (const String& mod_it : mods)
    {
      if (mod_it.empty())
        continue;
      //cout<<*mod_it<<endl;
      if (pnovo_modkey_to_mod_id.find(mod_it) != pnovo_modkey_to_mod_id.end())
      {
        //cout<<keys_to_id.find(*mod_it)->second<<endl;
        const ResidueModification* tmp_mod = ModificationsDB::getInstance()->getModification(pnovo_modkey_to_mod_id.find(mod_it)->second);
        if (mod_it.prefix(1) == "^" || mod_it.prefix(1) == "$")
        {
          mod_mask_map[mod_it] = "(" + tmp_mod->getId() + ")";
        }
        else
        {
          mod_mask_map[mod_it] = String(tmp_mod->getOrigin()) + "(" + tmp_mod->getId() + ")";
        }
      }
      else
      {
        if (mod_it.prefix(1) != "^" && mod_it.prefix(1) != "$")
        {
          mod_mask_map[mod_it] = mod_it.prefix(1) + "[" + mod_it.substr(1) + "]";
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
      line.trim();
      ++line_number;
      if (line.hasPrefix(">> "))         // >> 1 /home/shared/pepnovo/4611_raw_ms2_picked.mzXML.1001.2.dta
      {
        ++id_count;
        if (!peptide_identification.empty() && !peptide_identification.getHits().empty())
        {
          peptide_identifications.push_back(peptide_identification);
        }

        line.split(' ', substrings);
        //String index = File::basename(line.substr(line.find(' ', strlen(">> ")) + 1));
        if (substrings.size() < 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns (spectrum Id) in file in line " + String(line_number) + String(" (should be 2 or more)!"), result_filename);
        }

        try
        {
          index = substrings[2].trim().toInt();
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected an index number in line " + String(line_number) + String(" at position 2 (line was: '" + line + "')!"), result_filename);
        }

        //cout<<"INDEX: "<<index<<endl;
        peptide_identification = PeptideIdentification();
        bool success = false;
        if (!index_to_precursor.empty())
        {
          if (index_to_precursor.find(index) != index_to_precursor.end())
          {
            peptide_identification.setRT(index_to_precursor.find(index)->second.first);
            peptide_identification.setMZ(index_to_precursor.find(index)->second.second);
            success = true;
          }
          else
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Index '" + String(index) + String("' in line '" + line + "' not found in index table (line was: '" + line + "')!"), result_filename);
          }
        }

        if (!success)
        { // try to reconstruct from title entry (usually sensible when MGF is supplied to PepNovo)
          try
          {
            if (substrings.size() >= 4)
            {
              StringList parts = ListUtils::create<String>(substrings[3], '_');
              if (parts.size() >= 2)
              {
                peptide_identification.setRT(parts[1].toDouble());
                peptide_identification.setMZ(parts[0].toDouble());
                success = true;
              }
            }
          }
          catch (...)
          {

          }
          if (!success)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor could not be reconstructed from title '" + substrings[3] + String("' in line '" + line + "' (line was: '" + line + "')!"), result_filename);
          }
        }
        peptide_identification.setSignificanceThreshold(score_threshold);
        peptide_identification.setScoreType(score_type);
        peptide_identification.setIdentifier(identifier);
      }
      else if (line.hasPrefix("#Index"))         // #Index  Prob    Score   N-mass  C-Mass  [M+H]   Charge  Sequence
      {
        if (columns.empty())           // map the column names to their column number
        {
          line.split('\t', substrings);
          for (vector<String>::const_iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i)
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
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns in file in line " + String(line_number) + String(" (should be 8)!"), result_filename);
          }
        }
        while (getline(result_file, line))
        {
          ++line_number;
          if (!line.empty() && (line[line.length() - 1] < 33))
          {
            line.resize(line.length() - 1);
          }
          line.trim();

          if (line.empty())
          {
            break;
          }
          line.split('\t', substrings);
          if (!substrings.empty())
          {
            if (substrings.size() != 8)
            {
              result_file.close();
              result_file.clear();
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough columns in file in line " + String(line_number) + String(" (should be 8)!"), result_filename);
            }
            if (substrings[columns["RnkScr"]].toFloat() >= score_threshold)
            {
              peptide_hit = PeptideHit();
              peptide_hit.setCharge(substrings[columns["Charge"]].toInt());
              peptide_hit.setRank(substrings[columns["Index"]].toInt() + 1);
              peptide_hit.setScore(substrings[columns["RnkScr"]].toFloat());
              peptide_hit.setMetaValue("PnvScr", substrings[columns["PnvScr"]].toFloat());
              peptide_hit.setMetaValue("N-Gap", substrings[columns["N-Gap"]].toFloat());
              peptide_hit.setMetaValue("C-Gap", substrings[columns["C-Gap"]].toFloat());
              peptide_hit.setMetaValue("MZ", substrings[columns["[M+H]"]].toFloat());
              sequence = substrings[columns["Sequence"]];


              for (map<String, String>::iterator mask_it = mod_mask_map.begin(); mask_it != mod_mask_map.end(); ++mask_it)
              {
                if (mask_it->first.hasPrefix("^") && sequence.hasSubstring(mask_it->first))
                {
                  sequence.substitute(mask_it->first, "");
                  sequence = mask_it->second + sequence;
                }
                //cout<<mask_it->first<<" "<<mask_it->second<<endl;
                sequence.substitute(mask_it->first, mask_it->second);
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
    const String & pepnovo_output_without_parameters_filename,
    ProteinIdentification & protein_identification)
  {
    ifstream pepnovo_output_without_parameters(pepnovo_output_without_parameters_filename.c_str());
    if (!pepnovo_output_without_parameters)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pepnovo_output_without_parameters_filename);
    }

    ProteinIdentification::SearchParameters search_param;
    // searching for something like this: PepNovo v1.03
    String line;
    vector<String> substrings;
    while (getline(pepnovo_output_without_parameters, line))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      line.trim();
      if (line.empty())
      {
        continue;
      }
      if (line.hasPrefix("PepNovo"))
      {
        line.split(',', substrings);
        if (substrings.size() == 2) //previous version of PepNovo
        {
          protein_identification.setSearchEngine(substrings[0].trim());
          protein_identification.setSearchEngineVersion(substrings[1].trim()); //else something is strange and we use defaults later
        }
        else
        {
          line.split(' ', substrings);
          if (substrings.size() == 3)
          {
            protein_identification.setSearchEngine(substrings[0].trim());
            protein_identification.setSearchEngineVersion(substrings[2].trim()); //else something is strange and we use defaults later
          }
        }
      }
      if (line.hasPrefix("PM"))
      {
        line.split(' ', substrings);
        search_param.precursor_mass_tolerance = substrings.back().toFloat();
      }
      if (line.hasPrefix("Fragment"))
      {
        line.split(' ', substrings);
        search_param.fragment_mass_tolerance = substrings.back().toFloat();
      }
      if (line.hasPrefix("PTM"))
      {
        line.split(':', substrings);
        substrings.erase(substrings.begin());
        for (vector<String>::iterator ptm_it = substrings.begin(); ptm_it != substrings.end(); ++ptm_it)
        {
          ptm_it->trim();
        }
        if (!substrings.empty() && substrings[0] != "None")
        {
          search_param.variable_modifications = substrings;
        }
      }
      if (line.hasPrefix(">>"))
      {
        break;
      }
    }
    protein_identification.setSearchParameters(search_param);
  }

} //namespace OpenMS
