// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow, Hendrik Weisser, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{

    MascotXMLHandler::MascotXMLHandler(ProteinIdentification& protein_identification, vector<PeptideIdentification>& id_data, const String& filename, map<String, vector<AASequence> >& modified_peptides, const SpectrumMetaDataLookup& lookup):
      XMLHandler(filename, ""), protein_identification_(protein_identification),
      id_data_(id_data), peptide_identification_index_(0), actual_title_(""),
      modified_peptides_(modified_peptides), lookup_(lookup),
      no_rt_error_(false)
    {
    }

    MascotXMLHandler::~MascotXMLHandler()
    = default;

    void MascotXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
    {
      static const XMLCh* s_protein_accession = xercesc::XMLString::transcode("accession");
      static const XMLCh* s_queries_query_number = xercesc::XMLString::transcode("number");
      static const XMLCh* s_peptide_query = xercesc::XMLString::transcode("query");

      tag_ = String(sm_.convert(qname));
      // cerr << "open: " << tag_ << endl;

      tags_open_.push_back(tag_);

      if (tag_ == "mascot_search_results")
      {
        major_version_ = this->attributeAsString_(attributes, "majorVersion");
        minor_version_ = this->attributeAsString_(attributes, "minorVersion");
        no_rt_error_ = false; // reset for every new file
      }
      else if (tag_ == "protein")
      {
        String attribute_value = attributeAsString_(attributes, s_protein_accession);
        actual_protein_hit_.setAccession(attribute_value);
      }
      else if (tag_ == "query")
      {
        actual_query_ = attributeAsInt_(attributes, s_queries_query_number);
      }
      else if (tag_ == "peptide" || tag_ == "u_peptide" || tag_ == "q_peptide")
      {
        Int attribute_value = attributeAsInt_(attributes, s_peptide_query);
        peptide_identification_index_ = attribute_value - 1;

        if (peptide_identification_index_ > id_data_.size())
        {
          fatalError(LOAD, "No or conflicting header information present (make sure to use the 'show_header=1' option in the ./export_dat.pl script)");
        }
      }
    }

    void MascotXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      tag_ = String(sm_.convert(qname)).trim();
      // cerr << "close: " << tag_ << endl;

      if (tags_open_.empty())
      {
        fatalError(LOAD, String("Closing tag ") + tag_ + " not matched by opening tag", __LINE__);
      }

      tags_open_.pop_back();

      if (tag_ == "NumQueries")
      {
        id_data_.resize(character_buffer_.trim().toInt());
      }
      else if (tag_ == "prot_score")
      {
        actual_protein_hit_.setScore(character_buffer_.trim().toInt());
      }
      else if (tag_ == "pep_exp_mz")
      {
        id_data_[peptide_identification_index_].setMZ(
          character_buffer_.trim().toDouble());
      }
      else if (tag_ == "pep_scan_title")
      {
        // extract RT (and possibly m/z, if not already set) from title:
        String title = character_buffer_.trim();
        SpectrumMetaDataLookup::SpectrumMetaData meta;
        SpectrumMetaDataLookup::MetaDataFlags flags = SpectrumMetaDataLookup::MDF_RT;
        if (!id_data_[peptide_identification_index_].hasMZ())
        {
          flags |= SpectrumMetaDataLookup::MDF_PRECURSORMZ;
        }
        try
        {
          lookup_.getSpectrumMetaData(title, meta, flags);
          id_data_[peptide_identification_index_].setRT(meta.rt);
          // have we looked up the m/z value?
          if ((flags & SpectrumMetaDataLookup::MDF_PRECURSORMZ) == SpectrumMetaDataLookup::MDF_PRECURSORMZ)
          {
            id_data_[peptide_identification_index_].setMZ(meta.precursor_mz);
          }
        }
        catch (...)
        {
          String msg = "<pep_scan_title> element has unexpected format '" + title + "'. Could not extract spectrum meta data.";
          error(LOAD, msg);
        }
        // did it work?
        if (!id_data_[peptide_identification_index_].getRT())
        {
          if (!no_rt_error_) // report the error only the first time
          {
            String msg = "Could not extract RT value ";
            if (!lookup_.empty())
            {
              msg += "or a matching spectrum reference ";
            }
            msg += "from <pep_scan_title> element with format '" + title + "'. Try adjusting the 'scan_regex' parameter.";
            error(LOAD, msg);
          }
          no_rt_error_ = true;
        }
      }
      else if (tag_ == "pep_exp_z")
      {
        actual_peptide_hit_.setCharge(character_buffer_.trim().toInt());
      }
      else if (tag_ == "pep_score")
      {
        actual_peptide_hit_.setScore(character_buffer_.trim().toDouble());
      }
      else if (tag_ == "pep_expect")
      {
        // @todo what E-value flag? (andreas)
        actual_peptide_hit_.metaRegistry().registerName("EValue", "E-value of e.g. Mascot searches", "");
        actual_peptide_hit_.setMetaValue("EValue", character_buffer_.trim().toDouble());
      }
      else if (tag_ == "pep_homol")
      {
        id_data_[peptide_identification_index_].setSignificanceThreshold(character_buffer_.trim().toDouble());
      }
      else if (tag_ == "pep_ident")
      {
        double temp_homology = 0;
        double temp_identity = 0;

        // According to Matrix Science the homology threshold is only used if it
        // exists and is smaller than the identity threshold.
        temp_homology = id_data_[peptide_identification_index_].getSignificanceThreshold();
        temp_identity = character_buffer_.trim().toDouble();
        actual_peptide_hit_.setMetaValue("homology_threshold", temp_homology);
        actual_peptide_hit_.setMetaValue("identity_threshold", temp_identity);
        if (temp_homology > temp_identity || temp_homology == 0)
        {
          id_data_[peptide_identification_index_].setSignificanceThreshold(temp_identity);
        }
      }
      else if (tag_ == "pep_seq")
      {
        AASequence temp_aa_sequence = AASequence::fromString(character_buffer_.trim());
        
        // if everything is just read from the MascotXML file
        if (modified_peptides_.empty())
        {
          // fixed modifications
          for (vector<String>::const_iterator it = search_parameters_.fixed_modifications.begin(); it != search_parameters_.fixed_modifications.end(); ++it)
          {
            vector<String> mod_split;
            it->split(' ', mod_split);
            if (mod_split.size() < 2 || mod_split.size() > 3)
            {
              error(LOAD, String("Cannot parse fixed modification '") + *it + "'");
            }
            else
            {
              // C-term modification without specification or protein terminus
              if (mod_split[1] == "(C-term)" || (mod_split[1] == "(Protein" && mod_split[2] == "C-term)"))
              {
                temp_aa_sequence.setCTerminalModification(mod_split[0]);
              }
              // N-term modification without specification or protein terminus
              else if (mod_split[1] == "(N-term)" || (mod_split[1] == "(Protein" && mod_split[2] == "N-term)"))
              {
                temp_aa_sequence.setNTerminalModification(mod_split[0]);
              }
              // C-term modification for specific amino acid; e.g. <Modification> (N-term C)
              else if ((mod_split[1] == "(C-term") && (mod_split.size() == 3))
              {
                if ((temp_aa_sequence.end() - 1)->getOneLetterCode() == mod_split[2].remove(')'))
                {
                  temp_aa_sequence.setCTerminalModification(mod_split[0]);
                }
              }
              // N-term modification for specific amino acid; e.g. <Modification> (N-term C)
              else if ((mod_split[1] == "(N-term") && (mod_split.size() == 3))
              {
                if (temp_aa_sequence.begin()->getOneLetterCode() == mod_split[2].remove(')'))
                {
                  temp_aa_sequence.setNTerminalModification(mod_split[0]);
                }
              }
              else 
              { // e.g. Carboxymethyl (C)
                String AA = mod_split[1];
                AA.remove(')');
                AA.remove('(');
                for (Size i = 0; i != temp_aa_sequence.size(); ++i)
                {
                  if (AA == temp_aa_sequence[i].getOneLetterCode())
                  {
                    temp_aa_sequence.setModification(i, mod_split[0]);
                  }
                }
              }
            }
          }
        }
        actual_peptide_hit_.setSequence(temp_aa_sequence);
      }
      else if (tag_ == "pep_res_before")
      {
        String temp_string = character_buffer_.trim();
        if (!temp_string.empty())
        {
          actual_peptide_evidence_.setAABefore(temp_string[0]);
        }
      }
      else if (tag_ == "pep_res_after")
      {
        String temp_string = character_buffer_.trim();
        if (!temp_string.empty())
        {
          actual_peptide_evidence_.setAAAfter(temp_string[0]);
        }
      }
      else if (tag_ == "pep_var_mod_pos")
      {
        AASequence temp_aa_sequence = actual_peptide_hit_.getSequence();
        String temp_string = character_buffer_.trim();
        vector<String> parts;
        
        // E.g. seq: QKAAGSK, pos: 4.0000000.0 -> mod at position 4 in var_mods vector is n-terminal
        // therefore it is not possible to split Phospho (ST) to Phospho (S), Phospho (T) before this,
        // because the original order is required
        temp_string.split('.', parts);
        if (parts.size() == 3)
        {
          // handle internal modifications
          temp_string = parts[1];
          for (Size i = 0; i < temp_string.size(); ++i)
          {
            if (temp_string[i] != '0')
            {
              UInt temp_modification_index = String(temp_string[i]).toInt() - 1;
              OPENMS_PRECONDITION(temp_modification_index < search_parameters_.variable_modifications.size(), "Error when parsing variable modification string in <pep_var_mod_pos> (index too large)!");
              String& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);

              // e.g. "Carboxymethyl (C)"
              vector<String> mod_split;
              temp_modification.split(' ', mod_split);

              if (mod_split.size() >= 2)
              {
                // search this mod, if not directly use a general one
                temp_aa_sequence.setModification(i, mod_split[0]);
              }
              else
              {
                error(LOAD, String("Cannot parse variable modification '") + temp_modification  + "'");
              }
            }
          }

          temp_string = parts[0]; // N-term
          if (temp_string[0] != '0')
          {
            UInt temp_modification_index = String(temp_string[0]).toInt() - 1;
            String& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);
            vector<String> mod_split;
            temp_modification.split(' ', mod_split);

            if (mod_split.size() >= 2)
            {
              temp_aa_sequence.setNTerminalModification(mod_split[0]);
            }
            else
            {
              error(LOAD, String("Cannot parse variable N-term modification '") + temp_modification  + "'");
            }
          }
          temp_string = parts[2]; // C-term
          if (temp_string[0] != '0')
          {
            UInt temp_modification_index = String(temp_string[0]).toInt() - 1;
            String& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);
            vector<String> mod_split;
            temp_modification.split(' ', mod_split);

            if (mod_split.size() >= 2)
            {
              temp_aa_sequence.setCTerminalModification(mod_split[0]);
            }
            else
            {
              error(LOAD, String("Cannot parse variable C-term modification '") + temp_modification  + "'");
            }
          }

          actual_peptide_hit_.setSequence(temp_aa_sequence);
        }
      }
      else if (tag_ == "Date")
      {
        vector<String> parts;

        character_buffer_.trim().split('T', parts);
        if (parts.size() == 2)
        {
          date_.set(parts[0] + ' ' + parts[1].prefix('Z'));
          date_time_string_ = parts[0] + ' ' + parts[1].prefix('Z');
          identifier_ = "Mascot_" + date_time_string_;
        }
        protein_identification_.setDateTime(date_);
      }
      else if (tag_ == "StringTitle")
      {
        String title = character_buffer_.trim();
        vector<String> parts;

        actual_title_ = title;
        if (modified_peptides_.find(title) != modified_peptides_.end())
        {
          vector<AASequence>& temp_hits = modified_peptides_[title];
          vector<PeptideHit> temp_peptide_hits = id_data_[actual_query_ - 1].getHits();

          if (temp_hits.size() != temp_peptide_hits.size())
          {
            warning(LOAD, "pepXML hits and Mascot hits are not the same");
          }

          // pepXML can contain more hits than MascotXML; hence we try to match all of them...
          // run-time is O(n^2) in the number of peptide hits; should be a very small number

          for (Size i = 0; i < temp_peptide_hits.size(); ++i)
          {
            for (Size j = 0; j < temp_hits.size(); ++j)
            {
              if (temp_hits[j].isModified() && temp_hits[j].toUnmodifiedString() == temp_peptide_hits[i].getSequence().toUnmodifiedString())
              {
                temp_peptide_hits[i].setSequence(temp_hits[j]);
                break;
              }
            }
          }
          id_data_[actual_query_ - 1].setHits(temp_peptide_hits);
        }
        if (!id_data_[actual_query_ - 1].hasRT())
        {
          title.split('_', parts);
          if (parts.size() == 2)
          {
            id_data_[actual_query_ - 1].setRT(parts[1].toDouble());
          }
        }
      }
      else if (tag_ == "RTINSECONDS")
      {
        id_data_[actual_query_ - 1].setRT(character_buffer_.trim().toDouble());
      }
      else if (tag_ == "MascotVer")
      {
        protein_identification_.setSearchEngineVersion(character_buffer_.trim());
      }
      else if (tag_ == "DB")
      {
        search_parameters_.db = (character_buffer_.trim());
      }
      else if (tag_ == "FastaVer")
      {
        search_parameters_.db_version = (character_buffer_.trim());
      }
      else if (tag_ == "TAXONOMY")
      {
        search_parameters_.taxonomy = (character_buffer_.trim());
      }
      else if (tag_ == "CHARGE")
      {
        search_parameters_.charges = (character_buffer_.trim());
      }
      else if (tag_ == "PFA")
      {
        search_parameters_.missed_cleavages = character_buffer_.trim().toInt();
      }
      else if (tag_ == "MASS")
      {
        String temp_string = (character_buffer_.trim());
        if (temp_string == "Monoisotopic")
        {
          search_parameters_.mass_type = ProteinIdentification::MONOISOTOPIC;
        }
        else if (temp_string == "Average")
        {
          search_parameters_.mass_type = ProteinIdentification::AVERAGE;
        }
      }
      else if (tag_ == "MODS")
      {
        // if the modifications are listed in the "fixed_mods" section,
        // read from there; if <fixed_mods> was present it was already read
        if (search_parameters_.fixed_modifications.empty())
        {
          String temp_string = (character_buffer_.trim());
          vector<String> tmp_mods;
          temp_string.split(',', tmp_mods);
          
          for (vector<String>::const_iterator it = tmp_mods.begin(); it != tmp_mods.end(); ++it)
          {
            // check if modification is not on the remove list
            if (std::find(remove_fixed_mods_.begin(), remove_fixed_mods_.end(), *it) == remove_fixed_mods_.end())
            {
              // split because e.g. Phospho (ST)
              vector<String> mods_split = splitModificationBySpecifiedAA(*it);
              search_parameters_.fixed_modifications.insert(search_parameters_.fixed_modifications.end(), mods_split.begin(), mods_split.end());
            }
          }
        }
      }
      else if (tag_ == "IT_MODS")
      {
        // if the modifications are listed in the "variable_mods" section,
        // read from there, because sometimes mods are forced to be variable
        // (from user set fixed); if <variable_mods> was present it was already
        // read
        if (search_parameters_.variable_modifications.empty())
        {
          String temp_string = (character_buffer_.trim());
          vector<String> tmp_mods;
          temp_string.split(',', tmp_mods);
          
          for (vector<String>::const_iterator it = tmp_mods.begin(); it != tmp_mods.end(); ++it)
          {
            // split because e.g. Phospho (ST)
            vector<String> mods_split = splitModificationBySpecifiedAA(*it);
            search_parameters_.variable_modifications.insert(search_parameters_.variable_modifications.end(), mods_split.begin(), mods_split.end());
          }
        }
      }
      else if (tag_ == "CLE")
      {
        String temp_string = (character_buffer_.trim());
        if (ProteaseDB::getInstance()->hasEnzyme(temp_string))
        {
          search_parameters_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(temp_string));
        }
      }
      else if (tag_ == "TOL")
      {
        search_parameters_.precursor_mass_tolerance = (character_buffer_.trim()).toDouble();
      }
      else if (tag_ == "ITOL")
      {
        search_parameters_.fragment_mass_tolerance = (character_buffer_.trim()).toDouble();
      }
      else if (tag_ == "TOLU")
      {
        search_parameters_.precursor_mass_tolerance_ppm = (character_buffer_.trim()) == "ppm";
      }
      else if (tag_ == "ITOLU")
      {
        search_parameters_.fragment_mass_tolerance_ppm = (character_buffer_.trim()) == "ppm";
      }
      else if (tag_ == "name")
      {
        // cerr << "name tag: " << character_buffer_.trim() << "\n";
        if ((major_version_ == "1")
            // new since Mascot XML version 2.1 (at least): <fixed_mods> also have a subtag called <name>, thus we need to ensure we are in <variable_mods>
           || (tags_open_.size() >= 2 &&
               tags_open_[tags_open_.size() - 2] == "variable_mods"))
        {
          // e.g. Phospho (ST) cannot be split for variable modifications at this point, because the order of
          // variable modifications needs to be preserved. Split before search parameters are set.
          search_parameters_.variable_modifications.push_back(character_buffer_.trim());
          // cerr << "var. mod. added: " << search_parameters_.variable_modifications.back() << "\n";
        }
        else if (tags_open_.size() >= 2 &&
                 tags_open_[tags_open_.size() - 2] == "fixed_mods")
        {
          // check if modification is not on the remove list
          String fixed_mod = character_buffer_.trim();
          if (std::find(remove_fixed_mods_.begin(), remove_fixed_mods_.end(), fixed_mod) == remove_fixed_mods_.end())
          {
            // split because e.g. Phospho (ST)
            vector<String> mods_split = splitModificationBySpecifiedAA(character_buffer_.trim());
            search_parameters_.fixed_modifications.insert(search_parameters_.fixed_modifications.end(), mods_split.begin(), mods_split.end());
            // cerr << "fixed mod. added: " << search_parameters_.fixed_modifications.back() << "\n";
          }
          else
          {
            warning(LOAD, String("Modification removed as fixed modification: '") + character_buffer_.trim() + String("'"));
          }
        }
      }
      else if (tag_ == "warning")
      {
        warning(LOAD, String("Warnings were present: '") + character_buffer_ + String("'"));
        
        // check if fixed modification can only be used as variable modification
        if (character_buffer_.trim().hasSubstring("can only be used as a variable modification"))
        {
          vector<String> warn_split;
          character_buffer_.trim().split(';', warn_split);
          if (warn_split[0].hasPrefix("'"))
          {
            Size end_pos = warn_split[0].find("'", 1);
            if (end_pos < warn_split[0].size())
            {
              warn_split[0] = warn_split[0].substr(1, end_pos - 1);
              remove_fixed_mods_.push_back(warn_split[0]);
            }
          }
        }
      }
      else if (tag_ == "protein")
      {
        protein_identification_.setScoreType("Mascot");
        protein_identification_.insertHit(actual_protein_hit_);
        actual_protein_hit_ = ProteinHit();
      }
      else if (tag_ == "peptide")
      {
        bool already_stored(false);

        vector<PeptideHit> temp_peptide_hits = id_data_[peptide_identification_index_].getHits();

        vector<PeptideHit>::iterator it = temp_peptide_hits.begin();
        while (it != temp_peptide_hits.end())
        {
          if (it->getSequence() == actual_peptide_hit_.getSequence())
          {
            already_stored = true;
            break;
          }
          ++it;
        }

        if (!already_stored)
        {
          id_data_[peptide_identification_index_].setIdentifier(identifier_);
          id_data_[peptide_identification_index_].setScoreType("Mascot");
          actual_peptide_evidence_.setProteinAccession(actual_protein_hit_.getAccession());
          actual_peptide_hit_.addPeptideEvidence(actual_peptide_evidence_);
          id_data_[peptide_identification_index_].insertHit(actual_peptide_hit_);
        }
        else
        {
          actual_peptide_evidence_.setProteinAccession(actual_protein_hit_.getAccession());
          it->addPeptideEvidence(actual_peptide_evidence_);
          id_data_[peptide_identification_index_].setHits(temp_peptide_hits);
        }
        actual_peptide_evidence_ = PeptideEvidence();
        actual_peptide_hit_ = PeptideHit();
      }
      else if (tag_ == "u_peptide" || tag_ == "q_peptide")
      {
        id_data_[peptide_identification_index_].setIdentifier(identifier_);
        id_data_[peptide_identification_index_].setScoreType("Mascot");
        id_data_[peptide_identification_index_].insertHit(actual_peptide_hit_);
        actual_peptide_evidence_ = PeptideEvidence();
        actual_peptide_hit_ = PeptideHit();
      }
      else if (tag_ == "mascot_search_results")
      {
        protein_identification_.setSearchEngine("Mascot");
        protein_identification_.setIdentifier(identifier_);
        
        // split variable modifications e.g. Phospho (ST)
        //vector<String> var_mods;
        vector<String> var_mods;
        for (vector<String>::iterator it = search_parameters_.variable_modifications.begin(); it != search_parameters_.variable_modifications.end(); ++it)
        {
          vector<String> mods_split = splitModificationBySpecifiedAA(*it);
          var_mods.insert(var_mods.end(), mods_split.begin(), mods_split.end());
        }
        search_parameters_.variable_modifications = var_mods;
        protein_identification_.setSearchParameters(search_parameters_);        
      }

      tag_ = ""; // reset tag, for the following characters() call (due to line break) of the parent tag
      character_buffer_.clear();
    }

    void MascotXMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
    {
      // do not care about chars after internal tags, e.g.
      // <header>
      //   <COM>OpenMS_search</COM>
      //   <Date>
      // will trigger a characters() between </COM> and <Date>, which should be ignored
      if (tag_.empty())
      {
        return;
      }
      character_buffer_ += String(sm_.convert(chars));
    }
    
    vector<String> MascotXMLHandler::splitModificationBySpecifiedAA(const String& mod)
    {
      vector<String> mods;
      vector<String> parts;
      mod.split(' ', parts);
      
      // either format "Modification (Protein C-term)" or "Modification (C-term X)"
      if (parts.size() != 2)
      {
        mods.push_back(mod);
        return mods;
      }
      
      if (parts[1].hasPrefix("(N-term") || parts[1].hasPrefix("(C-term"))
      {
        mods.push_back(mod);
        return mods;
      }
      
      // format e.g. Phospho (ST)
      ModificationsDB* mod_db = ModificationsDB::getInstance();
      String AAs = parts[1];
      AAs.remove(')');
      AAs.remove('(');
      for (String::const_iterator it = AAs.begin(); it != AAs.end(); ++it)
      {
        String tmp_mod = parts[0] + " (" + *it + ")";
        if (mod_db->has(tmp_mod))
        {
          mods.push_back(tmp_mod);
        }
        else
        {
          throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tmp_mod);
        }
      }
      return mods;
    }
} // namespace OpenMS   // namespace Internal
