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

    MascotXMLHandler::MascotXMLHandler(ProteinIdentification& protein_identification, PeptideIdentificationList& id_data, const std::string& filename, map<std::string, vector<AASequence> >& modified_peptides, const SpectrumMetaDataLookup& lookup):
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

      tag_ =std::string(sm_.convert(qname));
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
        std::string attribute_value = attributeAsString_(attributes, s_protein_accession);
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
      tag_ =StringUtils::trimmed(std::string(sm_.convert(qname)));
      // cerr << "close: " << tag_ << endl;

      if (tags_open_.empty())
      {
        fatalError(LOAD,std::string("Closing tag ") + tag_ + " not matched by opening tag", __LINE__);
      }

      tags_open_.pop_back();

      if (tag_ == "NumQueries")
      {
        id_data_.resize(StringUtils::toInt32(StringUtils::trim(character_buffer_)));
      }
      else if (tag_ == "prot_score")
      {
        actual_protein_hit_.setScore(StringUtils::toInt32(StringUtils::trim(character_buffer_)));
      }
      else if (tag_ == "pep_exp_mz")
      {
        id_data_[peptide_identification_index_].setMZ(
          StringUtils::toDouble(StringUtils::trimmed(character_buffer_)));
      }
      else if (tag_ == "pep_scan_title")
      {
        // extract RT (and possibly m/z, if not already set) from title:
        std::string title = StringUtils::trim(character_buffer_);
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
          std::string msg = "<pep_scan_title> element has unexpected format '" + title + "'. Could not extract spectrum meta data.";
          error(LOAD, msg);
        }
        // did it work?
        if (!id_data_[peptide_identification_index_].getRT())
        {
          if (!no_rt_error_) // report the error only the first time
          {
            std::string msg = "Could not extract RT value ";
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
        actual_peptide_hit_.setCharge(StringUtils::toInt32(StringUtils::trim(character_buffer_)));
      }
      else if (tag_ == "pep_score")
      {
        actual_peptide_hit_.setScore(StringUtils::toDouble(StringUtils::trimmed(character_buffer_)));
      }
      else if (tag_ == "pep_expect")
      {
        // @todo what E-value flag? (andreas)
        actual_peptide_hit_.metaRegistry().registerName("EValue", "E-value of e.g. Mascot searches", "");
        actual_peptide_hit_.setMetaValue("EValue", StringUtils::toDouble(StringUtils::trimmed(character_buffer_)));
      }
      else if (tag_ == "pep_homol")
      {
        id_data_[peptide_identification_index_].setSignificanceThreshold(StringUtils::toDouble(StringUtils::trimmed(character_buffer_)));
      }
      else if (tag_ == "pep_ident")
      {
        double temp_homology = 0;
        double temp_identity = 0;

        // According to Matrix Science the homology threshold is only used if it
        // exists and is smaller than the identity threshold.
        temp_homology = id_data_[peptide_identification_index_].getSignificanceThreshold();
        temp_identity = StringUtils::toDouble(StringUtils::trimmed(character_buffer_));
        actual_peptide_hit_.setMetaValue("homology_threshold", temp_homology);
        actual_peptide_hit_.setMetaValue("identity_threshold", temp_identity);
        if (temp_homology > temp_identity || temp_homology == 0)
        {
          id_data_[peptide_identification_index_].setSignificanceThreshold(temp_identity);
        }
      }
      else if (tag_ == "pep_seq")
      {
        AASequence temp_aa_sequence = AASequence::fromString(StringUtils::trim(character_buffer_));
        
        // if everything is just read from the MascotXML file
        if (modified_peptides_.empty())
        {
          // fixed modifications
          for (vector<std::string>::const_iterator it = search_parameters_.fixed_modifications.begin(); it != search_parameters_.fixed_modifications.end(); ++it)
          {
            vector<std::string> mod_split;
            StringUtils::split(*it, ' ', mod_split);
            if (mod_split.size() < 2 || mod_split.size() > 3)
            {
              error(LOAD,std::string("Cannot parse fixed modification '") + *it + "'");
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
                if ((temp_aa_sequence.end() - 1)->getOneLetterCode() == StringUtils::remove(mod_split[2], ')'))
                {
                  temp_aa_sequence.setCTerminalModification(mod_split[0]);
                }
              }
              // N-term modification for specific amino acid; e.g. <Modification> (N-term C)
              else if ((mod_split[1] == "(N-term") && (mod_split.size() == 3))
              {
                if (temp_aa_sequence.begin()->getOneLetterCode() == StringUtils::remove(mod_split[2], ')'))
                {
                  temp_aa_sequence.setNTerminalModification(mod_split[0]);
                }
              }
              else 
              { // e.g. Carboxymethyl (C)
                std::string AA = mod_split[1];
                StringUtils::remove(AA, ')');
                StringUtils::remove(AA, '(');
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
        std::string temp_string = StringUtils::trim(character_buffer_);
        if (!temp_string.empty())
        {
          actual_peptide_evidence_.setAABefore(temp_string[0]);
        }
      }
      else if (tag_ == "pep_res_after")
      {
        std::string temp_string = StringUtils::trim(character_buffer_);
        if (!temp_string.empty())
        {
          actual_peptide_evidence_.setAAAfter(temp_string[0]);
        }
      }
      else if (tag_ == "pep_var_mod_pos")
      {
        AASequence temp_aa_sequence = actual_peptide_hit_.getSequence();
        std::string temp_string = StringUtils::trim(character_buffer_);
        vector<std::string> parts;
        
        // E.g. seq: QKAAGSK, pos: 4.0000000.0 -> mod at position 4 in var_mods vector is n-terminal
        // therefore it is not possible to split Phospho (ST) to Phospho (S), Phospho (T) before this,
        // because the original order is required
        StringUtils::split(temp_string, '.', parts);
        if (parts.size() == 3)
        {
          // handle internal modifications
          temp_string = parts[1];
          for (Size i = 0; i < temp_string.size(); ++i)
          {
            if (temp_string[i] != '0')
            {
              UInt temp_modification_index =StringUtils::toInt32(StringUtils::toStr(temp_string[i])) - 1;
              OPENMS_PRECONDITION(temp_modification_index < search_parameters_.variable_modifications.size(), "Error when parsing variable modification string in <pep_var_mod_pos> (index too large)!");
              std::string& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);

              // e.g. "Carboxymethyl (C)"
              vector<std::string> mod_split;
              StringUtils::split(temp_modification, ' ', mod_split);

              if (mod_split.size() >= 2)
              {
                // search this mod, if not directly use a general one
                temp_aa_sequence.setModification(i, mod_split[0]);
              }
              else
              {
                error(LOAD,std::string("Cannot parse variable modification '") + temp_modification  + "'");
              }
            }
          }

          temp_string = parts[0]; // N-term
          if (temp_string[0] != '0')
          {
            UInt temp_modification_index =StringUtils::toInt32(StringUtils::toStr(temp_string[0])) - 1;
            std::string& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);
            vector<std::string> mod_split;
            StringUtils::split(temp_modification, ' ', mod_split);

            if (mod_split.size() >= 2)
            {
              temp_aa_sequence.setNTerminalModification(mod_split[0]);
            }
            else
            {
              error(LOAD,std::string("Cannot parse variable N-term modification '") + temp_modification  + "'");
            }
          }
          temp_string = parts[2]; // C-term
          if (temp_string[0] != '0')
          {
            UInt temp_modification_index =StringUtils::toInt32(StringUtils::toStr(temp_string[0])) - 1;
            std::string& temp_modification = search_parameters_.variable_modifications.at(temp_modification_index);
            vector<std::string> mod_split;
            StringUtils::split(temp_modification, ' ', mod_split);

            if (mod_split.size() >= 2)
            {
              temp_aa_sequence.setCTerminalModification(mod_split[0]);
            }
            else
            {
              error(LOAD,std::string("Cannot parse variable C-term modification '") + temp_modification  + "'");
            }
          }

          actual_peptide_hit_.setSequence(temp_aa_sequence);
        }
      }
      else if (tag_ == "Date")
      {
        vector<std::string> parts;

        StringUtils::split(StringUtils::trimmed(character_buffer_), 'T', parts);
        if (parts.size() == 2)
        {
          date_.set(parts[0] + ' ' + StringUtils::prefix(parts[1], 'Z'));
          date_time_string_ = parts[0] + ' ' + StringUtils::prefix(parts[1], 'Z');
          identifier_ = "Mascot_" + date_time_string_;
        }
        protein_identification_.setDateTime(date_);
      }
      else if (tag_ == "StringTitle")
      {
        std::string title = StringUtils::trim(character_buffer_);
        vector<std::string> parts;

        actual_title_ = title;
        if (modified_peptides_.contains(title))
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
          StringUtils::split(title, '_', parts);
          if (parts.size() == 2)
          {
            id_data_[actual_query_ - 1].setRT(StringUtils::toDouble(parts[1]));
          }
        }
      }
      else if (tag_ == "RTINSECONDS")
      {
        id_data_[actual_query_ - 1].setRT(StringUtils::toDouble(StringUtils::trimmed(character_buffer_)));
      }
      else if (tag_ == "MascotVer")
      {
        protein_identification_.setSearchEngineVersion(StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "DB")
      {
        search_parameters_.db = (StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "FastaVer")
      {
        search_parameters_.db_version = (StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "TAXONOMY")
      {
        search_parameters_.taxonomy = (StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "CHARGE")
      {
        search_parameters_.charges = (StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "PFA")
      {
        search_parameters_.missed_cleavages = StringUtils::toInt32(StringUtils::trim(character_buffer_));
      }
      else if (tag_ == "MASS")
      {
        std::string temp_string = (StringUtils::trim(character_buffer_));
        if (temp_string == "Monoisotopic")
        {
          search_parameters_.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
        }
        else if (temp_string == "Average")
        {
          search_parameters_.mass_type = ProteinIdentification::PeakMassType::AVERAGE;
        }
      }
      else if (tag_ == "MODS")
      {
        // if the modifications are listed in the "fixed_mods" section,
        // read from there; if <fixed_mods> was present it was already read
        if (search_parameters_.fixed_modifications.empty())
        {
          std::string temp_string = (StringUtils::trim(character_buffer_));
          vector<std::string> tmp_mods;
          StringUtils::split(temp_string, ',', tmp_mods);
          
          for (vector<std::string>::const_iterator it = tmp_mods.begin(); it != tmp_mods.end(); ++it)
          {
            // check if modification is not on the remove list
            if (std::find(remove_fixed_mods_.begin(), remove_fixed_mods_.end(), *it) == remove_fixed_mods_.end())
            {
              // split because e.g. Phospho (ST)
              vector<std::string> mods_split = splitModificationBySpecifiedAA(*it);
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
          std::string temp_string = (StringUtils::trim(character_buffer_));
          vector<std::string> tmp_mods;
          StringUtils::split(temp_string, ',', tmp_mods);
          
          for (vector<std::string>::const_iterator it = tmp_mods.begin(); it != tmp_mods.end(); ++it)
          {
            // split because e.g. Phospho (ST)
            vector<std::string> mods_split = splitModificationBySpecifiedAA(*it);
            search_parameters_.variable_modifications.insert(search_parameters_.variable_modifications.end(), mods_split.begin(), mods_split.end());
          }
        }
      }
      else if (tag_ == "CLE")
      {
        std::string temp_string = (StringUtils::trim(character_buffer_));
        if (ProteaseDB::getInstance()->hasEnzyme(temp_string))
        {
          search_parameters_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(temp_string));
        }
      }
      else if (tag_ == "TOL")
      {
        search_parameters_.precursor_mass_tolerance = StringUtils::toDouble((StringUtils::trim(character_buffer_)));
      }
      else if (tag_ == "ITOL")
      {
        search_parameters_.fragment_mass_tolerance = StringUtils::toDouble((StringUtils::trim(character_buffer_)));
      }
      else if (tag_ == "TOLU")
      {
        search_parameters_.precursor_mass_tolerance_ppm = (StringUtils::trim(character_buffer_)) == "ppm";
      }
      else if (tag_ == "ITOLU")
      {
        search_parameters_.fragment_mass_tolerance_ppm = (StringUtils::trim(character_buffer_)) == "ppm";
      }
      else if (tag_ == "name")
      {
        // cerr << "name tag: " << StringUtils::trim(character_buffer_) << "\n";
        if ((major_version_ == "1")
            // new since Mascot XML version 2.1 (at least): <fixed_mods> also have a subtag called <name>, thus we need to ensure we are in <variable_mods>
           || (tags_open_.size() >= 2 &&
               tags_open_[tags_open_.size() - 2] == "variable_mods"))
        {
          // e.g. Phospho (ST) cannot be split for variable modifications at this point, because the order of
          // variable modifications needs to be preserved. Split before search parameters are set.
          search_parameters_.variable_modifications.push_back(StringUtils::trim(character_buffer_));
          // cerr << "var. mod. added: " << search_parameters_.variable_modifications.back() << "\n";
        }
        else if (tags_open_.size() >= 2 &&
                 tags_open_[tags_open_.size() - 2] == "fixed_mods")
        {
          // check if modification is not on the remove list
          std::string fixed_mod = StringUtils::trim(character_buffer_);
          if (std::find(remove_fixed_mods_.begin(), remove_fixed_mods_.end(), fixed_mod) == remove_fixed_mods_.end())
          {
            // split because e.g. Phospho (ST)
            vector<std::string> mods_split = splitModificationBySpecifiedAA(StringUtils::trim(character_buffer_));
            search_parameters_.fixed_modifications.insert(search_parameters_.fixed_modifications.end(), mods_split.begin(), mods_split.end());
            // cerr << "fixed mod. added: " << search_parameters_.fixed_modifications.back() << "\n";
          }
          else
          {
            warning(LOAD,std::string("Modification removed as fixed modification: '") + StringUtils::trim(character_buffer_) + std::string("'"));
          }
        }
      }
      else if (tag_ == "warning")
      {
        warning(LOAD,std::string("Warnings were present: '") + character_buffer_ + std::string("'"));
        
        // check if fixed modification can only be used as variable modification
        if (StringUtils::hasSubstring(StringUtils::trimmed(character_buffer_), "can only be used as a variable modification"))
        {
          vector<std::string> warn_split;
          StringUtils::split(StringUtils::trimmed(character_buffer_), ';', warn_split);
          if (StringUtils::hasPrefix(warn_split[0], "'"))
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
        //vector<std::string> var_mods;
        vector<std::string> var_mods;
        for (vector<std::string>::iterator it = search_parameters_.variable_modifications.begin(); it != search_parameters_.variable_modifications.end(); ++it)
        {
          vector<std::string> mods_split = splitModificationBySpecifiedAA(*it);
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
      character_buffer_ +=std::string(sm_.convert(chars));
    }
    
    vector<std::string> MascotXMLHandler::splitModificationBySpecifiedAA(const std::string& mod)
    {
      vector<std::string> mods;
      vector<std::string> parts;
      StringUtils::split(mod, ' ', parts);
      
      // either format "Modification (Protein C-term)" or "Modification (C-term X)"
      if (parts.size() != 2)
      {
        mods.push_back(mod);
        return mods;
      }
      
      if (StringUtils::hasPrefix(parts[1], "(N-term") || StringUtils::hasPrefix(parts[1], "(C-term"))
      {
        mods.push_back(mod);
        return mods;
      }
      
      // format e.g. Phospho (ST)
      const ModificationsDB* mod_db = ModificationsDB::getInstance();
      std::string AAs = parts[1];
      StringUtils::remove(AAs, ')');
      StringUtils::remove(AAs, '(');
      for (std::string::const_iterator it = AAs.begin(); it != AAs.end(); ++it)
      {
        std::string tmp_mod = parts[0] + " (" + *it + ")";
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
