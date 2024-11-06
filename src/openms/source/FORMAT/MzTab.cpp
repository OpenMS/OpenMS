// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzTab.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>


#include <tuple>

using namespace std;

namespace OpenMS
{
  MzTabModification::MzTabModification() = default;

  bool MzTabModification::isNull() const
  {
    return pos_param_pairs_.empty() && mod_identifier_.isNull();
  }

  void MzTabModification::setNull(bool b)
  {
    if (b)
    {
      pos_param_pairs_.clear();
      mod_identifier_.setNull(true);
    }
  }

  void MzTabModification::setPositionsAndParameters(const std::vector<std::pair<Size, MzTabParameter> >& ppp)
  {
    pos_param_pairs_ = ppp;
  }

  std::vector<std::pair<Size, MzTabParameter> > MzTabModification::getPositionsAndParameters() const
  {
    return pos_param_pairs_;
  }

  void MzTabModification::setModificationIdentifier(const MzTabString& mod_id)
  {
    mod_identifier_ = mod_id;
  }

  MzTabString MzTabModification::getModOrSubstIdentifier() const
  {
    assert(!isNull());
    return mod_identifier_;
  }

  String MzTabModification::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      String pos_param_string;

      for (Size i = 0; i != pos_param_pairs_.size(); ++i)
      {
        pos_param_string += pos_param_pairs_[i].first;

        // attach MzTabParameter if available
        if (!pos_param_pairs_[i].second.isNull())
        {
          pos_param_string += pos_param_pairs_[i].second.toCellString();
        }

        // add | as separator (except for last one)
        if (i < pos_param_pairs_.size() - 1)
        {
          pos_param_string += String("|");
        }
      }

      // quick sanity check
      if (mod_identifier_.isNull())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Modification or Substitution identifier MUST NOT be null or empty in MzTabModification"));
      }

      String res;
      // only add '-' if we have position information
      if (!pos_param_string.empty())
      {
        res = pos_param_string + "-" + mod_identifier_.toCellString();
      }
      else
      {
        res = mod_identifier_.toCellString();
      }
      return res;
    }
  }

  void MzTabModification::fromCellString(const String& s)
  {
    String lower = s;
    lower.toLower().trim();
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      if (!lower.hasSubstring("-")) // no positions? simply use s as mod identifier
      {
        mod_identifier_.set(String(s).trim());
      }
      else
      {
        String ss = s;
        ss.trim();
        std::vector<String> fields;
        ss.split("-", fields);

        if (fields.size() != 2)
        {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Can't convert to MzTabModification from '") + s);
        }
        mod_identifier_.fromCellString(fields[1].trim());

        std::vector<String> position_fields;
        fields[0].split("|", position_fields);

        for (Size i = 0; i != position_fields.size(); ++i)
        {
          Size spos = position_fields[i].find_first_of("[");

          if (spos == std::string::npos) // only position information and no parameter
          {
            pos_param_pairs_.emplace_back(position_fields[i].toInt(), MzTabParameter());
          }
          else
          {
            // extract position part
            Int pos = String(position_fields[i].begin(), position_fields[i].begin() + spos).toInt();

            // extract [,,,] part
            MzTabParameter param;
            param.fromCellString(position_fields[i].substr(spos));
            pos_param_pairs_.emplace_back(pos, param);
          }
        }
      }
    }
  }

  bool MzTabModificationList::isNull() const
  {
    return entries_.empty();
  }

  void MzTabModificationList::setNull(bool b)
  {
    if (b)
    {
      entries_.clear();
    }
  }

  String MzTabModificationList::toCellString() const
  {
    if (isNull())
    {
      return "null";
    }
    else
    {
      String ret;
      for (std::vector<MzTabModification>::const_iterator it = entries_.begin(); it != entries_.end(); ++it)
      {
        if (it != entries_.begin())
        {
          ret += ",";
        }
        ret += it->toCellString();
      }
      return ret;
    }
  }

  void MzTabModificationList::fromCellString(const String& s)
  {
    String lower = s;
    lower.toLower().trim();
    if (lower == "null")
    {
      setNull(true);
    }
    else
    {
      String ss = s;
      std::vector<String> fields;

      if (!ss.hasSubstring("[")) // no parameters
      {
        ss.split(",", fields);
        for (Size i = 0; i != fields.size(); ++i)
        {
          MzTabModification ms;
          ms.fromCellString(fields[i]);
          entries_.push_back(ms);
        }
      }
      else
      {
        // example string: 3|4[a,b,,v]|8[,,"blabla, [bla]",v],1|2|3[a,b,,v]-mod:123
        // we don't want to split at the , inside of [ ]  MzTabParameter brackets.
        // Additionally,  and we don't want to recognise quoted brackets inside the MzTabParameter where they can occur in quoted text (see example string)
        bool in_param_bracket = false;
        bool in_quotes = false;

        for (Size pos = 0; pos != ss.size(); ++pos)
        {
          // param_bracket state
          if (ss[pos] == '[' && !in_quotes)
          {
            in_param_bracket = true;
            continue;
          }

          if (ss[pos] == ']' && !in_quotes)
          {
            in_param_bracket = false;
            continue;
          }

          // quote state
          if (ss[pos] == '\"')
          {
            in_quotes = !in_quotes;
            continue;
          }

          // comma in param bracket
          if (ss[pos] == ',' && !in_quotes && in_param_bracket)
          {
            ss[pos] = ((char)007); // use ASCII bell as temporary separator
            continue;
          }
        }

        // now the split at comma is save
        ss.split(",", fields);

        for (Size i = 0; i != fields.size(); ++i)
        {
          fields[i].substitute(((char)007), ','); // resubstitute comma after split
          MzTabModification ms;
          ms.fromCellString(fields[i]);
          entries_.push_back(ms);
        }
      }
    }
  }

  std::vector<MzTabModification> MzTabModificationList::get() const
  {
    return entries_;
  }

  void MzTabModificationList::set(const std::vector<MzTabModification>& entries)
  {
    entries_ = entries;
  }

  MzTabProteinSectionRow::MzTabProteinSectionRow()
  {
    // use "," as list separator because "|" can be used for go terms and protein accessions
    go_terms.setSeparator(',');
    ambiguity_members.setSeparator(',');
  }

  MzTabMetaData::MzTabMetaData()
  {
    mz_tab_version.fromCellString(String("1.0.0"));
  }

  // static method remapping the target/decoy column from an opt_ to a standardized column
  static void remapTargetDecoyPSMAndPeptideSection_(std::vector<MzTabOptionalColumnEntry>& opt_entries)
  {
    const String old_header("opt_global_target_decoy");
    const String new_header("opt_global_cv_MS:1002217_decoy_peptide"); // for PRIDE
    for (auto &opt_entry : opt_entries)
    {
      if (opt_entry.first == old_header || opt_entry.first == new_header)
      {
	opt_entry.first = new_header;
        const String &current_value = opt_entry.second.get();
        if (current_value == "target" || current_value == "target+decoy")
        {
          opt_entry.second = MzTabString("0");
        }
        else if (current_value == "decoy")
        {
	  opt_entry.second = MzTabString("1");
        }
      }
    }
  }

  // static method remapping the target/decoy column from an opt_ to a standardized column
  static void remapTargetDecoyProteinSection_(std::vector<MzTabOptionalColumnEntry>& opt_entries)
  {
    const String old_header("opt_global_target_decoy");
    const String new_header("opt_global_cv_PRIDE:0000303_decoy_hit"); // for PRIDE
    for (auto &opt_entry : opt_entries)
    {
      if (opt_entry.first == old_header || opt_entry.first == new_header)
      {
	opt_entry.first = new_header;
        const String &current_value = opt_entry.second.get();
        if (current_value == "target" || current_value == "target+decoy")
        {
          opt_entry.second = MzTabString("0");
        }
        else if (current_value == "decoy")
        {
	  opt_entry.second = MzTabString("1");
        }
      }
    }
  }

  const MzTabMetaData& MzTab::getMetaData() const
  {
    return meta_data_;
  }

  void MzTab::setMetaData(const MzTabMetaData& md)
  {
    meta_data_ = md;
  }

  MzTabProteinSectionRows& MzTab::getProteinSectionRows()
  {
    return protein_data_;
  }

  const MzTabProteinSectionRows& MzTab::getProteinSectionRows() const
  {
    return protein_data_;
  }

  void MzTab::setProteinSectionRows(const MzTabProteinSectionRows& psd)
  {
    protein_data_ = psd;
  }

  MzTabPeptideSectionRows& MzTab::getPeptideSectionRows()
  {
    return peptide_data_;
  }

  const MzTabPeptideSectionRows& MzTab::getPeptideSectionRows() const
  {
    return peptide_data_;
  }

  void MzTab::setPeptideSectionRows(const MzTabPeptideSectionRows& psd)
  {
    peptide_data_ = psd;
  }

  const MzTabPSMSectionRows& MzTab::getPSMSectionRows() const
  {
    return psm_data_;
  }

  MzTabPSMSectionRows& MzTab::getPSMSectionRows()
  {
    return psm_data_;
  }

  size_t MzTab::getNumberOfPSMs() const
  {
    std::unordered_set<Int> psm_ids;
    for (const auto& psm : psm_data_)
    {
      psm_ids.insert(psm.PSM_ID.get());
    }
    return psm_ids.size();
  }

  void MzTab::setPSMSectionRows(const MzTabPSMSectionRows& psd)
  {
    psm_data_ = psd;
  }

  const MzTabSmallMoleculeSectionRows& MzTab::getSmallMoleculeSectionRows() const
  {
    return small_molecule_data_;
  }

  void MzTab::setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd)
  {
    small_molecule_data_ = smsd;
  }

  const MzTabNucleicAcidSectionRows& MzTab::getNucleicAcidSectionRows() const
  {
    return nucleic_acid_data_;
  }

  void MzTab::setNucleicAcidSectionRows(const MzTabNucleicAcidSectionRows& nasd)
  {
    nucleic_acid_data_ = nasd;
  }

  const MzTabOligonucleotideSectionRows& MzTab::getOligonucleotideSectionRows() const
  {
    return oligonucleotide_data_;
  }

  void MzTab::setOligonucleotideSectionRows(const MzTabOligonucleotideSectionRows& onsd)
  {
    oligonucleotide_data_ = onsd;
  }

  const MzTabOSMSectionRows& MzTab::getOSMSectionRows() const
  {
    return osm_data_;
  }

  void MzTab::setOSMSectionRows(const MzTabOSMSectionRows& osd)
  {
    osm_data_ = osd;
  }

  void MzTab::setCommentRows(const std::map<Size, String>& com)
  {
    comment_rows_ = com;
  }

  void MzTab::setEmptyRows(const std::vector<Size>& empty)
  {
    empty_rows_ = empty;
  }

  const std::vector<Size>& MzTab::getEmptyRows() const
  {
    return empty_rows_;
  }

  const std::map<Size, String>& MzTab::getCommentRows() const
  {
    return comment_rows_;
  }

  std::vector<String> MzTab::getProteinOptionalColumnNames() const
  {
    return getOptionalColumnNames_(protein_data_);
  }

  std::vector<String> MzTab::getPeptideOptionalColumnNames() const
  {
    return getOptionalColumnNames_(peptide_data_);
  }

  std::vector<String> MzTab::getPSMOptionalColumnNames() const
  {
    return getOptionalColumnNames_(psm_data_);
  }

  std::vector<String> MzTab::getSmallMoleculeOptionalColumnNames() const
  {
    return getOptionalColumnNames_(small_molecule_data_);
  }

  std::vector<String> MzTab::getNucleicAcidOptionalColumnNames() const
  {
    return getOptionalColumnNames_(nucleic_acid_data_);
  }

  std::vector<String> MzTab::getOligonucleotideOptionalColumnNames() const
  {
    return getOptionalColumnNames_(oligonucleotide_data_);
  }

  std::vector<String> MzTab::getOSMOptionalColumnNames() const
  {
    return getOptionalColumnNames_(osm_data_);
  }

  void MzTabPSMSectionRow::addPepEvidenceToRows(const vector<PeptideEvidence>& peptide_evidences)
  {
    if (peptide_evidences.empty())
    {
      // report without pep evidence information
      pre = MzTabString();
      post = MzTabString();
      start = MzTabString();
      end = MzTabString();
      return;
    }

    String pre, post, start, end, accession;
    for (Size i = 0; i != peptide_evidences.size(); ++i)
    {
      // get AABefore and AAAfter as well as start and end for all pep evidences

      // pre/post
      // from spec: Amino acid preceding the peptide (coming from the PSM) in the protein
      // sequence. If unknown “null” MUST be used, if the peptide is N-terminal “-“
      // MUST be used.
      if (peptide_evidences[i].getAABefore() == PeptideEvidence::UNKNOWN_AA)
      {
        pre += "null";
      }
      else if (peptide_evidences[i].getAABefore() == PeptideEvidence::N_TERMINAL_AA)
      {
        pre += "-";
      }
      else
      {
        pre += String(peptide_evidences[i].getAABefore());
      }

      if (peptide_evidences[i].getAAAfter() == PeptideEvidence::UNKNOWN_AA)
      {
        post += "null";
      }
      else if (peptide_evidences[i].getAAAfter() == PeptideEvidence::C_TERMINAL_AA)
      {
        post += "-";
      }
      else
      {
        post += String(peptide_evidences[i].getAAAfter());
      }

      // start/end
      if (peptide_evidences[i].getStart() == PeptideEvidence::UNKNOWN_POSITION)
      {
        start += "null";
      }
      else
      {
        start += String(peptide_evidences[i].getStart() + 1); // counting in mzTab starts at 1
      }

      if (peptide_evidences[i].getEnd() == PeptideEvidence::UNKNOWN_POSITION)
      {
        end += "null";
      }
      else
      {
        end += String(peptide_evidences[i].getEnd() + 1); // counting in mzTab starts at 1
      }

      accession += peptide_evidences[i].getProteinAccession();

      if (i < peptide_evidences.size() - 1) { pre += ','; post += ','; start += ','; end += ','; accession += ',';}
    }
    this->pre = MzTabString(pre);
    this->post = MzTabString(post);
    this->start = MzTabString(start);
    this->end = MzTabString(end);
    this->accession = MzTabString(accession);
  }


  void MzTab::addMetaInfoToOptionalColumns(
    const set<String>& keys,
    vector<MzTabOptionalColumnEntry>& opt,
    const String& id,
    const MetaInfoInterface& meta)
  {
    for (String const & key : keys)
    {
      MzTabOptionalColumnEntry opt_entry;
      // column names must not contain spaces
      opt_entry.first = "opt_" + id + "_" + String(key).substitute(' ','_');
      
      if (meta.metaValueExists(key))
      {
        opt_entry.second = MzTabString(meta.getMetaValue(key).toString());
      } // otherwise it is default ("null")

      opt.push_back(opt_entry);
    }
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromModifications(const vector<String>& mods)
  {
    map<Size, MzTabModificationMetaData> mods_mztab;
    Size index(1);
    for (String const & s : mods)
    {
      MzTabModificationMetaData mod;
      MzTabParameter mp;
      ModificationsDB* mod_db = ModificationsDB::getInstance();
      String unimod_accession;
      try
      {
        const ResidueModification* m = mod_db->getModification(s);
        unimod_accession = m->getUniModAccession();
        if (!unimod_accession.empty())
        {
          // MzTab standard is to report Unimod accession.
          mp.setCVLabel("UNIMOD");
          mp.setAccession(unimod_accession.toUpper());
        }
        mp.setName(m->getId());
        mod.modification = mp;

        if (m->getTermSpecificity() == ResidueModification::C_TERM)
        {
          mod.position = MzTabString("Any C-term");
        }
        else if (m->getTermSpecificity() == ResidueModification::N_TERM)
        {
          mod.position = MzTabString("Any N-term");
        }
        else if (m->getTermSpecificity() == ResidueModification::ANYWHERE)
        {
          mod.position = MzTabString("Anywhere");
        }
        else if (m->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
        {
          mod.position = MzTabString("Protein C-term");
        }
        else if (m->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
        {
          mod.position = MzTabString("Protein N-term");
        }
        mod.site = MzTabString(String(m->getOrigin()));
        mods_mztab[index] = mod;
      }
      catch(...)
      {
        OPENMS_LOG_WARN << "Skipping unknown residue modification: '" + s + "'" << endl; 
      }
      ++index;
    }
    return mods_mztab;
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromVariableModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002454, No variable modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  map<Size, MzTabModificationMetaData> MzTab::generateMzTabStringFromFixedModifications(const vector<String>& mods)
  {
    if (mods.empty())
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      MzTabModificationMetaData mod_mtd;
      mod_mtd.modification.fromCellString("[MS, MS:1002453, No fixed modifications searched, ]");
      mods_mztab.insert(make_pair(1, mod_mtd));
      return mods_mztab;
    }
    else
    {
      return generateMzTabStringFromModifications(mods);
    }
  }

  MzTab MzTab::exportFeatureMapToMzTab(
    const FeatureMap & feature_map,
    const String & filename)
  {
    OPENMS_LOG_INFO << "exporting feature map: \"" << filename << "\" to mzTab: " << std::endl;
    MzTab mztab;
    MzTabMetaData meta_data;

    const vector<ProteinIdentification> &prot_ids = feature_map.getProteinIdentifications();
    vector<String> var_mods, fixed_mods;
    MzTabString db, db_version;
    if (!prot_ids.empty())
    {
      const ProteinIdentification::SearchParameters &sp = prot_ids[0].getSearchParameters();
      var_mods = sp.variable_modifications;
      fixed_mods = sp.fixed_modifications;
      db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
      db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
    }

    meta_data.variable_mod = generateMzTabStringFromVariableModifications(var_mods);
    meta_data.fixed_mod = generateMzTabStringFromFixedModifications(fixed_mods);

    // mandatory meta values
    meta_data.mz_tab_type = MzTabString("Quantification");
    meta_data.mz_tab_mode = MzTabString("Summary");
    meta_data.description = MzTabString("OpenMS export from featureXML");

    MzTabMSRunMetaData ms_run;
    StringList spectra_data;
    feature_map.getPrimaryMSRunPath(spectra_data);

    if (!spectra_data.empty())
    {
      // prepend file:// if not there yet
      String m = spectra_data[0];
      if (!m.hasPrefix("file://")) {m = String("file://") + m; }
      ms_run.location = MzTabString(m);
    }
    else
    {
      ms_run.location = MzTabString();
    }

    meta_data.ms_run[1] = ms_run;
    meta_data.uri[1] = MzTabString(filename);
    meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO: we currently only support psm search engine scores annotated to the identification run
    meta_data.peptide_search_engine_score[1] = MzTabParameter();

    mztab.setMetaData(meta_data);

    // pre-analyze data for occurring meta values at feature and peptide hit level
    // these are used to build optional columns containing the meta values in internal data structures
    set<String> feature_user_value_keys;
    set<String> peptide_identifications_user_value_keys;
    set<String> peptide_hit_user_value_keys;
    MzTab::getFeatureMapMetaValues_(
      feature_map,     
      feature_user_value_keys,
      peptide_identifications_user_value_keys,
      peptide_hit_user_value_keys);

    for (Size i = 0; i < feature_map.size(); ++i)
    {
      const Feature& f = feature_map[i];
      auto row = peptideSectionRowFromFeature_(
        f,
        feature_user_value_keys, 
        peptide_identifications_user_value_keys, 
        peptide_hit_user_value_keys, 
        fixed_mods);
      mztab.getPeptideSectionRows().emplace_back(std::move(row));
    }

    return mztab;
  }

  MzTabPeptideSectionRow MzTab::peptideSectionRowFromFeature_(
    const Feature& f, 
    const set<String>& feature_user_value_keys,
    const set<String>& peptide_identifications_user_value_keys,
    const set<String>& peptide_hit_user_value_keys,
    const vector<String>& fixed_mods)
  {
    MzTabPeptideSectionRow row;
    row.mass_to_charge = MzTabDouble(f.getMZ());
    MzTabDoubleList rt_list;
    vector<MzTabDouble> rts;
    rts.emplace_back(f.getRT());
    rt_list.set(rts);
    row.retention_time = rt_list;

    // set rt window if a bounding box has been set
    vector<MzTabDouble> window;
    if (f.getConvexHull().getBoundingBox() != DBoundingBox<2>())
    {
      window.emplace_back(f.getConvexHull().getBoundingBox().minX());
      window.emplace_back(f.getConvexHull().getBoundingBox().maxX());
    }

    MzTabDoubleList rt_window;
    rt_window.set(window);
    row.retention_time_window = rt_window;
    row.charge = MzTabInteger(f.getCharge());
    row.peptide_abundance_stdev_study_variable[1];
    row.peptide_abundance_std_error_study_variable[1];
    row.peptide_abundance_study_variable[1] = MzTabDouble(f.getIntensity());
    row.best_search_engine_score[1] = MzTabDouble();
    row.search_engine_score_ms_run[1][1] = MzTabDouble();

    // create opt_ column for peptide sequence containing modification
    MzTabOptionalColumnEntry opt_global_modified_sequence;
    opt_global_modified_sequence.first = "opt_global_cv_MS:1000889_peptidoform_sequence";
    row.opt_.push_back(opt_global_modified_sequence);

    // create and fill opt_ columns for feature (peptide) user values
    addMetaInfoToOptionalColumns(feature_user_value_keys, row.opt_, String("global"), f);

    const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
    if (pep_ids.empty())
    {
      // still add empty opt_ columns before returning
      addMetaInfoToOptionalColumns(peptide_identifications_user_value_keys, row.opt_, "global", MetaInfoInterface());
      addMetaInfoToOptionalColumns(peptide_hit_user_value_keys, row.opt_, "global", MetaInfoInterface());      
      return row;
    }

    const PeptideIdentification& best_pid = f.getPeptideIdentifications()[0];
    addMetaInfoToOptionalColumns(peptide_identifications_user_value_keys, row.opt_, "global", best_pid);

    // TODO: here we assume that all have the same score type etc.
    vector<PeptideHit> all_hits;
    for (const PeptideIdentification& it : pep_ids)
    {
      all_hits.insert(all_hits.end(), it.getHits().begin(), it.getHits().end());
    }

    if (all_hits.empty())
    { 
      addMetaInfoToOptionalColumns(peptide_identifications_user_value_keys, row.opt_, "global", MetaInfoInterface());
      return row;
    }

    // create new peptide id object to assist in sorting
    PeptideIdentification new_pep_id = pep_ids[0];
    new_pep_id.setHits(all_hits);
    new_pep_id.assignRanks();

    const PeptideHit& best_ph = new_pep_id.getHits()[0];
    const AASequence& aas = best_ph.getSequence();
    row.sequence = MzTabString(aas.toUnmodifiedString());

    row.modifications = extractModificationList(best_ph, fixed_mods, vector<String>());

    const set<String>& accessions = best_ph.extractProteinAccessionsSet();
    const vector<PeptideEvidence>& peptide_evidences = best_ph.getPeptideEvidences();

    row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
    // select accession of first peptide_evidence as representative ("leading") accession
    row.accession = peptide_evidences.empty() ? MzTabString() : MzTabString(peptide_evidences[0].getProteinAccession());
    row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());
    row.search_engine_score_ms_run[1][1] = MzTabDouble(best_ph.getScore());

    // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
    for (Size j = 0; j != row.opt_.size(); ++j)
    {
      MzTabOptionalColumnEntry& opt_entry = row.opt_[j];

      if (opt_entry.first == "opt_global_cv_MS:1000889_peptidoform_sequence")
      {
        opt_entry.second = MzTabString(aas.toString());
      }
    }

    // create and fill opt_ columns for psm (PeptideHit) user values
    addMetaInfoToOptionalColumns(peptide_hit_user_value_keys, row.opt_, String("global"), best_ph);

    // remap the target/decoy column
    remapTargetDecoyPSMAndPeptideSection_(row.opt_);

    return row;
  }

  MzTabPeptideSectionRow MzTab::peptideSectionRowFromConsensusFeature_(
    const ConsensusFeature& c, 
    const ConsensusMap& consensus_map,
    const StringList& ms_runs,
    const Size n_study_variables,
    const set<String>& consensus_feature_user_value_keys,
    const set<String>& peptide_identifications_user_value_keys,
    const set<String>& peptide_hit_user_value_keys,
    const map<String, size_t>& idrun_2_run_index,
    const map<pair<size_t,size_t>,size_t>& map_run_fileidx_2_msfileidx,
    const std::map< std::pair< String, unsigned >, unsigned>& path_label_to_assay,
    const vector<String>& fixed_mods,
    bool export_subfeatures)
  {
    MzTabPeptideSectionRow row;

    const ConsensusMap::ColumnHeaders& cm_column_headers = consensus_map.getColumnHeaders();
    const String & experiment_type = consensus_map.getExperimentType();
    const vector<ProteinIdentification>& prot_id = consensus_map.getProteinIdentifications();

    // create opt_ column for peptide sequence containing modification
    MzTabOptionalColumnEntry opt_global_modified_sequence;
    opt_global_modified_sequence.first = "opt_global_cv_MS:1000889_peptidoform_sequence";
    row.opt_.push_back(opt_global_modified_sequence);

    // Defines how to consume user value keys for the upcoming keys
    const auto addUserValueToRowBy = [&row](const function<void(const String &s, MzTabOptionalColumnEntry &entry)>& f) -> function<void(const String &key)>
    {
      return [f,&row](const String &user_value_key)
        {
          MzTabOptionalColumnEntry opt_entry;
          opt_entry.first = "opt_global_" + user_value_key;
          f(user_value_key, opt_entry);

          // Use default column_header for target decoy
          row.opt_.push_back(opt_entry);
        };
    };

    // create opt_ columns for consensus map user values
    for_each(consensus_feature_user_value_keys.begin(), consensus_feature_user_value_keys.end(),
      addUserValueToRowBy([&c](const String &key, MzTabOptionalColumnEntry &opt_entry)
        {
          if (c.metaValueExists(key))
          {
            opt_entry.second = MzTabString(c.getMetaValue(key).toString());
          }
        })
      );

    // add optional columns for first peptide identification in consensus feature
    for_each(peptide_identifications_user_value_keys.begin(), peptide_identifications_user_value_keys.end(),
      addUserValueToRowBy([&c](const String &key, MzTabOptionalColumnEntry &opt_entry)
        {
          opt_entry.second = MzTabString(c.getMetaValue(key).toString());
        })
      );

    // create opt_ columns for psm (PeptideHit) user values
    for_each(peptide_hit_user_value_keys.begin(), peptide_hit_user_value_keys.end(),
      				addUserValueToRowBy([](const String&, MzTabOptionalColumnEntry&){}));

    row.mass_to_charge = MzTabDouble(c.getMZ());
    MzTabDoubleList rt_list;
    vector<MzTabDouble> rts;
    rts.emplace_back(c.getRT());
    rt_list.set(rts);
    row.retention_time = rt_list;
    MzTabDoubleList rt_window;
    row.retention_time_window = rt_window;
    row.charge = MzTabInteger(c.getCharge());
    row.best_search_engine_score[1] = MzTabDouble();

    // initialize columns
    OPENMS_LOG_DEBUG << "Initializing study variables:" << n_study_variables << endl;
    for (Size study_variable = 1; study_variable <= n_study_variables; ++study_variable)
    {
      row.peptide_abundance_stdev_study_variable[study_variable] = MzTabDouble();
      row.peptide_abundance_std_error_study_variable[study_variable] = MzTabDouble();
      row.peptide_abundance_study_variable[study_variable] = MzTabDouble();
    }

    for (Size ms_run = 1; ms_run <= ms_runs.size(); ++ms_run)
    {
      row.search_engine_score_ms_run[1][ms_run] = MzTabDouble();
    }

    ConsensusFeature::HandleSetType fs = c.getFeatures();
    for (auto fit = fs.begin(); fit != fs.end(); ++fit)
    {
      UInt study_variable{1};
      const int index = fit->getMapIndex();
      const ConsensusMap::ColumnHeader& ch = cm_column_headers.at(index);

      UInt label = ch.getLabelAsUInt(experiment_type);
      // convert from column index to study variable index
      auto pl = make_pair(ch.filename, label);
      study_variable = path_label_to_assay.at(pl) + 1; // for now, a study_variable is one assay (both 1-based). And pathLabelToSample mapping reports 0-based.

      //TODO implement aggregation in case we generalize study_variable to include multiple assays.
      row.peptide_abundance_stdev_study_variable[study_variable];
      row.peptide_abundance_std_error_study_variable[study_variable];
      row.peptide_abundance_study_variable[study_variable] = MzTabDouble(fit->getIntensity());

      if (export_subfeatures)
      {
        MzTabOptionalColumnEntry opt_global_mass_to_charge_study_variable;
        opt_global_mass_to_charge_study_variable.first = "opt_global_mass_to_charge_study_variable[" + String(study_variable) + "]";
        opt_global_mass_to_charge_study_variable.second = MzTabString(String(fit->getMZ()));
        row.opt_.push_back(opt_global_mass_to_charge_study_variable);

        MzTabOptionalColumnEntry opt_global_retention_time_study_variable;
        opt_global_retention_time_study_variable.first = "opt_global_retention_time_study_variable[" + String(study_variable) + "]";
        opt_global_retention_time_study_variable.second = MzTabString(String(fit->getRT()));
        row.opt_.push_back(opt_global_retention_time_study_variable);
        }
      }

    const vector<PeptideIdentification>& curr_pep_ids = c.getPeptideIdentifications();
    if (!curr_pep_ids.empty())
    {
      checkSequenceUniqueness_(curr_pep_ids);

      // Overall information for this feature in PEP section
      // Features need to be resolved for this. First is not necessarily the best since ids were resorted by map_index.
      const PeptideIdentification& best_id = curr_pep_ids[0];
      const PeptideHit& best_ph = best_id.getHits()[0];
      const AASequence& aas = best_ph.getSequence();
      row.sequence = MzTabString(aas.toUnmodifiedString());

      // annotate variable modifications (no fixed ones)
      row.modifications = extractModificationList(best_ph, fixed_mods, vector<String>());

      const set<String>& accessions = best_ph.extractProteinAccessionsSet();
      const vector<PeptideEvidence> &peptide_evidences = best_ph.getPeptideEvidences();

      row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
      // select accession of first peptide_evidence as representative ("leading") accession
      row.accession = peptide_evidences.empty() ? MzTabString() : MzTabString(peptide_evidences[0].getProteinAccession());

      // fill opt_ columns based on best ID in the feature
      vector<String> id_keys;
      best_id.getKeys(id_keys);

      for (Size k = 0; k != id_keys.size(); ++k)
      {
        String mztabstyle_key = id_keys[k];
        std::replace(mztabstyle_key.begin(), mztabstyle_key.end(), ' ', '_');

        // find matching entry in opt_ (TODO: speed this up)
        for (Size i = 0; i != row.opt_.size(); ++i)
        {
          MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

          if (opt_entry.first == String("opt_global_") + mztabstyle_key)
          {
            opt_entry.second = MzTabString(best_id.getMetaValue(id_keys[k]).toString());
          }
        }
      }

      // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
      for (Size i = 0; i != row.opt_.size(); ++i)
      {
        MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

        if (opt_entry.first == "opt_global_cv_MS:1000889_peptidoform_sequence")
        {
          opt_entry.second = MzTabString(aas.toString());
        }
      }

      // fill opt_ column of psm
      vector<String> ph_keys;
      best_ph.getKeys(ph_keys);

      for (Size k = 0; k != ph_keys.size(); ++k)
      {
        String mztabstyle_key = ph_keys[k];
        std::replace(mztabstyle_key.begin(), mztabstyle_key.end(), ' ', '_');

        // find matching entry in opt_ (TODO: speed this up)
        for (Size i = 0; i != row.opt_.size(); ++i)
        {
          MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

          if (opt_entry.first == String("opt_global_") + mztabstyle_key)
          {
            opt_entry.second = MzTabString(best_ph.getMetaValue(ph_keys[k]).toString());
          }
        }
      }

      // get msrun indices for each ID and insert best search_engine_score for this run
      // for the best run we also annotate the spectra_ref (since it is not designed to be a list)
      double best_score = best_ph.getScore();
      for (const auto& pep : curr_pep_ids)
      {
        size_t spec_run_index = idrun_2_run_index.at(pep.getIdentifier());
        StringList filenames;
        prot_id[spec_run_index].getPrimaryMSRunPath(filenames);
        size_t msfile_index(0);
        size_t id_merge_index(0);
        //TODO synchronize information from ID structures and quant structures somehow.
        // e.g. this part of the code now parses the ID information.
        // This is done because in IsobaricLabelling there is only one ID Run for the different labels
        if (filenames.size() <= 1) //either none or only one file for this run
        {
          msfile_index = map_run_fileidx_2_msfileidx.at({spec_run_index, 0});
        }
        else
        {
          if (pep.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
          {
            id_merge_index = pep.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
            msfile_index = map_run_fileidx_2_msfileidx.at({spec_run_index, id_merge_index});
          }
          else
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "Multiple files in a run, but no id_merge_index in PeptideIdentification found.");
          }
        }

        double curr_score = pep.getHits()[0].getScore();
        auto sit = row.search_engine_score_ms_run[1].find(msfile_index);
        if (sit == row.search_engine_score_ms_run[1].end())
        {
          String ref = "";
          if (pep.metaValueExists("spectrum_reference"))
          {
            ref = pep.getMetaValue("spectrum_reference");
          }
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "PSM " + ref + " does not map to an MS file registered in the quantitative metadata. "
                                              "Check your merging and filtering steps and/or report the issue, please.");
        }
        sit->second = MzTabDouble(curr_score);

        //TODO assumes same scores & score types
        if ((pep.isHigherScoreBetter() && curr_score >= best_score)
        || (!pep.isHigherScoreBetter() && curr_score <= best_score))
        {
          best_score = curr_score;
          if (pep.metaValueExists("spectrum_reference"))
          {
            row.spectra_ref.setSpecRef(pep.getSpectrumReference());
            row.spectra_ref.setMSFile(msfile_index);
          }
        }
      }
      row.best_search_engine_score[1] = MzTabDouble(best_score);
    }

    remapTargetDecoyPSMAndPeptideSection_(row.opt_);
    return row;
  }

  std::optional<MzTabPSMSectionRow> MzTab::PSMSectionRowFromPeptideID_(
     const PeptideIdentification& pid,
     const vector<const ProteinIdentification*>& prot_ids,
     map<String, size_t>& idrun_2_run_index,
     map<pair<size_t,size_t>,size_t>& map_run_fileidx_2_msfileidx,
     map<Size, vector<pair<String, String>>>& run_to_search_engines,
     const Size current_psm_idx,
     const Size psm_id,
     const MzTabString& db,
     const MzTabString& db_version,
     const bool export_empty_pep_ids,
     const bool export_all_psms)
  {
    // skip empty peptide identification objects, if they are not wanted
    if (pid.getHits().empty() && !export_empty_pep_ids)
    {
      return std::nullopt;
    }

    /////// Information that doesn't require a peptide hit ///////
    MzTabPSMSectionRow row;
    row.PSM_ID = MzTabInteger(int(psm_id));
    row.database = db;
    row.database_version = db_version;
    
    vector<MzTabDouble> rts_vector;
    rts_vector.emplace_back(pid.getRT());

    MzTabDoubleList rts;
    rts.set(rts_vector);
    row.retention_time = rts;

    row.exp_mass_to_charge = MzTabDouble(pid.getMZ());

    // meta data on peptide identifications
    vector<String> pid_keys;
    pid.getKeys(pid_keys);

    set<String> pid_key_set(pid_keys.begin(), pid_keys.end());
    addMetaInfoToOptionalColumns(pid_key_set, row.opt_, String("global"), pid);

    // link to spectrum in MS run
    String spectrum_nativeID = pid.getSpectrumReference();
    size_t run_index = idrun_2_run_index.at(pid.getIdentifier());
    StringList filenames;
    prot_ids[run_index]->getPrimaryMSRunPath(filenames);

    StringList localization_mods;
    if (prot_ids[run_index]->getSearchParameters().metaValueExists(Constants::UserParam::LOCALIZED_MODIFICATIONS_USERPARAM))
    {
      localization_mods = prot_ids[run_index]->getSearchParameters().getMetaValue(Constants::UserParam::LOCALIZED_MODIFICATIONS_USERPARAM);
    }

    size_t msfile_index(0);
    if (filenames.size() <= 1) //either none or only one file for this run
    {
      msfile_index = map_run_fileidx_2_msfileidx[{run_index, 0}];
    }
    else
    {
      if (pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
      {
        msfile_index = map_run_fileidx_2_msfileidx[{run_index, pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)}];
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Multiple files in a run, but no id_merge_index in PeptideIdentification found.");
      }
    }

    MzTabSpectraRef spec_ref;
    row.spectra_ref.setMSFile(msfile_index);
    if (spectrum_nativeID.empty())
    {
      OPENMS_LOG_WARN << "spectrum_reference not set in ID with precursor (RT, m/z) " << pid.getRT() << ", " << pid.getMZ() << endl;
    }
    else
    {
      row.spectra_ref.setSpecRef(spectrum_nativeID);
    }

    const vector<PeptideHit>& phs = pid.getHits();

    // add the row and continue to next PepID, if the current one was an empty one
    if (phs.empty())
    { 
      return row;
    }

    /////// Information that does require a peptide hit ///////
    PeptideHit current_ph;
    if (!export_all_psms)
    {
      // only consider best peptide hit for export
      vector<PeptideIdentification> dummy;
      dummy.push_back(pid);
      IDFilter::getBestHit<PeptideIdentification>(dummy, false, current_ph); // TODO: add getBestHit for PeptideHits so no copying to dummy is needed
    }
    else
    {
      current_ph = phs.at(current_psm_idx);
    }

    const AASequence& aas = current_ph.getSequence();
    row.sequence = MzTabString(aas.toUnmodifiedString());

    // extract all modifications in the current sequence for reporting.
    // In contrast to peptide and protein section where fixed modifications are not reported we now report all modifications.
    // If localization mods are specified we add localization scores
    row.modifications = extractModificationList(current_ph, vector<String>(), localization_mods);
    
    MzTabParameterList search_engines;

    //TODO support columns for multiple search engines/scores    
    pair<String, String> name_version = *run_to_search_engines[run_index].begin();
    search_engines.fromCellString("[,," + name_version.first + "," + name_version.second + "]");
    row.search_engine = search_engines;

    row.search_engine_score[1] = MzTabDouble(current_ph.getScore());
    
    row.charge = MzTabInteger(current_ph.getCharge());
    row.calc_mass_to_charge = current_ph.getCharge() != 0 ? MzTabDouble(aas.getMZ(current_ph.getCharge())) : MzTabDouble();

    // add opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
    MzTabOptionalColumnEntry opt_entry;
    opt_entry.first = "opt_global_cv_MS:1000889_peptidoform_sequence";
    opt_entry.second = MzTabString(aas.toString());
    row.opt_.push_back(opt_entry);

    // meta data on PSMs
    vector<String> ph_keys;
    current_ph.getKeys(ph_keys);

    set<String> ph_key_set(ph_keys.begin(), ph_keys.end());
    addMetaInfoToOptionalColumns(ph_key_set, row.opt_, String("global"), current_ph);

    // TODO Think about if the uniqueness can be determined by # of peptide evidences
    //  b/c this would only differ when evidences come from different DBs
    // TODO This also does not consider protein groups but this might be fine here
    const set<String>& accessions = current_ph.extractProteinAccessionsSet();
    row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);

    // create row for every PeptideEvidence entry (mapping to a protein)
    const vector<PeptideEvidence>& peptide_evidences = current_ph.getPeptideEvidences();

    // add peptide evidences to Rows
    row.addPepEvidenceToRows(peptide_evidences);
  
    // remap target/decoy column
    remapTargetDecoyPSMAndPeptideSection_(row.opt_);

    return row;
  }

  size_t MzTab::getQuantStudyVariables_(const ProteinIdentification& pid)
  {
    size_t quant_study_variables(0);
    for (auto & p : pid.getIndistinguishableProteins())
    {
      if (p.getFloatDataArrays().empty()
        || p.getFloatDataArrays()[0].getName() != "abundances")
      {
        quant_study_variables = 0;
        break;
      }
      quant_study_variables = p.getFloatDataArrays()[0].size();
    }
    return quant_study_variables; 
  }

  MzTabParameter MzTab::getProteinScoreType_(const ProteinIdentification& prot_id)
  {
    MzTabParameter protein_score_type;

    bool has_inference_data = prot_id.hasInferenceData();

    if (has_inference_data)
    {
      protein_score_type.fromCellString("[,," + prot_id.getInferenceEngine() + " " + prot_id.getScoreType() 
         + ",]"); // TODO: check if we need one for every run (should not be redundant!)
    }
    else
    {
      // if there was no inference all proteins just come from PeptideIndexer which kind of does a one-peptide rule
      // TODO actually: where are scores coming from in this case. Better to just not write any proteins IMHO
      protein_score_type.fromCellString("[,,one-peptide-rule " + prot_id.getScoreType() + ",]");
    }
    return protein_score_type;    
  }

  void MzTab::mapIDRunFileIndex2MSFileIndex_(
    const vector<const ProteinIdentification*>& prot_ids, 
    const map<String, size_t>& msfilename_2_msrunindex,
    bool skip_first_run, 
    std::map<std::pair<size_t,size_t>,size_t>& map_run_fileidx_2_msfileidx)
  {
    //TODO the following assumes that every file occurs max. once in all runs
    size_t run_index(0);
    for (const auto& run : prot_ids)
    {
      // First entry might be the inference result without (single) associated ms_run. We skip it.
      if (skip_first_run && run_index == 0)
      {
        ++run_index;
        continue;
      }

      StringList files;
      run->getPrimaryMSRunPath(files);
      if (!files.empty())
      {
        size_t file_index(0);
        for (const String& file : files)
        {
          map_run_fileidx_2_msfileidx[{run_index,file_index}] = msfilename_2_msrunindex.at(file);
          file_index++;
        }
      }
      else
      {
        OPENMS_LOG_WARN << "No MS file associated (primary MS run path)." << endl;
        map_run_fileidx_2_msfileidx[{run_index,0}] = run_index;
      }
      ++run_index;
    }
  }

  void MzTab::mapBetweenRunAndSearchEngines_(
    const vector<const ProteinIdentification*>& prot_ids,
    const vector<const PeptideIdentification*>& pep_ids,
    bool skip_first_run,
    map<tuple<String, String, String>, set<Size>>& search_engine_to_runs,
    map<Size, vector<pair<String, String>>>& run_to_search_engines,
    map<Size, vector<vector<pair<String, String>>>>& run_to_search_engine_settings,
    map<String, vector<pair<String, String>>>& search_engine_to_settings)
  {    
    size_t run_index(0);
    for (auto it = prot_ids.cbegin(); it != prot_ids.cend(); ++it)
    {
      // First entry might be the inference result without (single) associated ms_run. We skip it.
      if (skip_first_run && it == prot_ids.cbegin())
      {
        run_index++;
        continue;
      }

      const String &search_engine_name = prot_ids[run_index]->getSearchEngine();
      const String &search_engine_version = prot_ids[run_index]->getSearchEngineVersion();

      String search_engine_score_type = "unknown_score";

      // this is very inefficient.. but almost the only way
      for (const auto& pep : pep_ids)
      {
        if (pep->getIdentifier() == (*it)->getIdentifier())
        {
          search_engine_score_type = pep->getScoreType();
          break;
        }
      }

      search_engine_to_runs[make_tuple(search_engine_name, search_engine_version, search_engine_score_type)].insert(run_index);

      // store main search engine as first entry in run_to_search_engines
      run_to_search_engines[run_index].push_back(make_pair(search_engine_name, search_engine_version));

      vector<String> mvkeys;
      const ProteinIdentification::SearchParameters& sp2 = prot_ids[run_index]->getSearchParameters();
      sp2.getKeys(mvkeys);

      for (const String & mvkey : mvkeys)
      {
        // this is how search engines get overwritten by PercolatorAdapter or ConsensusID
        if (mvkey.hasPrefix("SE:"))
        {
          String se_name = mvkey.substr(3);
          String se_ver = sp2.getMetaValue(mvkey);
          run_to_search_engines[run_index].emplace_back(se_name, se_ver);
          // TODO conserve score_type of underlying search engines (currently always "")
          // TODO for now we only save the MAIN search engine in the SE_to_runs, to only have
          //  those in the meta_data later. -> No discrepancy with the rows (where we also only use
          //  the main search engine
          //search_engine_to_runs[make_tuple(se_name, se_ver, "")].insert(run_index);
        }
      }

      for (const auto& run_se_ver : run_to_search_engines[run_index])
      {
        const auto& se_setting_pairs = prot_ids[run_index]->getSearchEngineSettingsAsPairs(run_se_ver.first);
        // we currently only record the first occurring settings for each search engine.
        search_engine_to_settings.emplace(run_se_ver.first, se_setting_pairs);
        auto it_inserted = run_to_search_engine_settings.emplace(run_index, vector<vector<pair<String, String>>>{se_setting_pairs});
        if (!it_inserted.second)
        {
          it_inserted.first->second.emplace_back(se_setting_pairs);
        }
      }

      ++run_index;
    }
  }

  // static
  MzTabString MzTab::getModificationIdentifier_(const ResidueModification& r)
  {
    String unimod = r.getUniModAccession();
    unimod.toUpper();
    if (!unimod.empty())
    {
      return MzTabString(unimod);
    }
    else
    {
      MzTabString non_unimod_accession = MzTabString("CHEMMOD:" + String(r.getDiffMonoMass()));
      return non_unimod_accession;
    }
  }

  MzTabProteinSectionRow MzTab::proteinSectionRowFromProteinHit_(
    const ProteinHit& hit,
    const MzTabString& db,
    const MzTabString& db_version,
    const set<String>& protein_hit_user_value_keys)
  {
    MzTabProteinSectionRow protein_row;
    protein_row.accession = MzTabString(hit.getAccession());
    protein_row.description = MzTabString(hit.getDescription());
 // protein_row.taxid = hit.getTaxonomyID(); // TODO maybe add as meta value to protein hit NEWT taxonomy for the species.
 // MzTabString species = hit.getSpecies(); // Human readable name of the species
    protein_row.database = db; // Name of the protein database.
    protein_row.database_version = db_version; // String Version of the protein database.
    protein_row.best_search_engine_score[1] = MzTabDouble(hit.getScore());
 // MzTabParameterList search_engine; // Search engine(s) identifying the protein.
 // std::map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
 // std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
 // MzTabInteger reliability;
 // std::map<Size, MzTabInteger> num_psms_ms_run;
 // std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
 // std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
    MzTabModificationList modifications; // Modifications identified in the protein.
    const std::set<pair<Size, ResidueModification>>& leader_mods = hit.getModifications();
    for (auto const & m : leader_mods)
    {
      MzTabModification mztab_mod;
      mztab_mod.setModificationIdentifier(MzTab::getModificationIdentifier_(m.second));
      vector<std::pair<Size, MzTabParameter> > pos;
      pos.emplace_back(make_pair(m.first, MzTabParameter())); // position, parameter pair (e.g. FLR)
      mztab_mod.setPositionsAndParameters(pos);
    }
    protein_row.modifications = modifications;

 // MzTabString uri; // Location of the protein’s source entry.
 // MzTabStringList go_terms; // List of GO terms for the protein.
    double coverage = hit.getCoverage() / 100.0; // convert percent to fraction
    protein_row.coverage = coverage >= 0 ? MzTabDouble(coverage) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
 // std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”

    // create and fill opt_ columns for protein hit user values
    addMetaInfoToOptionalColumns(protein_hit_user_value_keys, protein_row.opt_, String("global"), hit);

    // optional column for protein groups
    MzTabOptionalColumnEntry opt_column_entry;
    opt_column_entry.first = "opt_global_result_type";
    opt_column_entry.second = MzTabString("protein_details");
    protein_row.opt_.push_back(opt_column_entry);

    remapTargetDecoyProteinSection_(protein_row.opt_);

    return protein_row;
  }

  MzTabProteinSectionRow MzTab::nextProteinSectionRowFromProteinGroup_(
    const ProteinIdentification::ProteinGroup& group,
    const MzTabString& db,
    const MzTabString& db_version)
  {
    MzTabProteinSectionRow protein_row;
    protein_row.database = db; // Name of the protein database.
    protein_row.database_version = db_version; // String Version of the protein database.

    MzTabStringList ambiguity_members;
    ambiguity_members.setSeparator(',');
    vector<MzTabString> entries;
    for (Size j = 0; j != group.accessions.size() ; ++j)
    {
      // set accession and description to first element of group
      if (j == 0)
      {
        protein_row.accession = MzTabString(group.accessions[j]);
        // protein_row.description  // TODO: how to set description? information not contained in group
      }
      entries.emplace_back(group.accessions[j]);
    }
    ambiguity_members.set(entries);
    protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
    protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);

    protein_row.coverage = MzTabDouble();

    MzTabOptionalColumnEntry opt_column_entry;
    opt_column_entry.first = "opt_global_result_type";
    opt_column_entry.second = MzTabString("general_protein_group");
    protein_row.opt_.push_back(opt_column_entry);

    remapTargetDecoyProteinSection_(protein_row.opt_);

    return protein_row;
  }

  MzTabProteinSectionRow MzTab::nextProteinSectionRowFromIndistinguishableGroup_(
    const std::vector<ProteinHit>& protein_hits,
    const ProteinIdentification::ProteinGroup& group,
    const size_t g,
    const map<Size, set<Size>>& ind2prot,
    const MzTabString& db,
    const MzTabString& db_version)
  {
    MzTabProteinSectionRow protein_row;

    // get references (indices) into proteins vector
    const set<Size> & protein_hits_idx = ind2prot.at(g);

    // determine group leader
    const ProteinHit& leader_protein = protein_hits[*protein_hits_idx.begin()];

    protein_row.database = db; // Name of the protein database.
    protein_row.database_version = db_version; // String Version of the protein database.

    // column: accession and ambiguity_members
    MzTabStringList ambiguity_members;
    ambiguity_members.setSeparator(',');
    vector<MzTabString> entries;

    // set accession and description to first element of group
    protein_row.accession = MzTabString(leader_protein.getAccession());

    // TODO: check with standard if it is important to also place leader at first position
    //       (because order in set and vector may differ)
    for (Size j = 0; j != group.accessions.size() ; ++j)
    {
      entries.emplace_back(group.accessions[j]);
    }
    ambiguity_members.set(entries);
    protein_row.ambiguity_members = ambiguity_members; // set of indistinguishable proteins

    // annotate if group contains only one or multiple proteins
    MzTabOptionalColumnEntry opt_column_entry;
    opt_column_entry.first = "opt_global_result_type";

    // TODO: we could count the number of targets or set it to target if at least one target is inside the group

    const String col_name = entries.size() == 1 ? "single_protein" :  "indistinguishable_protein_group";
    opt_column_entry.second = MzTabString(col_name);
    protein_row.opt_.push_back(opt_column_entry);

    // column: coverage
    // calculate mean coverage from individual protein coverages
    double coverage{0};
    for (const Size & prot_idx : protein_hits_idx)
    {
      coverage += (1.0 / (double)protein_hits_idx.size()) * 0.01 * protein_hits[prot_idx].getCoverage();
    }
    if (coverage >= 0) { protein_row.coverage = MzTabDouble(coverage); }

    // Store quantitative value attached to abundances in study variables
    if (!group.getFloatDataArrays().empty()
      && group.getFloatDataArrays()[0].getName() == "abundances")
    {
      const ProteinIdentification::ProteinGroup::FloatDataArray & fa = group.getFloatDataArrays()[0];
      Size s(1);
      for (float f : fa)
      {
        protein_row.protein_abundance_assay[s] = MzTabDouble(f); // assay has same information as SV (without design)
        protein_row.protein_abundance_study_variable[s] = MzTabDouble(f);
        protein_row.protein_abundance_stdev_study_variable[s] = MzTabDouble();
        protein_row.protein_abundance_std_error_study_variable[s] = MzTabDouble();
        ++s;
      }
    }

    // add protein description of first (leader) protein
    protein_row.description = MzTabString(leader_protein.getDescription());
    protein_row.taxid = (leader_protein.metaValueExists("TaxID")) ?
      MzTabInteger(static_cast<int>(leader_protein.getMetaValue("TaxID"))) :
      MzTabInteger();

    protein_row.species = (leader_protein.metaValueExists("Species")) ?
      MzTabString(leader_protein.getMetaValue("Species")) :
      MzTabString();

    protein_row.uri = (leader_protein.metaValueExists("URI")) ?
      MzTabString(leader_protein.getMetaValue("URI")) :
      MzTabString();

    if (leader_protein.metaValueExists("GO"))
    {
      StringList sl = leader_protein.getMetaValue("GO");
      String s{};
      s.concatenate(sl.begin(), sl.end(), ",");
      protein_row.go_terms.fromCellString(s);
    }

    protein_row.best_search_engine_score[1] = MzTabDouble(group.probability); // TODO: group probability or search engine score?

    protein_row.reliability = MzTabInteger();

    MzTabParameterList search_engine; // Search engine(s) identifying the protein.
    protein_row.search_engine = search_engine;

    MzTabModificationList modifications; // Modifications identified in the protein.
    const std::set<pair<Size, ResidueModification>>& leader_mods = leader_protein.getModifications();
    for (auto const & m : leader_mods)
    {
      MzTabModification mztab_mod;
      mztab_mod.setModificationIdentifier(MzTab::getModificationIdentifier_(m.second));
      vector<std::pair<Size, MzTabParameter> > pos;

      // mzTab position is one-based, internal is 0-based so we need to +1
      pos.emplace_back(make_pair(m.first + 1, MzTabParameter())); // position, parameter pair (e.g. FLR)
      mztab_mod.setPositionsAndParameters(pos);
      vector<MzTabModification> mztab_mods(1, mztab_mod);
      modifications.set(mztab_mods);
    }
    protein_row.modifications = modifications;

    if (leader_protein.metaValueExists("num_psms_ms_run"))
    {
      const IntList& il = leader_protein.getMetaValue("num_psms_ms_run");
      for (Size ili = 0; ili != il.size(); ++ili)
      {
        protein_row.num_psms_ms_run[ili+1] = MzTabInteger(il[ili]);
      }
    }

    if (leader_protein.metaValueExists("num_peptides_distinct_ms_run"))
    {
      const IntList& il = leader_protein.getMetaValue("num_peptides_distinct_ms_run");
      for (Size ili = 0; ili != il.size(); ++ili)
      {
        protein_row.num_peptides_distinct_ms_run[ili+1] = MzTabInteger(il[ili]);
      }
    }

    if (leader_protein.metaValueExists("num_peptides_unique_ms_run"))
    {
      const IntList& il = leader_protein.getMetaValue("num_peptides_unique_ms_run");
      for (Size ili = 0; ili != il.size(); ++ili)
      {
        protein_row.num_peptides_unique_ms_run[ili+1] = MzTabInteger(il[ili]);
      }
    }

/*
TODO:
Not sure how to handle these:
       // std::map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
       // std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
*/

    remapTargetDecoyProteinSection_(protein_row.opt_);

    // Add protein(group) row to MzTab
    return protein_row;
  }

  map<String, Size> MzTab::mapIDRunIdentifier2IDRunIndex_(const vector<const ProteinIdentification*>& prot_ids)
  {
    map<String, Size> idrunid_2_idrunindex;
    size_t current_idrun_index(0);
    for (auto const& pid : prot_ids)
    {
      idrunid_2_idrunindex[pid->getIdentifier()] = current_idrun_index;
      ++current_idrun_index;
    }
    return idrunid_2_idrunindex;
  }

  void MzTab::mapBetweenMSFileNameAndMSRunIndex_(
    const vector<const ProteinIdentification*>& prot_ids, 
    bool skip_first, 
    map<String, size_t>& msfilename_2_msrunindex,
    map<size_t, String>& msrunindex_2_msfilename)
  {
    size_t current_ms_run_index(1);
    bool first = true;
    for (auto const& pid : prot_ids)
    {
      if (skip_first && first)
      {
        first = false;
        continue;
      }

      StringList ms_run_in_data;
      pid->getPrimaryMSRunPath(ms_run_in_data);

      if (!ms_run_in_data.empty())
      {
        // prepend file:// if not there yet
        for (const String& s : ms_run_in_data)
        {
          // use the string without file: prefix for the map
          msrunindex_2_msfilename.emplace(current_ms_run_index, s);
          const auto& msfileidxpair_success = msfilename_2_msrunindex.emplace(s, current_ms_run_index);
          if (msfileidxpair_success.second) // newly inserted
          {
            current_ms_run_index++;
          }
        }
      }
      else
      {
        // next line is a hack. In case we would ever have some idXML where some runs are annotated
        // and others are not. If a run is not annotated use its index as a String key.        
        msrunindex_2_msfilename.emplace(current_ms_run_index, String(current_ms_run_index));
        msfilename_2_msrunindex.emplace(String(current_ms_run_index), current_ms_run_index);
        current_ms_run_index++;
      }
    }
  }

  void MzTab::addMSRunMetaData_(
    const map<size_t, String>& msrunindex_2_msfilename,
    MzTabMetaData& meta_data)
  {
    for (const auto& r2f : msrunindex_2_msfilename)
    {
      MzTabMSRunMetaData ms_run;
      String m = r2f.second;
      if (!m.hasPrefix("file://")) m = String("file://") + m;
      ms_run.location = MzTabString(m);
      meta_data.ms_run[r2f.first] = ms_run;
    }
  }

  map<Size, set<Size>> MzTab::mapGroupsToProteins_(
    const vector<ProteinIdentification::ProteinGroup>& groups, 
    const vector<ProteinHit>& proteins)
  {
    map<Size, set<Size>> group2prot;
    Size idx{0};
    for (const ProteinIdentification::ProteinGroup & p : groups)
    {
      for (const String & a : p.accessions)
      {
        // find protein corresponding to accession stored in group
        auto it = std::find_if(proteins.begin(), proteins.end(), [&a](const ProteinHit & ph)
          {
            return ph.getAccession() == a;
          }
        );
        if (it == proteins.end())
        { 
          continue;
        }
        Size protein_index = std::distance(proteins.begin(), it);
        group2prot[idx].insert(protein_index);
      }
      ++idx;
    }
    return group2prot;
  }

  void MzTab::addSearchMetaData_(
    const vector<const ProteinIdentification*>& prot_ids,
    const map<tuple<String, String, String>, set<Size>>& search_engine_to_runs,
    const map<String, vector<pair<String,String>>>& search_engine_to_settings,
    MzTabMetaData& meta_data,
    bool first_run_inference_only)
  {
    set<String> protein_scoretypes;
    map<pair<String, String>, vector<pair<String,String>>> protein_settings;
    for (const auto& prot_run : prot_ids)
    {
      //TODO this is a little hack to convert back and from
      protein_scoretypes.insert(getProteinScoreType_(*prot_run).toCellString());
      if (prot_run->hasInferenceData())
      {
        String eng = prot_run->getInferenceEngine();
        String ver = prot_run->getInferenceEngineVersion();
        protein_settings.emplace(make_pair(std::move(eng), std::move(ver)),vector<pair<String,String>>{});
        // TODO add settings for inference tools?
      }
      if (first_run_inference_only) break;
    }

    Size cnt(1);
    for (const auto& mztpar : protein_scoretypes)
    {
      MzTabParameter p{};
      //TODO actually we should make a distinction between protein and protein group-level FDRs
      if (mztpar.hasSubstring("q-value"))
      {
        p.fromCellString("[MS,MS:1003117,OpenMS:Target-decoy protein q-value, ]");
      }
      else if (mztpar.hasSubstring("Epifany"))
      {
        p.fromCellString("[MS,MS:1003119,EPIFANY:Protein posterior probability,]");
      }
      else
      {
        p.fromCellString(mztpar);
      }

      meta_data.protein_search_engine_score[cnt] = p;
      cnt++;
    }

    for (const auto& eng_ver_settings : protein_settings)
    {
      MzTabSoftwareMetaData sesoftwaremd;
      MzTabParameter sesoftware;
      if (eng_ver_settings.first.first == "Epifany")
      {
        sesoftware.fromCellString("[MS,MS:1003118,EPIFANY," + eng_ver_settings.first.second + "]");
      }
      else if (eng_ver_settings.first.first == "TOPPProteinInference")
      {
        sesoftware.fromCellString("[MS,MS:1002203,TOPP ProteinInference," + eng_ver_settings.first.second + "]");
      }
      else
      {
        sesoftware.fromCellString("[,," + eng_ver_settings.first.first + "," + eng_ver_settings.first.second + "]");
      }
      sesoftwaremd.software = sesoftware;
      meta_data.software[meta_data.software.size() + 1] = sesoftwaremd;
    }

    //TODO make software a list?? super weird to fill it like this.
    Size sw_idx(meta_data.software.size() + 1); //+1 since we always start with 1 anyway.

    // Print unique search engines as collected in the input map (globally with first setting encountered)
    for (auto const & name_ver_score_to_runs : search_engine_to_runs)
    {
      MzTabSoftwareMetaData sesoftwaremd;
      MzTabParameter sesoftware;
      //TODO decide if we should use the original search engine ontology entries or the TOPPAdapter entries
      if (get<0>(name_ver_score_to_runs.first) == "Comet")
      {
        sesoftware.fromCellString("[MS,MS:1002251,Comet," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      else if (get<0>(name_ver_score_to_runs.first) == "Percolator")
      {
        sesoftware.fromCellString("[MS,MS:1001490,Percolator," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      else if (get<0>(name_ver_score_to_runs.first) == "MS-GF+" || get<0>(name_ver_score_to_runs.first) == "MSGFPlus")
      {
        sesoftware.fromCellString("[MS,MS:1002048,MS-GF+," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      else if (get<0>(name_ver_score_to_runs.first) == "XTandem")
      {
        sesoftware.fromCellString("[MS,MS:1001476,X!Tandem," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      else if (get<0>(name_ver_score_to_runs.first).hasSubstring("ConsensusID"))
      {
        sesoftware.fromCellString("[MS,MS:1002188,TOPP ConsensusID," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      else
      {
        sesoftware.fromCellString("[,," + get<0>(name_ver_score_to_runs.first) + "," + get<1>(name_ver_score_to_runs.first) + "]");
      }
      sesoftwaremd.software = sesoftware;

      Size cnt2(1);
      for (auto const & sesetting : search_engine_to_settings.at(get<0>(name_ver_score_to_runs.first)))
      {
        sesoftwaremd.setting[cnt2] = MzTabString(sesetting.first + ":" + (!sesetting.second.empty() ? sesetting.second : "null"));
        cnt2++;
      }
      meta_data.software[sw_idx] = sesoftwaremd;
      sw_idx++;
    }

    //TODO actually when filling search_engine_to_runs we would need to go through all MetaValues
    // in the PeptideHits to get all scores and list them.
    // But we do not really annotate which meta values are scores, so we would need to do best guesses
    // Therefore we currently only list the main score type and potentially report additional scores as
    // opt_ columns
    Size psm_search_engine_index(1);
    for (auto const & se : search_engine_to_runs) // loop over (unique) search engine names
    {
      //Huge TODO: we need to somehow correctly support peptide-level scores
      MzTabParameter psm_score_type;
      MzTabParameter pep_score_type;
      const tuple<String, String, String>& name_version_score = se.first;

      psm_score_type.fromCellString("[,," + get<0>(name_version_score) + " " + get<2>(name_version_score) + ",]");
      pep_score_type.fromCellString("[MS,MS:1003114,OpenMS:Best PSM Score,]");

      //TODO we also should consider the different q-value calculation methods of Percolator....
      if (get<2>(name_version_score) == "peptide-level q-value")
      {
        if (get<0>(name_version_score) == "Percolator")
        {
          pep_score_type.fromCellString("[MS,MS:1002360,distinct peptide-level FDRScore,Percolator]");
          psm_score_type = pep_score_type; // since we have no way to have two types
        }
        else
        {
          pep_score_type.fromCellString("[MS,MS:1003116,OpenMS:Target-decoy peptide q-value,]");
          psm_score_type = pep_score_type; // since we have no way to have two types
        }
      }
      else if (get<2>(name_version_score).hasSubstring("q-value"))
      {
        if (get<0>(name_version_score) == "Percolator")
        {
          psm_score_type.fromCellString("[MS,MS:1001491,percolator:Q value,]");
        }
        else
        {
          psm_score_type.fromCellString("[MS,MS:1003115,OpenMS:Target-decoy PSM q-value,]");
        }
      }
      else if (get<2>(name_version_score) == "Posterior Error Probability" || get<2>(name_version_score) == "pep")
      {
        if (get<0>(name_version_score).hasSubstring("ConsensusID"))
        {
          const String& name = get<0>(name_version_score);
          String algo = name.suffix('_');
          psm_score_type.fromCellString("[MS,MS:1003113,OpenMS:ConsensusID PEP," + algo + "]");
        }
        else if (get<0>(name_version_score).hasSubstring("Percolator"))
        {
          psm_score_type.fromCellString("[MS,MS:1001493,percolator:PEP,]");
        }
      }
      else if (get<0>(name_version_score) == "Comet")
      {
        if (get<2>(name_version_score) == "expect")
        {
          psm_score_type.fromCellString("[MS,MS:1002257,Comet:expectation value,]");
        }
        //TODO other Comet scores
      }
      else if (get<0>(name_version_score) == "MSGFPlus" || get<0>(name_version_score) == "MS-GF+")
      {
        if (get<2>(name_version_score) == "SpecEValue")
        {
          psm_score_type.fromCellString("[MS,MS:1002052,MS-GF:SpecEValue,]");
        }
        //TODO other MSGF scores
      }
      //TODO all the other dozens of search engines

      meta_data.psm_search_engine_score[psm_search_engine_index] = psm_score_type;
      meta_data.peptide_search_engine_score[psm_search_engine_index] = pep_score_type; // same score type for peptides
      psm_search_engine_index++;
    }
  }

  MzTabParameter MzTab::getMSRunSpectrumIdentifierType_(const vector<const PeptideIdentification*>& peptide_ids)
  {
    MzTabParameter p;
    p.fromCellString("[MS,MS:1001530,mzML unique identifier,]");
    for (const auto& pid : peptide_ids)
    {
      String spec_ref = pid->getMetaValue("spectrum_reference", "");
      // note: don't change order as some may contain the other terms as well. Taken from mzTab specification document
      if (spec_ref.hasSubstring("controllerNumber=")) { p.fromCellString("[MS,MS:1000768,Thermo nativeID format,]"); return p; }
      if (spec_ref.hasSubstring("process=")) { p.fromCellString("[MS,MS:1000769,Waters nativeID format,]"); return p; }
      if (spec_ref.hasSubstring("cycle=")) { p.fromCellString("[MS,MS:1000770,WIFF nativeID format,]"); return p; }
      if (spec_ref.hasSubstring("scan=")) { p.fromCellString("[MS,MS:1000776,scan number only nativeID format,]"); return p; }
      if (spec_ref.hasSubstring("spectrum=")) { p.fromCellString("[MS,MS:1000777,spectrum identifier nativeID format,]"); return p; }
      return p;
    }
    return p;
  }

  MzTab::IDMzTabStream::IDMzTabStream(
    const std::vector<const ProteinIdentification*>& prot_ids,
    const std::vector<const PeptideIdentification*>& peptide_ids,
    const String& filename,
    bool first_run_inference_only,
    bool export_empty_pep_ids,
    bool export_all_psms,
    const String& title):
      prot_ids_(prot_ids),
      peptide_ids_(peptide_ids),
      filename_(filename),
      export_empty_pep_ids_(export_empty_pep_ids),
      export_all_psms_(export_all_psms)
  {
    ////////////////////////////////////////////////
    // create some lookup structures and precalculate some values
    idrunid_2_idrunindex_ = MzTab::mapIDRunIdentifier2IDRunIndex_(prot_ids_);

    bool has_inference_data = prot_ids_.empty() ? false : prot_ids_[0]->hasInferenceData();

    first_run_inference_ = has_inference_data && first_run_inference_only;
    if (first_run_inference_)
    {
      OPENMS_LOG_INFO << "MzTab: Inference data provided. Considering first run only for inference data." << std::endl;
    }

    map<String, size_t> msfilename_2_msrunindex;
    map<size_t, String> msrunindex_2_msfilename;
    MzTab::mapBetweenMSFileNameAndMSRunIndex_(prot_ids_, first_run_inference_, msfilename_2_msrunindex, msrunindex_2_msfilename);

    // MS runs of a peptide identification object is stored in
    // the protein identification object with the same "identifier".
    // Thus, we build a map from psm_idx->run_index (aka index of PeptideHit -> run index)
    MzTab::mapIDRunFileIndex2MSFileIndex_(prot_ids_, msfilename_2_msrunindex, first_run_inference_, map_id_run_fileidx_2_msfileidx_);

    // collect variable and fixed modifications from different runs
    StringList var_mods;
    MzTab::getSearchModifications_(prot_ids_, var_mods, fixed_mods_);

    // Determine search engines used in the different MS runs.
    map<tuple<String, String, String>, set<Size>> search_engine_to_runs;
    map<String, vector<pair<String,String>>> search_engine_to_settings;

    // search engine and version <-> MS runs index
    MzTab::mapBetweenRunAndSearchEngines_(
      prot_ids_,
      peptide_ids_,
      first_run_inference_,
      search_engine_to_runs,
      run_to_search_engines_,
      run_to_search_engines_settings_,
      search_engine_to_settings);

    ///////////////////////////////////////
    // create column names from meta values
    MzTab::getIdentificationMetaValues_(
      prot_ids, 
      peptide_ids_,
      protein_hit_user_value_keys_,
      peptide_id_user_value_keys_,
      peptide_hit_user_value_keys_);

    // determine nativeID format
    MzTabParameter msrun_spectrum_identifier_type = MzTab::getMSRunSpectrumIdentifierType_(peptide_ids);      

    // filter out redundant meta values
    protein_hit_user_value_keys_.erase("Description"); // already used in Description column

    // construct optional column names
    for (const auto& k : protein_hit_user_value_keys_) prt_optional_column_names_.emplace_back("opt_global_" + k);
    for (const auto& k : peptide_id_user_value_keys_) psm_optional_column_names_.emplace_back("opt_global_" + k);
    for (const auto& k : peptide_hit_user_value_keys_) psm_optional_column_names_.emplace_back("opt_global_" + k);
    
    // rename some of them to be compatible with PRIDE
    std::replace(prt_optional_column_names_.begin(), prt_optional_column_names_.end(), String("opt_global_target_decoy"), String("opt_global_cv_PRIDE:0000303_decoy_hit")); // for PRIDE
    prt_optional_column_names_.emplace_back("opt_global_result_type");
    std::replace(psm_optional_column_names_.begin(), psm_optional_column_names_.end(), String("opt_global_target_decoy"), String("opt_global_cv_MS:1002217_decoy_peptide")); // for PRIDE
    psm_optional_column_names_.emplace_back("opt_global_cv_MS:1000889_peptidoform_sequence");
 
    ///////////////////////////////////////////////////////////////////////
    // Export protein/-group quantifications (stored as meta value in protein IDs)
    // In this case, the first run is only for inference, get peptide info from the rest of the runs.


    // Check if abundances are annotated to the ind. protein groups
    // if so, we will output the abundances as in a quantification file
    // TODO: we currently assume groups are only in the first run, if at all
    //  if we add a field to an ProtIDRun to specify to which condition it belongs,
    //  a vector of ProtIDRuns can potentially hold multiple groupings with quants
    quant_study_variables_ = prot_ids_.empty() ? 0 : getQuantStudyVariables_(*prot_ids_[0]);

    // export PSMs of peptide identifications
    MzTab mztab;

    // mandatory meta values
    meta_data_.mz_tab_type = MzTabString("Identification");
    meta_data_.mz_tab_mode = MzTabString("Summary");
    meta_data_.description = MzTabString("OpenMS export from ID data");
    meta_data_.title = MzTabString(title);

    meta_data_.variable_mod = generateMzTabStringFromModifications(var_mods);
    meta_data_.fixed_mod = generateMzTabStringFromModifications(fixed_mods_);

    MzTabSoftwareMetaData sw;
    sw.software.fromCellString("[MS,MS:1000752,TOPP software," + VersionInfo::getVersion() + "]");
    meta_data_.software[std::max<size_t>(1u, meta_data_.software.size()+1)] = sw;

    if (!prot_ids_.empty())
    {
      // add filenames to the MSRuns in the metadata section
      MzTab::addMSRunMetaData_(msrunindex_2_msfilename, meta_data_);

      // add search settings to software meta data
      MzTab::addSearchMetaData_(
          prot_ids_,
          search_engine_to_runs,
          search_engine_to_settings,
          meta_data_,
          first_run_inference_);

      // trim db name for rows (full name already stored in meta data)
      const ProteinIdentification::SearchParameters & sp = prot_ids_[0]->getSearchParameters();
      String db_basename = File::basename(sp.db);
      db_ = MzTabString(FileHandler::stripExtension(db_basename));
      db_version_ = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
    }

    // condense consecutive unique MS runs to get the different MS files
    auto it = std::unique(ms_runs_.begin(), ms_runs_.end());
    ms_runs_.resize(std::distance(ms_runs_.begin(), it));
    // TODO according to the mzTab standard an MS run can or should be multiple files, when they are coming from
    //  a pre-fractionated sample -> this sounds more like our fraction groups ?!

    // set run meta data
    Size run_index{1};
    for (String m : ms_runs_)
    {
      MzTabMSRunMetaData mztab_run_metadata;
      mztab_run_metadata.format.fromCellString("[MS,MS:1000584,mzML file,]");
      mztab_run_metadata.id_format = msrun_spectrum_identifier_type;

      // prepend file:// if not there yet
      if (!m.hasPrefix("file://")) {m = String("file://") + m; }

      mztab_run_metadata.location = MzTabString(m);

      meta_data_.ms_run[run_index] = mztab_run_metadata;
      OPENMS_LOG_DEBUG << "Adding MS run for file: " << m << endl;
      ++run_index;
    }
  }

  const MzTabMetaData& MzTab::IDMzTabStream::getMetaData() const
  { 
    return meta_data_; 
  }

  const vector<String>& MzTab::IDMzTabStream::getProteinOptionalColumnNames() const
  {
    return prt_optional_column_names_;
  }

  const vector<String>& MzTab::IDMzTabStream::getPeptideOptionalColumnNames() const
  {
    return pep_optional_column_names_;
  }

  const vector<String>& MzTab::IDMzTabStream::getPSMOptionalColumnNames() const
  {
    return psm_optional_column_names_;
  }

  bool MzTab::IDMzTabStream::nextPRTRow(MzTabProteinSectionRow& row)
  {
    if (prot_ids_.empty()) return false;

    // simple state machine to write out 1. all proteins, 2. all general groups and 3. all indistinguishable groups
state0:
    // done if all protein information is contained in run zero and we enter run 1
    if (first_run_inference_ && prt_run_id_ > 0) return false; // done
    if (prt_run_id_ >= prot_ids_.size()) return false; // done for the first_run_inference_ == false case

    const ProteinIdentification& pid = *prot_ids_[prt_run_id_];
    const std::vector<ProteinHit>& protein_hits = pid.getHits();

    // We only report quantitative data for indistinguishable groups (which may be composed of single proteins).
    // We skip the more extensive reporting of general groups with complex shared peptide relations.
    const std::vector<ProteinIdentification::ProteinGroup>& protein_groups2 = quant_study_variables_ == 0 ? pid.getProteinGroups() : std::vector<ProteinIdentification::ProteinGroup>();
    const std::vector<ProteinIdentification::ProteinGroup>& indist_groups2 = pid.getIndistinguishableProteins();

    if (prt_hit_id_ == 0 && PRT_STATE_ == 0) 
    { // Processing new protein identification run?
      // Map (indist.)protein groups to their protein hits (by index) in this run.
      ind2prot_ = MzTab::mapGroupsToProteins_(pid.getIndistinguishableProteins(), protein_hits);
      pg2prot_ = MzTab::mapGroupsToProteins_(pid.getProteinGroups(), protein_hits);
    }

    if (PRT_STATE_ == 0) // write protein hits
    {
      if (prt_hit_id_ >= protein_hits.size())
      {
        prt_hit_id_ = 0;
        PRT_STATE_ = 1; // continue with next state (!)
      }
      else
      {
        const ProteinHit& protein = protein_hits[prt_hit_id_];
        auto prt_row = MzTab::proteinSectionRowFromProteinHit_(
          protein,
          db_,
          db_version_,
          protein_hit_user_value_keys_);
        ++prt_hit_id_;
        std::swap(row, prt_row);
        return true;
      }
    } 

    if (PRT_STATE_ == 1) // write general groups
    {
      if (prt_group_id_ >= protein_groups2.size())
      {
         prt_group_id_ = 0;
         PRT_STATE_ = 2;
      }
      else
      {
        const ProteinIdentification::ProteinGroup& group = protein_groups2[prt_group_id_];
        auto prt_row = MzTab::nextProteinSectionRowFromProteinGroup_(
          group,
          db_,
          db_version_);
        ++prt_group_id_;
        std::swap(row, prt_row);
        return true;
      }      
    }

    // PRT_STATE_ == 2
    if (prt_indistgroup_id_ >= indist_groups2.size()) 
    {
      prt_indistgroup_id_ = 0;
      prt_hit_id_ = 0;
      PRT_STATE_ = 0;
      ++prt_run_id_; // next protein run
      goto state0;
    }
    else
    {
      const ProteinIdentification::ProteinGroup& group = indist_groups2[prt_indistgroup_id_];
      auto prt_row = MzTab::nextProteinSectionRowFromIndistinguishableGroup_(
        protein_hits,
        group,
        prt_indistgroup_id_,
        ind2prot_,
        db_,
        db_version_);
      ++prt_indistgroup_id_;

      std::swap(row, prt_row);
      return true;
    }
    return false; // should not be reached
  }

  bool MzTab::IDMzTabStream::nextPEPRow(MzTabPeptideSectionRow&)
  {
    return false; // no linked feature information 
  }

  bool MzTab::IDMzTabStream::nextPSMRow(MzTabPSMSectionRow& row)
  {
    if (pep_id_ >= peptide_ids_.size()) 
    {
      return false;
    }
    const PeptideIdentification* pid = peptide_ids_[pep_id_];

    auto psm_row = MzTab::PSMSectionRowFromPeptideID_(
        *pid,
        prot_ids_,
        idrunid_2_idrunindex_,
        map_id_run_fileidx_2_msfileidx_,
        run_to_search_engines_,
        current_psm_idx_,
        psm_id_,
        db_,
        db_version_,
        export_empty_pep_ids_,
        export_all_psms_);

    if (!export_all_psms_ || current_psm_idx_ == pid->getHits().size()-1)
    {
      ++pep_id_;
      current_psm_idx_ = 0;
    }
    else
    {
      ++current_psm_idx_;
    }
    ++psm_id_; //global psm id "counter"

    if (psm_row) // valid row?
    {
      row = *psm_row;
    }
    else
    {
      row = MzTabPSMSectionRow();
    }
    return true;
  }

  MzTab MzTab::exportIdentificationsToMzTab(
    const vector<ProteinIdentification>& prot_ids,
    const vector<PeptideIdentification>& peptide_ids,
    const String& filename,
    bool first_run_inference_only,
    bool export_empty_pep_ids,
    bool export_all_psms,
    const String& title)
  {
    vector<const PeptideIdentification*> pep_ids_ptr;
    pep_ids_ptr.reserve(peptide_ids.size());
    for (const PeptideIdentification& pi : peptide_ids) { pep_ids_ptr.push_back(&pi); }

    vector<const ProteinIdentification*> prot_ids_ptr;
    prot_ids_ptr.reserve(prot_ids.size());
    for (const ProteinIdentification& pi : prot_ids) { prot_ids_ptr.push_back(&pi); }

    IDMzTabStream s(prot_ids_ptr, pep_ids_ptr, filename, first_run_inference_only, export_empty_pep_ids, export_all_psms, title);

    MzTab m;
    m.setMetaData(s.getMetaData());

    MzTabProteinSectionRow prot_row;
    while (s.nextPRTRow(prot_row))
    {
      m.getProteinSectionRows().emplace_back(std::move(prot_row));
    }
    
    MzTabPSMSectionRow psm_row;
    while (s.nextPSMRow(psm_row))
    {
      m.getPSMSectionRows().emplace_back(std::move(psm_row));
    }

    return m;
  }

  MzTabModificationList MzTab::extractModificationList(const PeptideHit& pep_hit, const vector<String>& fixed_mods, const vector<String>& localization_mods)
  {
    const AASequence& aas = pep_hit.getSequence();
    MzTabModificationList mod_list;
    vector<MzTabModification> mods;

    bool has_loc_mods = !localization_mods.empty();
    MzTabParameter localization_score;
    if (has_loc_mods && pep_hit.metaValueExists("Luciphor_global_flr"))
    {
      localization_score.fromCellString("[MS,MS:1002380,false localization rate," + String(pep_hit.getMetaValue("Luciphor_global_flr"))+"]");
    }

    if (aas.isModified())
    {
      if (aas.hasNTerminalModification())
      {
        MzTabModification mod;
        const ResidueModification& res_mod = *(aas.getNTerminalModification());

        bool is_fixed = std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) != fixed_mods.end();
        if (!is_fixed)
        {
          mod.setModificationIdentifier(MzTab::getModificationIdentifier_(res_mod));
          vector<std::pair<Size, MzTabParameter> > pos;
          pos.emplace_back(0, MzTabParameter());
          mod.setPositionsAndParameters(pos);
          mods.push_back(mod);
        }
      }

      for (Size ai = 0; ai != aas.size(); ++ai)
      {
        if (aas[ai].isModified())
        {
          MzTabModification mod;
          const ResidueModification& res_mod = *(aas[ai].getModification());
          bool is_fixed = std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) != fixed_mods.end();
          if (!is_fixed)
          {
            // MzTab standard is to just report Unimod accession.
            vector<std::pair<Size, MzTabParameter> > pos;
            if (has_loc_mods && std::find(localization_mods.begin(), localization_mods.end(), res_mod.getFullId()) != localization_mods.end())
            { // store localization score for this mod
              pos.emplace_back(ai + 1, localization_score);
            }
            else
            {
              pos.emplace_back(ai + 1, MzTabParameter());
            }
            mod.setPositionsAndParameters(pos);
            mod.setModificationIdentifier(MzTab::getModificationIdentifier_(res_mod));
            mods.push_back(mod);
          }
        }
      }

      if (aas.hasCTerminalModification())
      {
        MzTabModification mod;
        const ResidueModification& res_mod = *(aas.getCTerminalModification());
        if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
        {
          vector<std::pair<Size, MzTabParameter> > pos;
          pos.emplace_back(aas.size() + 1, MzTabParameter());
          mod.setPositionsAndParameters(pos);
          mod.setModificationIdentifier(MzTab::getModificationIdentifier_(res_mod));
          mods.push_back(mod);
        }
      }
    }
    mod_list.set(mods);
    return mod_list;
  }

  void MzTab::getFeatureMapMetaValues_(const FeatureMap& feature_map,
    set<String>& feature_user_value_keys, 
    set<String>& peptide_identification_user_value_keys, 
    set<String>& peptide_hit_user_value_keys)
  {
    for (Size i = 0; i < feature_map.size(); ++i)
    {
      const Feature& f = feature_map[i];
      vector<String> keys;
      f.getKeys(keys); //TODO: why not just return it?

      feature_user_value_keys.insert(keys.begin(), keys.end());

      const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
      for (PeptideIdentification const & pep_id : pep_ids)
      {
        vector<String> pep_keys;
        pep_id.getKeys(pep_keys);
        peptide_identification_user_value_keys.insert(pep_keys.begin(), pep_keys.end());

        for (PeptideHit const & hit : pep_id.getHits())
        {
          vector<String> ph_keys;
          hit.getKeys(ph_keys);
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }
    // we don't want spectrum reference to show up as meta value (already in dedicated column)
    peptide_hit_user_value_keys.erase("spectrum_reference");
  }

  // local helper to extract meta values with space substituted with '_'
  void extractMetaValuesFromIDs(const vector<PeptideIdentification> & curr_pep_ids, 
    set<String>& peptide_identification_user_value_keys,
    set<String>& peptide_hit_user_value_keys)
    {
      for (auto const & pep_id : curr_pep_ids)
      {      
        vector<String> pep_keys;
        pep_id.getKeys(pep_keys);
        // replace whitespaces with underscore
        std::transform(pep_keys.begin(), pep_keys.end(), pep_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
        peptide_identification_user_value_keys.insert(pep_keys.begin(), pep_keys.end());

        for (auto const & hit : pep_id.getHits())
        {
          vector<String> ph_keys;
          hit.getKeys(ph_keys);
          // replace whitespaces with underscore
          std::transform(ph_keys.begin(), ph_keys.end(), ph_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

  // local helper to extract meta values with space substituted with '_'
  void extractMetaValuesFromIDPointers(const vector<const PeptideIdentification*> & curr_pep_ids, 
    set<String>& peptide_identification_user_value_keys,
    set<String>& peptide_hit_user_value_keys)
    {
      for (auto const * pep_id : curr_pep_ids)
      {      
        vector<String> pep_keys;
        pep_id->getKeys(pep_keys);
        // replace whitespaces with underscore
        std::transform(pep_keys.begin(), pep_keys.end(), pep_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
        peptide_identification_user_value_keys.insert(pep_keys.begin(), pep_keys.end());

        for (auto const & hit : pep_id->getHits())
        {
          vector<String> ph_keys;
          hit.getKeys(ph_keys);
          // replace whitespaces with underscore
          std::transform(ph_keys.begin(), ph_keys.end(), ph_keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
          peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
        }
      }
    }

  // extract *all* meta values stored at consensus feature, peptide id and peptide hit level
  void MzTab::getConsensusMapMetaValues_(const ConsensusMap& consensus_map,
    set<String>& consensus_feature_user_value_keys,
    set<String>& peptide_identification_user_value_keys,
    set<String>& peptide_hit_user_value_keys)
  {
    // extract meta values from unassigned peptide identifications
    const vector<PeptideIdentification> & curr_pep_ids = consensus_map.getUnassignedPeptideIdentifications();
    extractMetaValuesFromIDs(curr_pep_ids, peptide_identification_user_value_keys, peptide_hit_user_value_keys);

    for (ConsensusFeature const & c : consensus_map)
    {
      vector<String> keys;
      c.getKeys(keys);
      // replace whitespaces with underscore
      std::transform(keys.begin(), keys.end(), keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });

      consensus_feature_user_value_keys.insert(keys.begin(), keys.end());

      // extract meta values from assigned peptide identifications
      const vector<PeptideIdentification> & curr_pep_ids = c.getPeptideIdentifications();
      extractMetaValuesFromIDs(curr_pep_ids, peptide_identification_user_value_keys, peptide_hit_user_value_keys);
    }

    // we don't want spectrum reference to show up as meta value (already in dedicated column)
    peptide_identification_user_value_keys.erase("spectrum_reference");
  }

  void MzTab::getIdentificationMetaValues_(
    const std::vector<const ProteinIdentification*>& prot_ids, 
    std::vector<const PeptideIdentification*>& peptide_ids_,
    std::set<String>& protein_hit_user_value_keys,
    std::set<String>& peptide_id_user_value_keys,
    std::set<String>& peptide_hit_user_value_keys)
  {
    for (auto const & pid : prot_ids)
    {
      for (auto const & hit : pid->getHits())
      {
        vector<String> keys;
        hit.getKeys(keys);
        // replace whitespaces with underscore
        std::transform(keys.begin(), keys.end(), keys.begin(), [&](String& s) { return s.substitute(' ', '_'); });
        protein_hit_user_value_keys.insert(keys.begin(), keys.end());
      }
    }

    extractMetaValuesFromIDPointers(peptide_ids_, peptide_id_user_value_keys, peptide_hit_user_value_keys);
  }

  void MzTab::getSearchModifications_(const vector<const ProteinIdentification*>& prot_ids, StringList& var_mods, StringList& fixed_mods)
  {
    for (auto const & pid : prot_ids)
    {
      const ProteinIdentification::SearchParameters & sp = pid->getSearchParameters();
      var_mods.insert(std::end(var_mods), std::begin(sp.variable_modifications), std::end(sp.variable_modifications));
      fixed_mods.insert(std::end(fixed_mods), std::begin(sp.fixed_modifications), std::end(sp.fixed_modifications));
    }

    // make mods unique
    std::sort(var_mods.begin(), var_mods.end());
    auto v_it = std::unique(var_mods.begin(), var_mods.end());
    var_mods.resize(std::distance(var_mods.begin(), v_it));
    std::sort(fixed_mods.begin(), fixed_mods.end());
    auto f_it = std::unique(fixed_mods.begin(), fixed_mods.end());
    fixed_mods.resize(std::distance(fixed_mods.begin(), f_it));
  }


  MzTab::CMMzTabStream::CMMzTabStream(
    const ConsensusMap& consensus_map,
    const String& filename,
    const bool first_run_inference_only,
    const bool export_unidentified_features,
    const bool export_unassigned_ids,
    const bool export_subfeatures,
    const bool export_empty_pep_ids,
    const bool export_all_psms,
    const String& title) 
  :
    consensus_map_(consensus_map),
    filename_(filename), 
    export_unidentified_features_(export_unidentified_features),
    export_subfeatures_(export_subfeatures),
    export_empty_pep_ids_(export_empty_pep_ids),
    export_all_psms_(export_all_psms)
  {
    // fill ID datastructure without copying
    const vector<ProteinIdentification>& prot_id = consensus_map.getProteinIdentifications();
    prot_ids_.reserve(prot_id.size());
    for (const auto & i : prot_id)
    {
      prot_ids_.push_back(&i);
    }

    // pre-reserve to prevent reallocations
    size_t new_size = peptide_ids_.size();

    for (const auto& elem : consensus_map) new_size += elem.getPeptideIdentifications().size();

    if (export_unassigned_ids)
    {
      new_size += consensus_map.getUnassignedPeptideIdentifications().size();
    }

    peptide_ids_.reserve(new_size);

    // extract mapped IDs
    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      const ConsensusFeature& c = consensus_map[i];
      const vector<PeptideIdentification>& p = c.getPeptideIdentifications();
      for (const PeptideIdentification& pi : p) { peptide_ids_.push_back(&pi); }
    }

    // also export PSMs of unassigned peptide identifications
    if (export_unassigned_ids)
    {
      const vector<PeptideIdentification>& up = consensus_map.getUnassignedPeptideIdentifications();
      for (const PeptideIdentification& pi : up) { peptide_ids_.push_back(&pi); }
    }

    ////////////////////////////////////////////////
    // create some lookup structures and precalculate some values
    idrunid_2_idrunindex_ = MzTab::mapIDRunIdentifier2IDRunIndex_(prot_ids_);

    bool has_inference_data = !prot_ids_.empty() && prot_ids_[0]->hasInferenceData();

    first_run_inference_ = has_inference_data && first_run_inference_only;
    if (first_run_inference_)
    {
      OPENMS_LOG_INFO << "MzTab: Inference data provided. Considering first run only for inference data." << std::endl;
    }

    map<String, size_t> msfilename_2_msrunindex;
    map<size_t, String> msrunindex_2_msfilename;
    MzTab::mapBetweenMSFileNameAndMSRunIndex_(prot_ids_, first_run_inference_, msfilename_2_msrunindex, msrunindex_2_msfilename);

    // MS runs of a peptide identification object is stored in
    // the protein identification object with the same "identifier".
    // Thus, we build a map from psm_idx->run_index (aka index of PeptideHit -> run index)
    MzTab::mapIDRunFileIndex2MSFileIndex_(prot_ids_, msfilename_2_msrunindex, first_run_inference_, map_id_run_fileidx_2_msfileidx_);

    // collect variable and fixed modifications from different runs
    StringList var_mods;
    MzTab::getSearchModifications_(prot_ids_, var_mods, fixed_mods_);

    // determine nativeID format
    MzTabParameter msrun_spectrum_identifier_type = MzTab::getMSRunSpectrumIdentifierType_(peptide_ids_);

    // Determine search engines used in the different MS runs.
    map<tuple<String, String, String>, set<Size>> search_engine_to_runs;
    map<String, vector<pair<String,String>>> search_engine_to_settings;

    // search engine and version <-> MS runs index
    MzTab::mapBetweenRunAndSearchEngines_(
      prot_ids_,
      peptide_ids_,
      first_run_inference_,
      search_engine_to_runs,
      run_to_search_engines_,
      run_to_search_engines_settings_,
      search_engine_to_settings);

    // Pre-analyze data for re-occurring meta values at consensus feature and contained peptide id and hit level.
    // These are stored in optional columns of the PEP section.
    MzTab::getConsensusMapMetaValues_(consensus_map, 
      consensus_feature_user_value_keys_,
      consensus_feature_peptide_identification_user_value_keys_,
      consensus_feature_peptide_hit_user_value_keys_);

    // create column names from meta values
    // feature meta values
    for (const auto& k : consensus_feature_user_value_keys_) 
    {
      pep_optional_column_names_.emplace_back("opt_global_" + k);
    }      

    /* 
      Note: In the PEP section, meta values in peptide identifications are better omitted as they can be easily looked up from the PSM-level are otherwise duplicates.
      For debug purposes one can enable this line:
      for (const auto& k : consensus_feature_peptide_identification_user_value_keys_) pep_optional_column_names_.emplace_back("opt_global_" + k);
    */

    // peptide hit (PSM) meta values
    
    // whitelisted meta values "target_decoy". Expose in the PEP section with special CV term (for PRIDE compatibility):
    if (auto it = consensus_feature_peptide_hit_user_value_keys_.find("target_decoy");
      it != consensus_feature_peptide_hit_user_value_keys_.end())
    {
      pep_optional_column_names_.emplace_back("opt_global_cv_MS:1002217_decoy_peptide");
    }
    // added during export (for PRIDE compatibility):
    pep_optional_column_names_.emplace_back("opt_global_cv_MS:1000889_peptidoform_sequence");

    // PSM optional columns: also from meta values in consensus features
    for (const auto& k : consensus_feature_peptide_hit_user_value_keys_) psm_optional_column_names_.emplace_back("opt_global_" + k);
    std::replace(psm_optional_column_names_.begin(), psm_optional_column_names_.end(), String("opt_global_target_decoy"), String("opt_global_cv_MS:1002217_decoy_peptide")); // for PRIDE
    psm_optional_column_names_.emplace_back("opt_global_cv_MS:1000889_peptidoform_sequence");

    ///////////////////////////////////////////////////////////////////////
    // Export protein/-group quantifications (stored as meta value in protein IDs)
    // In this case, the first run is only for inference, get peptide info from the rest of the runs.


    // Check if abundances are annotated to the ind. protein groups
    // if so, we will output the abundances as in a quantification file
    // TODO: we currently assume groups are only in the first run, if at all
    //  if we add a field to an ProtIDRun to specify to which condition it belongs,
    //  a vector of ProtIDRuns can potentially hold multiple groupings with quants
    quant_study_variables_ = prot_ids_.empty() ? 0 : getQuantStudyVariables_(*prot_ids_[0]);

    // export PSMs of peptide identifications
    MzTab mztab;

    // mandatory meta values
    meta_data_.mz_tab_type = MzTabString("Quantification");
    meta_data_.mz_tab_mode = MzTabString("Summary");
    meta_data_.description = MzTabString("OpenMS export from consensusXML");

    meta_data_.variable_mod = generateMzTabStringFromModifications(var_mods);
    meta_data_.fixed_mod = generateMzTabStringFromModifications(fixed_mods_);

    MzTabSoftwareMetaData sw;
    sw.software.fromCellString("[MS,MS:1000752,TOPP software," + VersionInfo::getVersion() + "]");
    meta_data_.software[std::max<size_t>(1u, meta_data_.software.size()+1)] = sw;

    if (!prot_ids_.empty())
    {
      // add filenames to the MSRuns in the metadata section
      MzTab::addMSRunMetaData_(msrunindex_2_msfilename, meta_data_);

      // add search settings to software meta data
      MzTab::addSearchMetaData_(
          prot_ids_,
          search_engine_to_runs,
          search_engine_to_settings,
          meta_data_,
          first_run_inference_);

      // trim db name for rows (full name already stored in meta data)
      const ProteinIdentification::SearchParameters & sp = prot_ids_[0]->getSearchParameters();
      String db_basename = File::basename(sp.db);
      db_ = MzTabString(FileHandler::stripExtension(db_basename));

      db_version_ = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);

      ////////////////////////////////////////////////////////////////
      // generate protein section
      for (auto it = prot_ids_.cbegin(); it != prot_ids_.cend(); ++it)
      {
        const std::vector<ProteinHit>& protein_hits = (*it)->getHits();

        // TODO: add processing information that this file has been exported from "filename"

        // pre-analyze data for occurring meta values at protein hit level
        // these are used to build optional columns containing the meta values in internal data structures
        
        set<String> protein_hit_user_value_keys_tmp =
          MetaInfoInterfaceUtils::findCommonMetaKeys<vector<ProteinHit>, set<String> >(protein_hits.begin(), protein_hits.end(), 100.0);

        // we do not want descriptions twice
        protein_hit_user_value_keys_tmp.erase("Description");

        protein_hit_user_value_keys_.insert(protein_hit_user_value_keys_tmp.begin(), protein_hit_user_value_keys_tmp.end());
      }
    }

    // column headers may not contain spaces
    set<String> protein_hit_user_value_keys_tmp_2;
    // replace whitespaces with underscore
    std::transform(protein_hit_user_value_keys_.begin(),
                   protein_hit_user_value_keys_.end(),
                   std::inserter(protein_hit_user_value_keys_tmp_2, protein_hit_user_value_keys_tmp_2.begin()),
                   [](String s) { return s.substitute(' ', '_'); });

    std::swap(protein_hit_user_value_keys_, protein_hit_user_value_keys_tmp_2);

    // PRT optional columns
    for (const auto& k : protein_hit_user_value_keys_) prt_optional_column_names_.emplace_back("opt_global_" + k);
    std::replace(prt_optional_column_names_.begin(), prt_optional_column_names_.end(), String("opt_global_target_decoy"), String("opt_global_cv_PRIDE:0000303_decoy_hit")); // for PRIDE
    prt_optional_column_names_.emplace_back("opt_global_result_type");

    // determine number of samples
    ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(consensus_map);

    Size n_assays = ed.getNumberOfSamples();

    // TODO for now every assay is a study variable since we do not aggregate across e.g. replicates.
    n_study_variables_ = n_assays;

    ///////////////////////////////////////////////////////////////////////
    // MetaData section

    meta_data_.title = MzTabString(title);

    MzTabParameter quantification_method;
    const String & experiment_type = consensus_map.getExperimentType();
    if (experiment_type == "label-free")
    {
      quantification_method.fromCellString("[MS,MS:1001834,LC-MS label-free quantitation analysis,]");
    }
    else if (experiment_type == "labeled_MS1")
    {
      quantification_method.fromCellString("[PRIDE,PRIDE_0000316,MS1 based isotope labeling,]");
    }
    else if (experiment_type == "labeled_MS2")
    {
      quantification_method.fromCellString("[PRIDE,PRIDE_0000317,MS2 based isotope labeling,]");
    }

    meta_data_.quantification_method = quantification_method;
    MzTabParameter protein_quantification_unit;
    // TODO: add better term to obo: Would need to be a combination of settings:
    //  "sum/mean/max/..fancy", "top3,all,..", "shared,unique,razor+unique", potentially make clear that it is usually
    //  on group level here (even if singleton group).
    protein_quantification_unit.fromCellString("[,,Abundance,]");

    meta_data_.protein_quantification_unit = protein_quantification_unit;
    MzTabParameter peptide_quantification_unit;
    // TODO: I think we could use MS1 feature area or feature height or spectral count here.
    peptide_quantification_unit.fromCellString("[,,Abundance,]");
    meta_data_.peptide_quantification_unit = peptide_quantification_unit;

    consensus_map.getPrimaryMSRunPath(ms_runs_);

    // condense consecutive unique MS runs to get the different MS files
    auto it = std::unique(ms_runs_.begin(), ms_runs_.end());
    ms_runs_.resize(std::distance(ms_runs_.begin(), it));
    // TODO according to the mzTab standard an MS run can or should be multiple files, when they are coming from
    //  a pre-fractionated sample -> this sounds more like our fraction groups ?!

    // set run meta data
    Size run_index{1};
    for (String m : ms_runs_)
    {
      MzTabMSRunMetaData mztab_run_metadata;
      mztab_run_metadata.format.fromCellString("[MS,MS:1000584,mzML file,]");
      mztab_run_metadata.id_format = msrun_spectrum_identifier_type;

      // prepend file:// if not there yet
      if (!m.hasPrefix("file://")) {m = String("file://") + m; }

      mztab_run_metadata.location = MzTabString(m);

      meta_data_.ms_run[run_index] = mztab_run_metadata;
      OPENMS_LOG_DEBUG << "Adding MS run for file: " << m << endl;
      ++run_index;
    }

    // assay index (and sample index) must be unique numbers 1..n
    // fraction_group + label define the quant. values of an assay (which currently corresponds to our Sample ID)
    path_label_to_assay_ = ed.getPathLabelToSampleMapping(false);

    // assay meta data
    for (auto const & c : consensus_map.getColumnHeaders())
    {
      Size assay_index{1};

      MzTabAssayMetaData assay;
      MzTabParameter quantification_reagent;
      Size label = c.second.getLabelAsUInt(experiment_type);
      auto pl = make_pair(c.second.filename, label);
      assay_index = path_label_to_assay_[pl] + 1; // sample rows are a vector and therefore their IDs zero-based, mzTab assays 1-based

      if (experiment_type == "label-free")
      {
        quantification_reagent.fromCellString("[MS,MS:1002038,unlabeled sample,]");
      }
      else if (experiment_type == "labeled_MS1")
      {
        quantification_reagent.fromCellString("[PRIDE,PRIDE:0000316,MS1 based isotope labeling," + c.second.label + "]");
      }
      else if (experiment_type == "labeled_MS2")
      {       
        quantification_reagent.fromCellString("[PRIDE,PRIDE:0000317,MS2 based isotope labeling," + c.second.label + "]");
      }
      
      // look up run index by filename
      //TODO again, check if we rather want fraction groups instead of individual files.
      auto md_it = find_if(meta_data_.ms_run.begin(), meta_data_.ms_run.end(),
        [&c] (const pair<Size, MzTabMSRunMetaData>& m) {
          return m.second.location.toCellString().hasSuffix(c.second.filename);
        } );
      Size curr_run_index = md_it->first;

      meta_data_.assay[assay_index].quantification_reagent = quantification_reagent;
      meta_data_.assay[assay_index].ms_run_ref.push_back(curr_run_index);

      // study variable meta data
      MzTabString sv_description;
      // TODO how would we represent study variables? = Collection of sample rows that are equal except for replicate
      //  columns?
      meta_data_.study_variable[assay_index].description.fromCellString("no description given");
      IntList al;
      al.push_back(assay_index);
      meta_data_.study_variable[assay_index].assay_refs = al;
    }
  }

  const MzTabMetaData& MzTab::CMMzTabStream::getMetaData() const
  { 
    return meta_data_; 
  }

  const vector<String>& MzTab::CMMzTabStream::getProteinOptionalColumnNames() const
  {
    return prt_optional_column_names_;
  }

  const vector<String>& MzTab::CMMzTabStream::getPeptideOptionalColumnNames() const
  {
    return pep_optional_column_names_;
  }

  const vector<String>& MzTab::CMMzTabStream::getPSMOptionalColumnNames() const
  {
    return psm_optional_column_names_;
  }

  bool MzTab::CMMzTabStream::nextPRTRow(MzTabProteinSectionRow& row)
  {
    if (prot_ids_.empty()) return false;

    // simple state machine to write out 1. all proteins, 2. all general groups and 3. all indistinguishable groups
state0:
    // done if all protein information is contained in run zero and we enter run 1
    if (first_run_inference_ && prt_run_id_ > 0) return false; // done
    if (prt_run_id_ >= prot_ids_.size()) return false; // done for the first_run_inference_ == false case

    const ProteinIdentification& pid = *prot_ids_[prt_run_id_];
    const std::vector<ProteinHit>& protein_hits = pid.getHits();

    // We only report quantitative data for indistinguishable groups (which may be composed of single proteins).
    // We skip the more extensive reporting of general groups with complex shared peptide relations.
    const std::vector<ProteinIdentification::ProteinGroup>& protein_groups2 = quant_study_variables_ == 0 ? pid.getProteinGroups() : std::vector<ProteinIdentification::ProteinGroup>();
    const std::vector<ProteinIdentification::ProteinGroup>& indist_groups2 = pid.getIndistinguishableProteins();

    if (prt_hit_id_ == 0 && PRT_STATE_ == 0) 
    { // Processing new protein identification run?
      // Map (indist.)protein groups to their protein hits (by index) in this run.
      ind2prot_ = MzTab::mapGroupsToProteins_(pid.getIndistinguishableProteins(), protein_hits);
      pg2prot_ = MzTab::mapGroupsToProteins_(pid.getProteinGroups(), protein_hits);
    }

    if (PRT_STATE_ == 0) // write protein hits
    {
      if (prt_hit_id_ >= protein_hits.size())
      {
        prt_hit_id_ = 0;
        PRT_STATE_ = 1; // continue with next state (!)
      }
      else
      {
        const ProteinHit& protein = protein_hits[prt_hit_id_];
        auto prt_row = MzTab::proteinSectionRowFromProteinHit_(
          protein,
          db_,
          db_version_,
          protein_hit_user_value_keys_);
        ++prt_hit_id_;
        std::swap(row, prt_row);
        return true;
      }
    } 

    if (PRT_STATE_ == 1) // write general groups
    {
      if (prt_group_id_ >= protein_groups2.size())
      {
         prt_group_id_ = 0;
         PRT_STATE_ = 2;
      }
      else
      {
        const ProteinIdentification::ProteinGroup& group = protein_groups2[prt_group_id_];
        auto prt_row = MzTab::nextProteinSectionRowFromProteinGroup_(
          group,
          db_,
          db_version_);
        ++prt_group_id_;
        std::swap(row, prt_row);
        return true;
      }      
    }

    // PRT_STATE_ == 2
    if (prt_indistgroup_id_ >= indist_groups2.size()) 
    {
      prt_indistgroup_id_ = 0;
      prt_hit_id_ = 0;
      PRT_STATE_ = 0;
      ++prt_run_id_; // next protein run
      goto state0;
    }
    else
    {
      const ProteinIdentification::ProteinGroup& group = indist_groups2[prt_indistgroup_id_];
      auto prt_row = MzTab::nextProteinSectionRowFromIndistinguishableGroup_(
        protein_hits,
        group,
        prt_indistgroup_id_,
        ind2prot_,
        db_,
        db_version_);
      ++prt_indistgroup_id_;

      std::swap(row, prt_row);
      return true;
    }
    return false; // should not be reached
  }

  bool MzTab::CMMzTabStream::nextPEPRow(MzTabPeptideSectionRow& row)
  {
    if (pep_id_ >= consensus_map_.size()) return false; 

    auto c = std::cref(consensus_map_[pep_id_]);

    auto has_peptide_hits = [&](const ConsensusFeature& c) 
      { 
        for (const auto& pid : c.getPeptideIdentifications())
        {
          if (!pid.getHits().empty()) return true;
        }
        return false; 
      };
    
    // skip unidentified features
    while (!export_unidentified_features_ && !has_peptide_hits(c.get()))
    {
      ++pep_id_;
      if (pep_id_ >= consensus_map_.size()) return false;      
      c = consensus_map_[pep_id_];
    }

    auto pep_row = MzTab::peptideSectionRowFromConsensusFeature_(
     c.get(), 
     consensus_map_, 
     ms_runs_,
     n_study_variables_, 
     consensus_feature_user_value_keys_, 
     consensus_feature_peptide_identification_user_value_keys_, 
     consensus_feature_peptide_hit_user_value_keys_,
     idrunid_2_idrunindex_,
     map_id_run_fileidx_2_msfileidx_,
     path_label_to_assay_,
     fixed_mods_,
     export_subfeatures_);

    ++pep_id_;

    std::swap(row, pep_row);
    return true;
  }
  
  bool MzTab::CMMzTabStream::nextPSMRow(MzTabPSMSectionRow& row)
  {
    if (pep_counter_ >= peptide_ids_.size()) return false;
    const PeptideIdentification* pid = peptide_ids_[pep_counter_];
    auto psm_row = MzTab::PSMSectionRowFromPeptideID_(
        *pid,
        prot_ids_,
        idrunid_2_idrunindex_,
        map_id_run_fileidx_2_msfileidx_,
        run_to_search_engines_,
        current_psm_idx_,
        psm_id_,
        db_,
        db_version_,
        export_empty_pep_ids_,
        export_all_psms_);

    if (!export_all_psms_ || current_psm_idx_ == pid->getHits().size()-1)
    {
      ++pep_counter_;
      current_psm_idx_ = 0;
    }
    else
    {
      ++current_psm_idx_;
    }
    ++psm_id_; //global psm id "counter"

    if (psm_row) // valid row?
    {
      row = *psm_row;
    }
    else
    {
      row = MzTabPSMSectionRow();
    }
    return true;
  }

  MzTab MzTab::exportConsensusMapToMzTab(
    const ConsensusMap& consensus_map,
    const String& filename,
    const bool first_run_inference_only,
    const bool export_unidentified_features,
    const bool export_unassigned_ids,
    const bool export_subfeatures,
    const bool export_empty_pep_ids,
    const bool export_all_psms,
    const String& title)
  {  
    OPENMS_LOG_INFO << "exporting consensus map: \"" << filename << "\" to mzTab: " << std::endl;

    CMMzTabStream s(consensus_map,
      filename,
      first_run_inference_only,
      export_unidentified_features,
      export_unassigned_ids,
      export_subfeatures,
      export_empty_pep_ids,
      export_all_psms,
      title);

    MzTab m;
    m.setMetaData(s.getMetaData());

    MzTabProteinSectionRow prot_row;
    while (s.nextPRTRow(prot_row))
    {
      m.getProteinSectionRows().emplace_back(std::move(prot_row));
    }
    
    MzTabPeptideSectionRow pep_row;
    while (s.nextPEPRow(pep_row))
    {
      m.getPeptideSectionRows().emplace_back(std::move(pep_row));
    }

    MzTabPSMSectionRow psm_row;
    while (s.nextPSMRow(psm_row))
    {
      // TODO better return a State enum instead of relying on some uninitialized
      // parts of a row..
      if (!psm_row.sequence.isNull())
      {
        m.getPSMSectionRows().push_back(psm_row);
      }
    }

    return m;
  }

  void MzTab::checkSequenceUniqueness_(const vector<PeptideIdentification>& curr_pep_ids)
  {
    const auto& refseq = curr_pep_ids[0].getHits()[0].getSequence();
    for (const auto& pep : curr_pep_ids)
    {
      if (pep.getHits()[0].getSequence() != refseq)
      {
        throw OpenMS::Exception::IllegalArgument(
            __FILE__
            , __LINE__
            , __FUNCTION__
            , "Consensus features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!");
      }
    }
  }
}
