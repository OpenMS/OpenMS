// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MzTabExporter MzTabExporter

   @brief This application converts several %OpenMS XML formats (featureXML, consensusXML, and idXML) to mzTab.

  <CENTER>
    <table>
     <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MzTabExporter \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any tool producing one of the input formats </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> External tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

  @experimental This algorithm and underlying format is work in progress and might change.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MzTabExporter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_MzTabExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

namespace OpenMS
{
  class TOPPMzTabExporter :
    public TOPPBase
  {
public:
    TOPPMzTabExporter() :
      TOPPBase("MzTabExporter", "Exports various XML formats to an mzTab file.")
    {
    }

protected:

    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "Input files used to generate the mzTab file.", false);
      setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML,idXML,mzid"));
      registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
      setValidFormats_("out", ListUtils::create<String>("tsv"));
    }


    /**
      @brief Gets peptide_evidences with data from internal structures adds their info to an MzTabPSMSectionRow (pre- or unfilled)

      @param peptide_evidences Vector of PeptideEvidence holding internal data.
      @param row Pre- or unfilled MzTabPSMSectionRow to be filled with the data.
      @param rows Vector of MzTabPSMSectionRow to add the differently updated rows to.
    */
    static void addPepEvidenceToRows(const vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows)
    {
      if (!peptide_evidences.empty())
      {
        for (Size i = 0; i != peptide_evidences.size(); ++i)
        {
          // get AABefore and AAAfter as well as start and end for all pep evidences

          // pre/post
          // from spec: Amino acid preceding the peptide (coming from the PSM) in the protein
          // sequence. If unknown “null” MUST be used, if the peptide is N-terminal “-“
          // MUST be used.
          if (peptide_evidences[i].getAABefore() == PeptideEvidence::UNKNOWN_AA)
          {
            row.pre = MzTabString("null");
          }
          else if (peptide_evidences[i].getAABefore() == PeptideEvidence::N_TERMINAL_AA)
          {
            row.pre = MzTabString("-");
          }
          else
          {
            row.pre = MzTabString(String(peptide_evidences[i].getAABefore()));
          }

          if (peptide_evidences[i].getAAAfter() == PeptideEvidence::UNKNOWN_AA)
          {
            row.post = MzTabString("null");
          }
          else if (peptide_evidences[i].getAAAfter() == PeptideEvidence::C_TERMINAL_AA)
          {
            row.post = MzTabString("-");
          }
          else
          {
            row.post = MzTabString(String(peptide_evidences[i].getAAAfter()));
          }

          // start/end
          if (peptide_evidences[i].getStart() == PeptideEvidence::UNKNOWN_POSITION)
          {
            row.start = MzTabString("null");
          }
          else
          {
            row.start = MzTabString(String(peptide_evidences[i].getStart() + 1)); // counting in mzTab starts at 1
          }

          if (peptide_evidences[i].getEnd() == PeptideEvidence::UNKNOWN_POSITION)
          {
            row.end = MzTabString("null");
          }
          else
          {
            row.end = MzTabString(String(peptide_evidences[i].getEnd() + 1)); // counting in mzTab starts at 1
          }

          row.accession = MzTabString(peptide_evidences[i].getProteinAccession());

          rows.push_back(row);
        }
      }
      else
      { // report without pep evidence information
        row.pre = MzTabString("null");
        row.post = MzTabString("null");
        row.start = MzTabString("null");
        row.end = MzTabString("null");
        rows.push_back(row);
      }
    }

    /**
      @brief Inserts values from MetaInfoInterface objects matching a (precalculated or filtered) set of keys to optional columns of an MzTab row.

      @param keys Only values matching those keys will be extracted from the object inheriting from MetaInfoInterface.
      @param opt A vector of optional columns to add to.
      @param id The identifier for this optional value according to the mzTab standard (like global, MS_Run, assay, etc.)
      @param meta The object holding the MetaInformation (like PeptideHit, ProteinHit, etc.)
      @return void: Only updates the values of the columns in opt
    */
    static void addMetaInfoToOptionalColumns(const set<String>& keys, vector<MzTabOptionalColumnEntry>& opt, const String id, const MetaInfoInterface meta)
    {
      for (set<String>::const_iterator sit = keys.begin(); sit != keys.end(); ++sit)
      {
        const String& key = *sit;
        MzTabOptionalColumnEntry opt_entry;
        opt_entry.first = String("opt_") + id + String("_") + String(key).substitute(' ','_');
        if (meta.metaValueExists(key))
        {
          opt_entry.second = MzTabString(meta.getMetaValue(key).toString().substitute(' ','_'));
        } // otherwise it is default ("null")
        opt.push_back(opt_entry);
      }
    }

    static map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const vector<String>& mods)
    {
      map<Size, MzTabModificationMetaData> mods_mztab;
      for (vector<String>::const_iterator sit = mods.begin(); sit != mods.end(); ++sit)
      {
        Size index = (sit - mods.begin()) + 1;
        MzTabModificationMetaData mod;
        MzTabParameter mp;
        mp.setCVLabel("UNIMOD");
        ModificationsDB* mod_db = ModificationsDB::getInstance();
        // MzTab standard is to just report Unimod accession.
        ResidueModification m = mod_db->getModification(*sit);
        String unimod_accession = m.getUniModAccession();
        mp.setAccession(unimod_accession.toUpper());
        mp.setName(m.getId());
        mod.modification = mp;

        if (m.getTermSpecificity() == ResidueModification::C_TERM)
        {
          mod.position = MzTabString("Any C-term");
        }
        else if (m.getTermSpecificity() == ResidueModification::N_TERM)
        {
          mod.position = MzTabString("Any N-term");
        }
        else if (m.getTermSpecificity() == ResidueModification::ANYWHERE)
        {
          mod.position = MzTabString("Anywhere");
        }
        else if (m.getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
        {
          mod.position = MzTabString("Protein C-term");
        }
        else if (m.getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
        {
          mod.position = MzTabString("Protein N-term");
        }

        mod.site = MzTabString(m.getOrigin());
        mods_mztab[index] = mod;
      }
      return mods_mztab;
    }

    static map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const vector<String>& mods)
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

    static map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const vector<String>& mods)
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
 
    static MzTab exportFeatureMapToMzTab(const FeatureMap& feature_map, const String& filename)
    {
      LOG_INFO << "exporting feature map: \"" << filename << "\" to mzTab: " << std::endl;
      MzTab mztab;
      MzTabMetaData meta_data;

      vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        
      vector<String> var_mods, fixed_mods;
      MzTabString db, db_version;
      if (!prot_ids.empty())
      {
        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
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
      meta_data.description = MzTabString("Export from featureXML");

      MzTabMSRunMetaData ms_run;
      StringList spectra_data;
      feature_map.getPrimaryMSRunPath(spectra_data);
      ms_run.location = spectra_data.empty() ? MzTabString("null") : MzTabString(spectra_data[0]);
      meta_data.ms_run[1] = ms_run;
      meta_data.uri[1] = MzTabString(filename);
      meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO: we currently only support psm search engine scores annotated to the identification run
      meta_data.peptide_search_engine_score[1] = MzTabParameter();

      mztab.setMetaData(meta_data);

      // pre-analyze data for occuring meta values at feature and peptide hit level
      // these are used to build optional columns containing the meta values in internal data structures

      set<String> feature_user_value_keys;
      set<String> peptide_hit_user_value_keys;
      for (Size i = 0; i < feature_map.size(); ++i)
      {
        const Feature& f = feature_map[i];
        vector<String> keys;
        f.getKeys(keys); //TODO: why not just return it?
        feature_user_value_keys.insert(keys.begin(), keys.end());

        const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
        for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            vector<String> ph_keys;
            hit->getKeys(ph_keys);
            peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
          }
        }
      }

      MzTabPeptideSectionRows rows;
      for (Size i = 0; i < feature_map.size(); ++i)
      {
        MzTabPeptideSectionRow row;
        const Feature& f = feature_map[i];
        row.mass_to_charge = MzTabDouble(f.getMZ());
        MzTabDoubleList rt_list;
        vector<MzTabDouble> rts;
        rts.push_back(MzTabDouble(f.getRT()));
        rt_list.set(rts);
        row.retention_time = rt_list;

        // set rt window if a bounding box has been set
        vector<MzTabDouble> window;
        if (f.getConvexHull().getBoundingBox() != DBoundingBox<2>())
        {
          window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().minX()));
          window.push_back(MzTabDouble(f.getConvexHull().getBoundingBox().maxX()));
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
        opt_global_modified_sequence.first = String("opt_global_modified_sequence");
        row.opt_.push_back(opt_global_modified_sequence);

        // create and fill opt_ columns for feature (peptide) user values
        addMetaInfoToOptionalColumns(feature_user_value_keys, row.opt_, String("global"), f);

        const vector<PeptideIdentification>& pep_ids = f.getPeptideIdentifications();
        if (pep_ids.empty())
        {
          rows.push_back(row);
          continue;
        }

        // TODO: here we assume that all have the same score type etc.
        vector<PeptideHit> all_hits;
        for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
        }

        if (all_hits.empty())
        {
          rows.push_back(row);
          continue;
        }

        // create new peptide id object to assist in sorting
        PeptideIdentification new_pep_id = pep_ids[0];
        new_pep_id.setHits(all_hits);
        new_pep_id.assignRanks();

        const PeptideHit& best_ph = new_pep_id.getHits()[0];
        const AASequence& aas = best_ph.getSequence();
        row.sequence = MzTabString(aas.toUnmodifiedString());

        row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

        const set<String>& accessions = best_ph.extractProteinAccessionsSet();
        const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

        row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
        // select accession of first peptide_evidence as representative ("leading") accession
        row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());
        row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());
        row.search_engine_score_ms_run[1][1] = MzTabDouble(best_ph.getScore());

        // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
        for (Size i = 0; i != row.opt_.size(); ++i)
        {
          MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

          if (opt_entry.first == String("opt_global_modified_sequence"))
          {
            opt_entry.second = MzTabString(aas.toString());
          }
        }

        // create and fill opt_ columns for psm (PeptideHit) user values
        addMetaInfoToOptionalColumns(peptide_hit_user_value_keys, row.opt_, String("global"), best_ph);

        rows.push_back(row);
      }
      mztab.setPeptideSectionRows(rows);
      return mztab;
    }

    static MzTab exportIdentificationsToMzTab(const vector<ProteinIdentification>& prot_ids, const vector<PeptideIdentification>& peptide_ids, const String& filename)
    {
      LOG_INFO << "exporting identifications: \"" << filename << "\" to mzTab: " << std::endl;
      vector<PeptideIdentification> pep_ids = peptide_ids;
      MzTab mztab;
      MzTabMetaData meta_data;
      vector<String> var_mods, fixed_mods;
      MzTabString db, db_version;
      String search_engine;
      String search_engine_version;

      if (!prot_ids.empty())
      {
        search_engine = prot_ids[0].getSearchEngine();
        search_engine_version = prot_ids[0].getSearchEngineVersion();
      }

      if (!prot_ids.empty())
      {
        MzTabParameter protein_score_type;
        protein_score_type.fromCellString("[,,custom score,]"); // TODO at least it should be noted if higher score is better. Better document type of score
        meta_data.protein_search_engine_score[1] = protein_score_type; // TODO add meta value to ProteinIdentification
        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
        var_mods = sp.variable_modifications;
        fixed_mods = sp.fixed_modifications;
        db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
        db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);

        //sp.digestion_enzyme
        //sp.missed_cleavages
        // generate protein section
        MzTabProteinSectionRows protein_rows;

        Size current_run_index(1);
        for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin();
         it != prot_ids.end(); ++it, ++current_run_index)
        {
          const std::vector<ProteinIdentification::ProteinGroup> protein_groups = it->getProteinGroups();
          const std::vector<ProteinIdentification::ProteinGroup> indist_groups = it->getIndistinguishableProteins();
          const std::vector<ProteinHit> protein_hits = it->getHits();

          MzTabMSRunMetaData ms_run;
          StringList ms_run_in_data;
          it->getPrimaryMSRunPath(ms_run_in_data);
          ms_run.location = ms_run_in_data.empty() ? MzTabString("null") : MzTabString(ms_run_in_data[0]);
          // TODO: add processing information that this file has been exported from "filename"
          meta_data.ms_run[current_run_index] = ms_run;

          // pre-analyze data for occuring meta values at protein hit level
          // these are used to build optional columns containing the meta values in internal data structures
          set<String> protein_hit_user_value_keys = 
            MetaInfoInterfaceUtils::findCommonMetaKeys<vector<ProteinHit>, set<String> >(protein_hits.begin(), protein_hits.end(), 100.0);

          // we do not want descriptions twice
          protein_hit_user_value_keys.erase("Description");

          for (Size i = 0; i != protein_hits.size(); ++i)
          {
            const ProteinHit& hit = protein_hits[i];
            MzTabProteinSectionRow protein_row;

            protein_row.accession = MzTabString(hit.getAccession());
            protein_row.description = MzTabString(hit.getDescription()); 
         // protein_row.taxid = hit.getTaxonomyID(); // TODO add as meta value to protein hitNEWT taxonomy for the species.
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
         // MzTabModificationList modifications; // Modifications identified in the protein.
         // MzTabString uri; // Location of the protein’s source entry.
         // MzTabStringList go_terms; // List of GO terms for the protein.
            double coverage = hit.getCoverage();
            protein_row.protein_coverage = coverage >= 0 ? MzTabDouble(coverage) : MzTabDouble(); // (0-1) Amount of protein sequence identified.
         // std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”

            // create and fill opt_ columns for protein hit user values
            addMetaInfoToOptionalColumns(protein_hit_user_value_keys, protein_row.opt_, String("global"), hit);

            // optional column for protein groups
            MzTabOptionalColumnEntry opt_column_entry;
            opt_column_entry.first = "opt_global_protein_group_type";
            opt_column_entry.second = MzTabString("single_protein");
            protein_row.opt_.push_back(opt_column_entry);
             
            protein_rows.push_back(protein_row);
          }

          // Protein groups are currently simply PRT rows with extra opt columns
          for (Size i = 0; i != protein_groups.size(); ++i)
          {
            const ProteinIdentification::ProteinGroup& group = protein_groups[i];
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
              entries.push_back(MzTabString(group.accessions[j]));
            }
            ambiguity_members.set(entries);
            protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
            protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);

            MzTabOptionalColumnEntry opt_column_entry;
            opt_column_entry.first = "opt_global_protein_group_type";
            opt_column_entry.second = MzTabString("protein_group");
            protein_row.opt_.push_back(opt_column_entry);
            protein_rows.push_back(protein_row);
          }

          for (Size i = 0; i != indist_groups.size(); ++i)
          {
            const ProteinIdentification::ProteinGroup& group = indist_groups[i];
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
              }
              entries.push_back(MzTabString(group.accessions[j]));
            }
            ambiguity_members.set(entries);
            protein_row.ambiguity_members = ambiguity_members; // Alternative protein identifications.
            MzTabOptionalColumnEntry opt_column_entry;
            opt_column_entry.first = "opt_global_protein_group_type";
            opt_column_entry.second = MzTabString("indistinguishable_group");
            protein_row.opt_.push_back(opt_column_entry);
            protein_row.best_search_engine_score[1] = MzTabDouble(group.probability);
            //std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”
            protein_rows.push_back(protein_row);
          }

        }
        mztab.setProteinSectionRows(protein_rows);
      }
      // end protein groups

      // start PSMs

      // mandatory meta values
      meta_data.mz_tab_type = MzTabString("Identification");
      meta_data.mz_tab_mode = MzTabString("Summary");
      meta_data.description = MzTabString("Export from idXML");

      meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
      meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);
      MzTabParameter psm_search_engine_score;
      psm_search_engine_score.fromCellString("[,," + search_engine + "," + search_engine_version + "]");
      meta_data.psm_search_engine_score[1] = psm_search_engine_score;

      mztab.setMetaData(meta_data);

      MzTabPSMSectionRows rows;
      Size psm_id(0);
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it, ++psm_id)
      {
        // skip empty peptide identification objects
        if (it->getHits().empty())
        {
          continue;
        }

        // sort by rank
        it->assignRanks();

        MzTabPSMSectionRow row;

        // only consider best peptide hit for export
        const PeptideHit& best_ph = it->getHits()[0];
        const AASequence& aas = best_ph.getSequence();
        row.sequence = MzTabString(aas.toUnmodifiedString());

        // extract all modifications in the current sequence for reporting. In contrast to peptide and protein section all modifications are reported.
        row.modifications = extractModificationListFromAASequence(aas);

        row.PSM_ID = MzTabInteger(psm_id);
        row.database = db;
        row.database_version = db_version;
        MzTabParameterList search_engines;
        search_engines.fromCellString("[,," + search_engine + "," + search_engine_version + "]");
        row.search_engine = search_engines;

        row.search_engine_score[1] = MzTabDouble(best_ph.getScore());
        vector<MzTabDouble> rts_vector;
        rts_vector.push_back(MzTabDouble(it->getRT()));
        MzTabDoubleList rts;
        rts.set(rts_vector);
        row.retention_time = rts;
        row.charge = MzTabInteger(best_ph.getCharge());
        row.exp_mass_to_charge = MzTabDouble(it->getMZ());
        row.calc_mass_to_charge = best_ph.getCharge() != 0 ? MzTabDouble(aas.getMonoWeight(Residue::Full, best_ph.getCharge()) / best_ph.getCharge()) : MzTabDouble();

        // add opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
        MzTabOptionalColumnEntry opt_entry;
        opt_entry.first = String("opt_global_modified_sequence");
        opt_entry.second = MzTabString(aas.toString());
        row.opt_.push_back(opt_entry);

        // currently write all keys
        // TODO: percentage procedure with MetaInfoInterfaceUtils
        vector<String> ph_keys;
        best_ph.getKeys(ph_keys);
        // TODO: no conversion but make function on collections
        set<String> ph_key_set(ph_keys.begin(), ph_keys.end());
        addMetaInfoToOptionalColumns(ph_key_set, row.opt_, String("global"), best_ph);

        // TODO Think about if the uniqueness can be determined by # of peptide evidences
        // b/c this would only differ when evidences come from different DBs
        const set<String>& accessions = best_ph.extractProteinAccessionsSet();
        row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);

        // create row for every PeptideEvidence entry (mapping to a protein)
        const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

        // pass common row entries and create rows for all peptide evidences
        addPepEvidenceToRows(peptide_evidences, row, rows);
      }

      mztab.setPSMSectionRows(rows);

      return mztab;
    }

    // Generate MzTab style list of PTMs from AASequence object. 
    // All passed fixed modifications are not reported (as suggested by the standard for the PRT and PEP section).
    // In contrast, all modifications are reported in the PSM section (see standard document for details).
    static MzTabModificationList extractModificationListFromAASequence(const AASequence& aas, const vector<String>& fixed_mods = vector<String>())
    {
      MzTabModificationList mod_list;
      vector<MzTabModification> mods;

      if (aas.isModified())
      {
        if (aas.hasNTerminalModification())
        {
          MzTabModification mod;
          const ResidueModification& res_mod = *(aas.getNTerminalModification());
          if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
          {
            MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
            vector<std::pair<Size, MzTabParameter> > pos;
            pos.push_back(make_pair(0, MzTabParameter()));
            mod.setModificationIdentifier(unimod_accession);
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
            if (std::find(fixed_mods.begin(), fixed_mods.end(), res_mod.getId()) == fixed_mods.end())
            {
              // MzTab standard is to just report Unimod accession.
              MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
              vector<std::pair<Size, MzTabParameter> > pos;
              pos.push_back(make_pair(ai + 1, MzTabParameter()));
              mod.setPositionsAndParameters(pos);
              mod.setModificationIdentifier(unimod_accession);
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
            MzTabString unimod_accession = MzTabString(res_mod.getUniModAccession());
            vector<std::pair<Size, MzTabParameter> > pos;
            pos.push_back(make_pair(aas.size() + 1, MzTabParameter()));
            mod.setPositionsAndParameters(pos);
            mod.setModificationIdentifier(unimod_accession);
            mods.push_back(mod);
          }
        }
      }
      mod_list.set(mods);
      return mod_list;
    }

    static MzTab exportConsensusMapToMzTab(const ConsensusMap& consensus_map, const String& filename)
    {
      LOG_INFO << "exporting consensus map: \"" << filename << "\" to mzTab: " << std::endl;
      MzTab mztab;
      vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();
      vector<String> var_mods, fixed_mods;
      MzTabString db, db_version;
      if (!prot_ids.empty())
      {
        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
        var_mods = sp.variable_modifications;
        fixed_mods = sp.fixed_modifications;
        db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
        db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);
      }

      // determine number of channels
      Size n_study_variables = consensus_map.getFileDescriptions().size();

      MzTabMetaData meta_data;

      // mandatory meta values
      meta_data.mz_tab_type = MzTabString("Quantification");
      meta_data.mz_tab_mode = MzTabString("Summary");
      meta_data.description = MzTabString("Export from consensusXML");

      // For consensusXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
      // but only at the study variable variable level.

      meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
      meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);
      meta_data.peptide_search_engine_score[1] = MzTabParameter();
      meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO insert search engine information
      MzTabMSRunMetaData ms_run;
      StringList ms_runs;
      consensus_map.getPrimaryMSRunPath(ms_runs);
      for (Size i = 0; i != ms_runs.size(); ++i)
      {
        ms_run.location = MzTabString(ms_runs[i]);
        meta_data.ms_run[i + 1] = ms_run;
      }

      mztab.setMetaData(meta_data);

      // pre-analyze data for occurring meta values at consensus feature and peptide hit level
      // these are used to build optional columns containing the meta values in internal data structures
      set<String> consensus_feature_user_value_keys;
      set<String> peptide_hit_user_value_keys;
      for (Size i = 0; i < consensus_map.size(); ++i)
      {
        const ConsensusFeature& c = consensus_map[i];
        vector<String> keys;
        c.getKeys(keys);
        consensus_feature_user_value_keys.insert(keys.begin(), keys.end());

        const vector<PeptideIdentification>& pep_ids = c.getPeptideIdentifications();
        for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            vector<String> ph_keys;
            hit->getKeys(ph_keys);
            peptide_hit_user_value_keys.insert(ph_keys.begin(), ph_keys.end());
          }
        }
      }

      MzTabPeptideSectionRows rows;

      for (Size i = 0; i < consensus_map.size(); ++i)
      {
        MzTabPeptideSectionRow row;
        const ConsensusFeature& c = consensus_map[i];

        // create opt_ column for peptide sequence containing modification
        MzTabOptionalColumnEntry opt_global_modified_sequence;
        opt_global_modified_sequence.first = String("opt_global_modified_sequence");
        row.opt_.push_back(opt_global_modified_sequence);

        // create opt_ columns for consensus feature (peptide) user values
        for (set<String>::const_iterator mit = consensus_feature_user_value_keys.begin(); mit != consensus_feature_user_value_keys.end(); ++mit)
        {
          MzTabOptionalColumnEntry opt_entry;
          const String& key = *mit;
          opt_entry.first = String("opt_global_") + key;
          if (c.metaValueExists(key))
          {
            opt_entry.second = MzTabString(c.getMetaValue(key).toString());
          } // otherwise it is default ("null")
          row.opt_.push_back(opt_entry);
        }

        // create opt_ columns for psm (PeptideHit) user values
        for (set<String>::const_iterator mit = peptide_hit_user_value_keys.begin(); mit != peptide_hit_user_value_keys.end(); ++mit)
        {
          MzTabOptionalColumnEntry opt_entry;
          const String& key = *mit;
          opt_entry.first = String("opt_global_") + key;
          // leave value empty as we have to fill it with the value from the best peptide hit
          row.opt_.push_back(opt_entry);
        }

        row.mass_to_charge = MzTabDouble(c.getMZ());
        MzTabDoubleList rt_list;
        vector<MzTabDouble> rts;
        rts.push_back(MzTabDouble(c.getRT()));
        rt_list.set(rts);
        row.retention_time = rt_list;
        MzTabDoubleList rt_window;
        row.retention_time_window = rt_window;
        row.charge = MzTabInteger(c.getCharge());
        row.best_search_engine_score[1] = MzTabDouble();

        // initialize columns
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
        for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
        {
          Size study_variable = fit->getMapIndex() + 1;
          row.peptide_abundance_stdev_study_variable[study_variable];
          row.peptide_abundance_std_error_study_variable[study_variable];
          row.peptide_abundance_study_variable[study_variable] = MzTabDouble(fit->getIntensity());
        }

        vector<PeptideIdentification> pep_ids = c.getPeptideIdentifications();
        if (!pep_ids.empty())
        {
          if (pep_ids.size() != 1)
          {
            throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Consensus features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!");
          }

          pep_ids[0].assignRanks();
          const PeptideHit& best_ph = pep_ids[0].getHits()[0];
          const AASequence& aas = best_ph.getSequence();
          row.sequence = MzTabString(aas.toUnmodifiedString());

          row.modifications = extractModificationListFromAASequence(aas, fixed_mods);

          const set<String>& accessions = best_ph.extractProteinAccessionsSet();
          const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

          row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);
          // select accession of first peptide_evidence as representative ("leading") accession
          row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());

          row.best_search_engine_score[1] = MzTabDouble(best_ph.getScore());

          // TODO: support run level scores - for now we assume we got the same score from every ms run
          for (Size ms_run = 1; ms_run <= ms_runs.size(); ++ms_run)
          {
            row.search_engine_score_ms_run[1][ms_run] = MzTabDouble(best_ph.getScore());
          }

          // fill opt_ columns

          // find opt_global_modified_sequence in opt_ and set it to the OpenMS amino acid string (easier human readable than unimod accessions)
          for (Size i = 0; i != row.opt_.size(); ++i)
          {
            MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

            if (opt_entry.first == String("opt_global_modified_sequence"))
            {
              opt_entry.second = MzTabString(aas.toString());
            }
          }

          // fill opt_ column of psm
          vector<String> ph_keys;
          best_ph.getKeys(ph_keys);
          for (Size k = 0; k != ph_keys.size(); ++k)
          {
            const String& key = ph_keys[k];

            // find matching entry in opt_ (TODO: speed this up)
            for (Size i = 0; i != row.opt_.size(); ++i)
            {
              MzTabOptionalColumnEntry& opt_entry = row.opt_[i];

              if (opt_entry.first == String("opt_global_") + key)
              {
                opt_entry.second = MzTabString(best_ph.getMetaValue(key).toString());
              }
            }
          }
        }

        rows.push_back(row);
      }

      mztab.setPeptideSectionRows(rows);
      return mztab;
    }

    ExitCodes main_(int, const char**) override
    {
      // parameter handling
      String in = getStringOption_("in");
      FileTypes::Type in_type = FileHandler().getType(in);

      String out = getStringOption_("out");

      MzTab mztab;

      if (in_type == FileTypes::FEATUREXML)
      {
        // For featureXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
        // but only at the (single) study variable variable level.

        // load featureXML
        FeatureMap feature_map;
        FeatureXMLFile f;
        f.load(in, feature_map);

        // calculate coverage
        vector<PeptideIdentification> pep_ids;
        vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        

        for (Size i = 0; i < feature_map.size(); ++i) // collect all (assigned and unassigned to a feature) peptide ids
        {
          vector<PeptideIdentification> pep_ids_bf = feature_map[i].getPeptideIdentifications();
          pep_ids.insert(pep_ids.end(), pep_ids_bf.begin(), pep_ids_bf.end());
        }

        pep_ids.insert(pep_ids.end(), feature_map.getUnassignedPeptideIdentifications().begin(), feature_map.getUnassignedPeptideIdentifications().end());

        try // might throw Exception::MissingInformation()
        {
          for (Size i = 0; i < prot_ids.size(); ++i)
          {
            prot_ids[i].computeCoverage(pep_ids);
          }
        }
        catch (Exception::MissingInformation& e)
        {
          LOG_WARN << "Non-critical exception: " << e.what() << "\n";
        }
        feature_map.setProteinIdentifications(prot_ids);

        mztab = exportFeatureMapToMzTab(feature_map, in);
      }

      // export identification data from idXML
      if (in_type == FileTypes::IDXML)
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        IdXMLFile().load(in, prot_ids, pep_ids, document_id);
        mztab = exportIdentificationsToMzTab(prot_ids, pep_ids, in); 
      }

      // export identification data from mzIdentML
      if (in_type == FileTypes::MZIDENTML)
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        MzIdentMLFile().load(in, prot_ids, pep_ids);
        mztab = exportIdentificationsToMzTab(prot_ids, pep_ids, in); 
      }

      // export quantification data
      if (in_type == FileTypes::CONSENSUSXML)
      {
        ConsensusMap consensus_map;
        ConsensusXMLFile c;
        c.load(in, consensus_map);
        mztab = exportConsensusMapToMzTab(consensus_map, in);
      }

      MzTabFile().store(out, mztab);
      return EXECUTION_OK;
    }
  };
} //namespace OpenMS

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  TOPPMzTabExporter t;
  return t.main(argc, argv);
}

/// @endcond

