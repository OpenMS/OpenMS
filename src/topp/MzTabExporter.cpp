// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
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
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> ProteinQuantifier </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

  @experimental This algorithm and underlying format is work in progress and might change.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MzTabExporter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MzTabExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

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

    void registerOptionsAndFlags_()
    {
      registerInputFile_("in_feature", "<file>", "", "FeatureXMLs used to generate the mzTab file.", false);
      setValidFormats_("in_feature", ListUtils::create<String>("featureXML"));
      registerInputFile_("in_consensus", "<file>", "", "ConsensusXMLs used to generate the mzTab file.", false);
      setValidFormats_("in_consensus", ListUtils::create<String>("consensusXML"));
      registerInputFile_("in_id", "<file>", "", "Identifications used to generate the mzTab file.", false);
      setValidFormats_("in_id", ListUtils::create<String>("idXML"));
      registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
      setValidFormats_("out", ListUtils::create<String>("csv"));
    }

    map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const vector<String>& mods)
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

    map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const vector<String>& mods)
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

    map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const vector<String>& mods)
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


    ExitCodes main_(int, const char **)
    {
      // parameter handling
      String in_feature = getStringOption_("in_feature");
      String in_id = getStringOption_("in_id");
      String in_consensus = getStringOption_("in_consensus");
      String out = getStringOption_("out");

      if (!in_feature.empty())
      {
        MzTab mztab;
        MzTabMetaData meta_data;

        // For featureXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
        // but only at the (single) study variable variable level.

        // load featureXML
        FeatureMap feature_map;
        FeatureXMLFile f;
        f.load(in_feature, feature_map);

        vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();
        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
        vector<String> var_mods = sp.variable_modifications;
        vector<String> fixed_mods = sp.fixed_modifications;

        meta_data.variable_mod = generateMzTabStringFromVariableModifications(var_mods);
        meta_data.fixed_mod = generateMzTabStringFromFixedModifications(fixed_mods);

        vector<PeptideIdentification> pep_ids;

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
        catch (Exception::MissingInformation & e)
        {
          LOG_WARN << "Non-critical exception: " << e.what() << "\n";
        }
        feature_map.setProteinIdentifications(prot_ids);

        // mandatory meta values
        meta_data.mz_tab_type = MzTabString("Quantification");
        meta_data.mz_tab_mode = MzTabString("Summary");
        meta_data.description = MzTabString("Export from featureXML");

        MzTabMSRunMetaData ms_run;
        ms_run.location = MzTabString("null"); // TODO: file origin of ms run (e.g. mzML) not stored in featureXML so far
        meta_data.ms_run[1] = ms_run;
        meta_data.uri[1] = MzTabString(in_feature);
        meta_data.psm_search_engine_score[1] = MzTabParameter();  // TODO: we currently only support psm search engine scores annotated to the identification run
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
          f.getKeys(keys);  //TODO: why not just return it?
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

          // create opt_ columns for feature (peptide) user values
          for (set<String>::const_iterator mit = feature_user_value_keys.begin(); mit != feature_user_value_keys.end(); ++mit)
          {
            MzTabOptionalColumnEntry opt_entry;
            const String& key = *mit;
            opt_entry.first = String("opt_peptide_") + key;
            if (f.metaValueExists(key))
            {
              opt_entry.second = MzTabString(f.getMetaValue(key).toString());
            } // otherwise it is default ("null")
            row.opt_.push_back(opt_entry);
          }

          // create opt_ columns for psm (PeptideHit) user values
          for (set<String>::const_iterator mit = peptide_hit_user_value_keys.begin(); mit != peptide_hit_user_value_keys.end(); ++mit)
          {
            MzTabOptionalColumnEntry opt_entry;
            const String& key = *mit;
            opt_entry.first = String("opt_psm_") + key;
            // leave value empty as we have to fill it with the value from the best peptide hit
            row.opt_.push_back(opt_entry);
          }

          vector<PeptideIdentification> pep_ids = f.getPeptideIdentifications();
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

          MzTabModificationList mod_list;
          vector<MzTabModification> mods;
          if (aas.isModified())
          {
            for (Size ai = 0; ai != aas.size(); ++ai)
            {
              if (aas.isModified(ai))
              {
                MzTabModification mod;
                String mod_name = aas[ai].getModification();
                ModificationsDB* mod_db = ModificationsDB::getInstance();

                // MzTab standard is to just report Unimod accession.
                MzTabString unimod_accession = MzTabString(mod_db->getModification(mod_name).getUniModAccession());
                mod.setModificationIdentifier(unimod_accession);
                vector<std::pair<Size, MzTabParameter> > pos;
                pos.push_back(make_pair(ai + 1, MzTabParameter()));
                mod.setPositionsAndParameters(pos);
                mods.push_back(mod);
              }
            }
          }
          mod_list.set(mods);
          row.modifications = mod_list;

          const set<String>& accessions = best_ph.extractProteinAccessions();
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

              if (opt_entry.first == String("opt_psm_") + key)
              {
                opt_entry.second = MzTabString(best_ph.getMetaValue(key).toString());
              }
            }
          }
          rows.push_back(row);
        }
        mztab.setPeptideSectionRows(rows);
        MzTabFile().store(out, mztab);
        return EXECUTION_OK;
      }

      // export identification data
      if (!in_id.empty())
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        IdXMLFile().load(in_id, prot_ids, pep_ids, document_id);

        if (prot_ids.empty())
        {
          //TODO
        }

        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
        vector<String> var_mods = sp.variable_modifications;
        vector<String> fixed_mods = sp.fixed_modifications;
        MzTabString db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
        MzTabString db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);

        MzTab mztab;
        MzTabMetaData meta_data;

        // mandatory meta values
        meta_data.mz_tab_type = MzTabString("Identification");
        meta_data.mz_tab_mode = MzTabString("Summary");
        meta_data.description = MzTabString("Export from idXML");

        meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
        meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);

        meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO insert search engine information
        MzTabMSRunMetaData ms_run;
        ms_run.location = MzTabString(in_id);
        meta_data.ms_run[1] = ms_run;
        mztab.setMetaData(meta_data);

        MzTabPSMSectionRows rows;
        for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
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

          // extract all modifications in the current sequence for reporting
          MzTabModificationList mod_list;
          vector<MzTabModification> mods;
          if (aas.isModified())
          {
            for (Size ai = 0; ai != aas.size(); ++ai)
            {
              if (aas.isModified(ai))
              {
                MzTabModification mod;
                String mod_name = aas[ai].getModification();
                ModificationsDB* mod_db = ModificationsDB::getInstance();

                // MzTab standard is to just report Unimod accession.
                MzTabString unimod_accession = MzTabString(mod_db->getModification(mod_name).getUniModAccession());
                mod.setModificationIdentifier(unimod_accession);
                vector<std::pair<Size, MzTabParameter> > pos;
                pos.push_back(make_pair(ai + 1, MzTabParameter()));
                mod.setPositionsAndParameters(pos);
                mods.push_back(mod);
              }
            }
          }
          mod_list.set(mods);
          row.modifications = mod_list;

          const set<String>& accessions = best_ph.extractProteinAccessions();
          const vector<PeptideEvidence> peptide_evidences = best_ph.getPeptideEvidences();

          // determine if peptide unique (TODO: move to static helper)
          row.unique = accessions.size() == 1 ? MzTabBoolean(true) : MzTabBoolean(false);

          // TODO: add option to export all peptide evidences of a peptide (e.g. same sequence but different proteins)
          // select accession of first peptide_evidence as representative ("leading") accession
          row.accession = peptide_evidences.empty() ? MzTabString("null") : MzTabString(peptide_evidences[0].getProteinAccession());

          row.database = db;
          row.database_version = db_version;
          row.search_engine_score[1] = MzTabDouble(best_ph.getScore());
          vector<MzTabDouble> rts_vector;
          rts_vector.push_back(MzTabDouble(it->getRT()));
          MzTabDoubleList rts;
          rts.set(rts_vector);
          row.retention_time = rts;
          row.charge = MzTabInteger(best_ph.getCharge());
          row.exp_mass_to_charge = MzTabDouble(it->getMZ());
          row.calc_mass_to_charge = best_ph.getCharge() != 0 ? MzTabDouble(aas.getMonoWeight(Residue::Full, best_ph.getCharge()) / best_ph.getCharge()) : MzTabDouble();

          // get AABefore and AAAfter as well as start and end of first peptide_evidence as representative
          if (peptide_evidences.empty())
          {
            row.pre = MzTabString("-");
            row.post = MzTabString("-");
            row.start = MzTabString("-");
            row.end = MzTabString("-");
          }
          else
          {
            // pre/post
            // from spec: Amino acid preceding the peptide (coming from the PSM) in the protein
            // sequence. If unknown “null” MUST be used, if the peptide is N-terminal “-“
            // MUST be used.
            if (peptide_evidences[0].getAABefore() == PeptideEvidence::UNKNOWN_AA)
            {
              row.pre = MzTabString("null");
            }
            else if (peptide_evidences[0].getAABefore() == PeptideEvidence::N_TERMINAL_AA)
            {
              row.pre = MzTabString("-");
            }
            else
            {
              row.pre = MzTabString(String(peptide_evidences[0].getAABefore()));
            }

            if (peptide_evidences[0].getAAAfter() == PeptideEvidence::UNKNOWN_AA)
            {
              row.post = MzTabString("null");
            }
            else if (peptide_evidences[0].getAAAfter() == PeptideEvidence::C_TERMINAL_AA)
            {
              row.post = MzTabString("-");
            }
            else
            {
              row.post = MzTabString(String(peptide_evidences[0].getAAAfter()));
            }

            // start/end
            if (peptide_evidences[0].getStart() == PeptideEvidence::UNKNOWN_POSITION)
            {
              row.start = MzTabString("-");
            }
            else
            {
              row.start = MzTabString(String(peptide_evidences[0].getStart() + 1));  // counting in mzTab starts at 1
            }

            if (peptide_evidences[0].getEnd() == PeptideEvidence::UNKNOWN_POSITION)
            {
              row.end = MzTabString("-");
            }
            else
            {
              row.end = MzTabString(String(peptide_evidences[0].getEnd() + 1)); // counting in mzTab starts at 1
            }
          }

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

              if (opt_entry.first == String("opt_psm_") + key)
              {
                opt_entry.second = MzTabString(best_ph.getMetaValue(key).toString());
              }
            }
          }
          rows.push_back(row);
        }

        mztab.setPSMSectionRows(rows);

        MzTabFile().store(out, mztab);
        return EXECUTION_OK;
      }

      if (!in_consensus.empty())
      {
        // load consensus map
        ConsensusMap consensus_map;
        ConsensusXMLFile c;
        c.load(in_consensus, consensus_map);

        vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();
        ProteinIdentification::SearchParameters sp = prot_ids[0].getSearchParameters();
        vector<String> var_mods = sp.variable_modifications;
        vector<String> fixed_mods = sp.fixed_modifications;
//        MzTabString db = sp.db.empty() ? MzTabString() : MzTabString(sp.db);
//        MzTabString db_version = sp.db_version.empty() ? MzTabString() : MzTabString(sp.db_version);

        MzTab mztab;
        MzTabMetaData meta_data;

        // mandatory meta values
        meta_data.mz_tab_type = MzTabString("Quantification");
        meta_data.mz_tab_mode = MzTabString("Summary");
        meta_data.description = MzTabString("Export from consensusXML");

        // For consensusXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
        // but only at the study variable variable level.

        meta_data.variable_mod = generateMzTabStringFromModifications(var_mods);
        meta_data.fixed_mod = generateMzTabStringFromModifications(fixed_mods);

        meta_data.psm_search_engine_score[1] = MzTabParameter(); // TODO insert search engine information
        MzTabMSRunMetaData ms_run;
        ms_run.location = MzTabString(in_id);
        meta_data.ms_run[1] = ms_run;
        mztab.setMetaData(meta_data);

        MzTabPeptideSectionRows rows;

        for (Size i = 0; i < consensus_map.size(); ++i)
        {
          MzTabPeptideSectionRow row;
          const ConsensusFeature& c = consensus_map[i];
          row.mass_to_charge = MzTabDouble(c.getMZ());
          MzTabDoubleList rt_list;
          vector<MzTabDouble> rts;
          rts.push_back(MzTabDouble(c.getRT()));
          rt_list.set(rts);
          row.retention_time = rt_list;
          MzTabDoubleList rt_window;
          row.retention_time_window = rt_window;
          row.charge = MzTabInteger(c.getCharge());

          ConsensusFeature::HandleSetType fs = c.getFeatures();
          Size index = 0;
          for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit, ++index)
          {
//            Size index = std::distance(fit, fs.begin());
            row.peptide_abundance_stdev_study_variable[index];
            row.peptide_abundance_std_error_study_variable[index];
            row.peptide_abundance_study_variable[index] = MzTabDouble(fit->getIntensity());
            row.best_search_engine_score[index] = MzTabDouble();
            row.search_engine_score_ms_run[index][1] = MzTabDouble();
          }

          rows.push_back(row);
        }

        MzTabFile().store(out, mztab);
        return EXECUTION_OK;
      }

      return EXECUTION_OK;
    }

  };
}

int main(int argc, const char ** argv)
{
  TOPPMzTabExporter t;
  return t.main(argc, argv);
}

/// @endcond
