// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQEvidenceExporter.h>/*
#include <string>
#include <fstream>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
*/
using namespace std;
using namespace OpenMS;


MQEvidence::MQEvidence(const string &file)
{
    file_ = fstream(file, fstream::out);

    //TODO: Header Funktion , die spalten prüft
}

MQEvidence::~MQEvidence() {
    file_.close();
}
/*
std::fstream MQEvidence::getFile() const
{
    return file_;
}*/

void MQEvidence::export_header()
{
    if(!file_.good())
    {
        //TODO: exception
    }
    file_.clear();
    file_ << "Sequenz" << "\t";
}

void MQEvidence::f_export(const FeatureMap &fm) {
    int intensity;
    int charge;
    double rt;
    double mz;
    for (Feature f : fm) {
        intensity = f.getIntensity();
        charge = f.getCharge();
        mz = f.getMZ();
        rt = f.getRT();
        auto v_pepI = f.getPeptideIdentifications();
        auto pepH = v_pepI[0].getHits()[0];
        for (auto pep : v_pepI) {

            pep.sort();
            auto pH = pep.getHits();
            if (pH[0].getScore() > pepH.getScore()) {
                pepH = pH[0];
            }

        }
        double score_pepH = pepH.getScore();
        auto seq = pepH.getSequence();
        pepH.getSequence();


        file_ << seq << "\t" << seq.getNTerminalModificationName() << seq << seq.getCTerminalModificationName() << "\t"
              << charge << "\t" << mz << "\t" << rt << "\t" << score_pepH << "\t" << intensity << "\t";
    }
}
//////////////////////////////////////////////////////////////////////////////////////


void MQEvidence::exportFeatureMapTotxt(const FeatureMap & feature_map) {
    //go trough all features
    for (Feature f : feature_map)
    {
        //Todo: Export row über die Features
        //row funktion

    }
}


//////////////////////////////////////////////////////////////////////////////////////

MzTabPeptideSectionRow MzTab::peptideSectionRowFromConsensusFeature_(
        const ConsensusFeature& c,
        const ConsensusMap& consensus_map,
        const StringList& ms_runs,
        const Size n_study_variables,
        const set<String>& consensus_feature_user_value_keys,
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
    const auto addUserValueToRowBy = [&row](function<void(const String &s, MzTabOptionalColumnEntry &entry)> f) -> function<void(const String &key)>
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
        study_variable = path_label_to_assay.at(pl); // for now, a study_variable is one assay

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
        const PeptideHit& best_ph = curr_pep_ids[0].getHits()[0];
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
                if (pep.metaValueExists("id_merge_index"))
                {
                    id_merge_index = pep.getMetaValue("id_merge_index");
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
                    row.spectra_ref.setSpecRef(pep.getMetaValue("spectrum_reference").toString());
                    row.spectra_ref.setMSFile(msfile_index);
                }
            }
        }
        row.best_search_engine_score[1] = MzTabDouble(best_score);
    }

    remapTargetDecoyPSMAndPeptideSection_(row.opt_);
    return row;
}

