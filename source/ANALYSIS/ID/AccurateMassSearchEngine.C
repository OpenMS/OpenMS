// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>



#include <OpenMS/SYSTEM/File.h>

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <fstream>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{

/// default constructor
AccurateMassSearchResult::AccurateMassSearchResult() :
    adduct_mass_(),
    query_mass_(),
    found_mass_(),
    charge_(),
    error_ppm_(),
    observed_rt_(),
    observed_intensity_(),
    individual_intensities_(),
    matching_index_(),
    source_feature_index_(),
    found_adduct_(),
    empirical_formula_(),
    matching_hmdb_ids_(),
    isotopes_sim_score_(-1.0)
{

}

/// default destructor
AccurateMassSearchResult::~AccurateMassSearchResult()
{

}

/// copy constructor
AccurateMassSearchResult::AccurateMassSearchResult(const AccurateMassSearchResult& source) :
    adduct_mass_(source.adduct_mass_),
    query_mass_(source.query_mass_),
    found_mass_(source.found_mass_),
    charge_(source.charge_),
    error_ppm_(source.error_ppm_),
    observed_rt_(source.observed_rt_),
    observed_intensity_(source.observed_intensity_),
    individual_intensities_(source.individual_intensities_),
    source_feature_index_(source.source_feature_index_),
    matching_index_(source.matching_index_),
    found_adduct_(source.found_adduct_),
    empirical_formula_(source.empirical_formula_),
    matching_hmdb_ids_(source.matching_hmdb_ids_),
    isotopes_sim_score_(source.isotopes_sim_score_)
{

}

/// assignment operator
AccurateMassSearchResult& AccurateMassSearchResult::operator=(const AccurateMassSearchResult& rhs)
{
    if (this == &rhs)
        return *this;

    adduct_mass_ = rhs.adduct_mass_;
    query_mass_ = rhs.query_mass_;
    found_mass_ = rhs.found_mass_;
    charge_ = rhs.charge_;
    error_ppm_ = rhs.error_ppm_;
    observed_rt_ = rhs.observed_rt_;
    observed_intensity_ = rhs.observed_intensity_;
    individual_intensities_ = rhs.individual_intensities_;
    matching_index_ = rhs.matching_index_;
    source_feature_index_ = rhs.source_feature_index_;
    found_adduct_ = rhs.found_adduct_;
    empirical_formula_ = rhs.empirical_formula_;
    matching_hmdb_ids_ = rhs.matching_hmdb_ids_;
    isotopes_sim_score_ = rhs.isotopes_sim_score_;

    return *this;
}

DoubleReal AccurateMassSearchResult::getAdductMass() const
{
    return adduct_mass_;
}

void AccurateMassSearchResult::setAdductMass(const DoubleReal& m)
{
    adduct_mass_ = m;
}

DoubleReal AccurateMassSearchResult::getQueryMass() const
{
    return query_mass_;
}

void AccurateMassSearchResult::setQueryMass(const DoubleReal& m)
{
    query_mass_ = m;
}

DoubleReal AccurateMassSearchResult::getFoundMass() const
{
    return found_mass_;
}

void AccurateMassSearchResult::setFoundMass(const DoubleReal& m)
{
    found_mass_ = m;
}

Int AccurateMassSearchResult::getCharge() const
{
    return charge_;
}

void AccurateMassSearchResult::setCharge(const Int& ch)
{
    charge_ = ch;
}

DoubleReal AccurateMassSearchResult::getErrorPPM() const
{
    return error_ppm_;
}

void AccurateMassSearchResult::setErrorPPM(const DoubleReal& ppm)
{
    error_ppm_ = ppm;
}

DoubleReal AccurateMassSearchResult::getObservedRT() const
{
    return observed_rt_;
}

void AccurateMassSearchResult::setObservedRT(const DoubleReal& rt)
{
    observed_rt_ = rt;
}

DoubleReal AccurateMassSearchResult::getObservedIntensity() const
{
    return observed_intensity_;
}

void AccurateMassSearchResult::setObservedIntensity(const DoubleReal& intensity)
{
    observed_intensity_ = intensity;
}

std::vector<DoubleReal> AccurateMassSearchResult::getIndividualIntensities() const
{
    return individual_intensities_;
}

void AccurateMassSearchResult::setIndividualIntensities(const std::vector<DoubleReal>& indiv_ints)
{
    individual_intensities_ = indiv_ints;
}

Size AccurateMassSearchResult::getMatchingIndex() const
{
    return matching_index_;
}

void AccurateMassSearchResult::setMatchingIndex(const Size& idx)
{
    matching_index_ = idx;
}

Size AccurateMassSearchResult::getSourceFeatureIndex() const
{
    return source_feature_index_;
}

void AccurateMassSearchResult::setSourceFeatureIndex(const Size& idx)
{
    source_feature_index_ = idx;
}

const String& AccurateMassSearchResult::getFoundAdduct() const
{
    return found_adduct_;
}

void AccurateMassSearchResult::setFoundAdduct(const String& add)
{
    found_adduct_ = add;
}

const String& AccurateMassSearchResult::getFormulaString() const
{
    return empirical_formula_;
}

void AccurateMassSearchResult::setEmpiricalFormula(const String& ep)
{
    empirical_formula_ = ep;
}

const std::vector<String>& AccurateMassSearchResult::getMatchingHMDBids() const
{
    return matching_hmdb_ids_;
}

void AccurateMassSearchResult::setMatchingHMDBids(const std::vector<String>& match_ids)
{
    matching_hmdb_ids_ = match_ids;
}

DoubleReal AccurateMassSearchResult::getIsotopesSimScore() const
{
    return isotopes_sim_score_;
}

void AccurateMassSearchResult::setIsotopesSimScore(const DoubleReal& sim_score)
{
    isotopes_sim_score_ = sim_score;
}

void AccurateMassSearchResult::outputResults() const
{
    std::cout << "adduct_mass: " << std::setprecision(8) << adduct_mass_ << "\n";
    std::cout << "query_mass: " << query_mass_ << "\n";
    std::cout << "found_mass: " << found_mass_ << "\n";
    std::cout << "charge: " << charge_ << "\n";
    std::cout << "error ppm: " << error_ppm_ << "\n";
    std::cout << "observed rt: " << observed_rt_ << "\n";
    std::cout << "observed intensity: " << observed_intensity_ << "\n";


    std::cout << "matching idx: " << matching_index_ << "\n";

    std::cout << "found_adduct_: " << found_adduct_ << "\n";
    std::cout << "emp. formula: " << empirical_formula_ << "\n";
    std::cout << "matching HMDB ids:";

    for (Size i = 0; i < matching_hmdb_ids_.size(); ++i)
    {
        std::cout << " " << matching_hmdb_ids_[i];
    }

    std::cout << "\n";
    std::cout << "isocheck sim score: " << isotopes_sim_score_ << std::endl;    // ensure endl used at the end (but not before! performance!)
}

AccurateMassSearchEngine::AccurateMassSearchEngine() :
    DefaultParamHandler("AccurateMassSearchEngine"), 
    ProgressLogger(),
    is_initialized_(false)
{
    defaults_.setValue("mass_error_value", 5.0, "Tolerance allowed for accurate mass search.");

    defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
    defaults_.setValidStrings("mass_error_unit", StringList::create(("ppm,Da")));

    defaults_.setValue("ionization_mode", "positive", "Positive or negative ionization mode? If 'auto' is used, the first feature of the input map must contain the meta-value 'scan_polarity'. If its missing, the tool will exit with error.");
    defaults_.setValidStrings("ionization_mode", StringList::create(("positive,negative,auto")));

    defaults_.setValue("isotopic_similarity", "false", "Computes a similarity score for each hit (only if the feature exhibits at least two isotopic mass traces).");
    defaults_.setValidStrings("isotopic_similarity", StringList::create(("false,true")));

    defaults_.setValue("report_mode", "all", "Results are reported in one of several modes: Either (all) matching hits, the (top3) scoring hits, or the (best) scoring hit.");
    defaults_.setValidStrings("report_mode", StringList::create(("all,top3,best")));

    defaults_.setValue("db:mapping", "CHEMISTRY/HMDBMappingFile.tsv", "Database input file, containing three tab-separated columns of mass, formula, identifier. "
                                                                      "If 'mass' is 0, it is re-computed from the molecular sum formula. "
                                                                      "By default CHEMISTRY/HMDBMappingFile.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("db:struct", "CHEMISTRY/HMDB2StructMapping.tsv", "Database input file, containing four tab-separated columns of identifier, name, SMILES, INCHI."
                                                                        "The identifier should match with mapping file. SMILES and INCHI are reported in the output, but not used otherwise. "
                                                                        "By default CHEMISTRY/HMDB2StructMapping.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("positive_adducts_file", "CHEMISTRY/PositiveAdducts.tsv", "This file contains the list of potential positive adducts that will be looked for in the database. "
                                                                                 "Edit the list if you wish to exclude/include adducts. "
                                                                                 "By default CHEMISTRY/PositiveAdducts.tsv in OpenMS/share is used! If empty, the default will be used.", StringList::create("advanced"));
    defaults_.setValue("negative_adducts_file", "CHEMISTRY/NegativeAdducts.tsv", "This file contains the list of potential negative adducts that will be looked for in the database. "
                                                                                 "Edit the list if you wish to exclude/include adducts. "
                                                                                 "By default CHEMISTRY/NegativeAdducts.tsv in OpenMS/share is used! If empty, the default will be used.", StringList::create("advanced"));


    defaultsToParam_();

    this->setLogType(CMD);
}

AccurateMassSearchEngine::~AccurateMassSearchEngine()
{
}

/// public methods

void AccurateMassSearchEngine::queryByMass(const DoubleReal& adduct_mass, const Int& adduct_charge, std::vector<AccurateMassSearchResult>& results)
{
    if (!is_initialized_) init_(); // parse DB

    // Depending on ion_mode_internal_, the file containing the rules for positive or negative adducts is loaded
    StringList::ConstIterator it_s, it_e;
    if (ion_mode_internal_ == "positive")
    {
      it_s = pos_adducts_.begin();  
      it_e = pos_adducts_.end();  
    }
    else if (ion_mode_internal_ == "negative")
    {
      it_s = neg_adducts_.begin();  
      it_e = neg_adducts_.end();  
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Ion mode cannot be set to '") + ion_mode_ + "'!");
    }
    
    
    for (StringList::ConstIterator it = it_s; it != it_e; ++it)
    {
        DoubleReal query_mass;
        Int charge;
        String adduct_name(*it);

        computeNeutralMassFromAdduct_(adduct_mass, adduct_name, query_mass, charge);

        // std::cout << "looking for " << pos_adducts_[adduct_idx] << std::endl;
        //if ((adduct_charge > 0) && (charge != adduct_charge))
        if (adduct_charge != 0 && (std::abs(adduct_charge) != std::abs(charge)))
        { // we ignore 0 charge, but anything else must match in absolute terms (absolute, since any FeatureFinder gives only positive charges, even for negative-mode spectra)
            continue;
        }

        // get potential hits as indices in masskey_table
        std::vector<Size> hit_idx;
        searchMass_(query_mass, hit_idx);


        //std::cerr << ion_mode_internal_ << " adduct: " << adduct_name << ", " << adduct_mass << " Da, " << query_mass << " qm(against DB), " << charge << " q\n"; 

        // store information from query hits in AccurateMassSearchResult objects
        for (Size i = 0; i < hit_idx.size(); ++i)
        {
            DoubleReal found_mass(mass_mappings_[hit_idx[i]].mass);
            DoubleReal found_error_ppm(((query_mass - found_mass) / query_mass) * 1e6);
            String found_formula(mass_mappings_[hit_idx[i]].formula);

            AccurateMassSearchResult ams_result;
            ams_result.setAdductMass(adduct_mass);
            ams_result.setQueryMass(query_mass);
            ams_result.setFoundMass(found_mass);
            ams_result.setCharge(adduct_charge);
            ams_result.setErrorPPM(found_error_ppm);
            ams_result.setMatchingIndex(hit_idx[i]);
            ams_result.setFoundAdduct(adduct_name);
            ams_result.setEmpiricalFormula(found_formula);
            ams_result.setMatchingHMDBids(mass_mappings_[hit_idx[i]].massIDs);

            results.push_back(ams_result);

            // ams_result.outputResults();
            // std::cout << "****************************************************" << std::endl;
        }

    }

    return;
}

void AccurateMassSearchEngine::queryByFeature(const Feature& feature, const Size& feature_index, std::vector<AccurateMassSearchResult>& results)
{
    if (!is_initialized_) init_(); // parse DB

    std::vector<AccurateMassSearchResult> results_part;

    queryByMass(feature.getMZ(), feature.getCharge(), results_part);

    for (Size hit_idx = 0; hit_idx < results_part.size(); ++hit_idx)
    {
        results_part[hit_idx].setObservedRT(feature.getRT());
        results_part[hit_idx].setSourceFeatureIndex(feature_index);
        results_part[hit_idx].setObservedIntensity(feature.getIntensity());
    }

    std::copy(results_part.begin(), results_part.end(), std::back_inserter(results));
}

void AccurateMassSearchEngine::queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, std::vector<AccurateMassSearchResult>& results)
{
    if (!is_initialized_) init_(); // parse DB

    std::vector<AccurateMassSearchResult> results_part;

    queryByMass(cfeat.getMZ(), cfeat.getCharge(), results_part);

    ConsensusFeature::HandleSetType ind_feats(cfeat.getFeatures());


    //    for ( ; f_it != ind_feats.end(); ++f_it)
    //    {
    //        std::cout << f_it->getRT() << "\t" << f_it->getMZ() << "\t" << f_it->getIntensity() << std::endl;
    //    }

    ConsensusFeature::const_iterator f_it = ind_feats.begin();
    std::vector<DoubleReal> tmp_f_ints;
    for (Size map_idx = 0; map_idx < number_of_maps; ++map_idx)
    {
        // std::cout << "map idx: " << f_it->getMapIndex() << std::endl;
        if (map_idx == f_it->getMapIndex())
        {
            tmp_f_ints.push_back(f_it->getIntensity());
            ++f_it;
        }
        else
        {
            tmp_f_ints.push_back(0.0);
        }
    }


    for (Size hit_idx = 0; hit_idx < results_part.size(); ++hit_idx)
    {
        results_part[hit_idx].setObservedRT(cfeat.getRT());
        results_part[hit_idx].setSourceFeatureIndex(cf_index);
        // results_part[hit_idx].setObservedIntensity(cfeat.getIntensity());
        results_part[hit_idx].setIndividualIntensities(tmp_f_ints);
    }

    std::copy(results_part.begin(), results_part.end(), std::back_inserter(results));
}


void AccurateMassSearchEngine::init_()
{
     // Loads the default mapping file (chemical formulas -> HMDB IDs)
    parseMappingFile_(db_mapping_file_);
    // This loads additional properties like common name, smiles, and inchi key for each HMDB id
    parseStructMappingFile_(db_struct_file_);

    parseAdductsFile_(pos_adducts_fname_, pos_adducts_);
    parseAdductsFile_(neg_adducts_fname_, neg_adducts_);
    
    is_initialized_ = true;
}

const String& AccurateMassSearchEngine::getInternalIonMode()
{
  return ion_mode_internal_;
}

void AccurateMassSearchEngine::run(const FeatureMap<>& fmap, MzTab& mztab_out)
{
    if (!is_initialized_) init_(); // parse DB

    if (ion_mode_ == "auto")
    {
       resolveAutoMode_(fmap);
    }
    else
    { // just copy
      ion_mode_internal_ = ion_mode_;
    }
    
    // map for storing overall results
    QueryResultsTable overall_results;

    for (Size i = 0; i < fmap.size(); ++i)
    {
        // Feature().getMetaValue(3)
        std::vector<AccurateMassSearchResult> query_results;

        // std::cout << i << ": " << fmap[i].getMetaValue(3) << " mass: " << fmap[i].getMZ() << " num_traces: " << fmap[i].getMetaValue("num_of_masstraces") << " charge: " << fmap[i].getCharge() << std::endl;
        queryByFeature(fmap[i], i, query_results);

        if (query_results.size() == 0) continue;

        if (iso_similarity_ && (Size)fmap[i].getMetaValue("num_of_masstraces") > 1)
        {
            // compute isotope pattern similarities and determine best matching one
            DoubleReal best_iso_sim(std::numeric_limits<DoubleReal>::max());
            Size best_iso_idx(0);

            for (Size hit_idx = 0; hit_idx < query_results.size(); ++hit_idx)
            {
                String emp_formula(query_results[hit_idx].getFormulaString());
                DoubleReal iso_sim(computeIsotopePatternSimilarity_(fmap[i], emp_formula));
                query_results[hit_idx].setIsotopesSimScore(iso_sim);

                if (iso_sim > best_iso_sim)
                {
                    best_iso_sim = iso_sim;
                    best_iso_idx = hit_idx;
                }
            }

            std::vector<AccurateMassSearchResult> tmp_results;
            tmp_results.push_back(query_results[best_iso_idx]);

            // keep the best AccurateMassSearchResult, drop all other hits
            query_results = tmp_results;
        }

        // debug output
        //        for (Size hit_idx = 0; hit_idx < query_results.size(); ++hit_idx)
        //        {
        //            query_results[hit_idx].outputResults();
        //        }

        // String feat_label(fmap[i].getMetaValue(3));
        overall_results.push_back(query_results);
    }

    LOG_INFO << "Found "<< overall_results.size() << " matched masses (with at least one hit each) from " << fmap.size() << " features." << std::endl;

    exportMzTab_(overall_results, mztab_out);

    return;
}

void AccurateMassSearchEngine::run(const ConsensusMap& cmap, MzTab& mztab_out)
{
    if (!is_initialized_) init_(); // parse DB
    
    if (ion_mode_ == "auto")
    {
       resolveAutoMode_(cmap);
    }
    else
    { // just copy
      ion_mode_internal_ = ion_mode_;
    }
    ConsensusMap::FileDescriptions fd_map = cmap.getFileDescriptions();
    Size num_of_maps = fd_map.size();

    // map for storing overall results
    QueryResultsTable overall_results;

    for (Size i = 0; i < cmap.size(); ++i)
    {
        std::vector<AccurateMassSearchResult> query_results;
        // std::cout << i << ": " << cmap[i].getMetaValue(3) << " mass: " << cmap[i].getMZ() << " num_traces: " << cmap[i].getMetaValue("num_of_masstraces") << " charge: " << cmap[i].getCharge() << std::endl;
        queryByConsensusFeature(cmap[i], i, num_of_maps, query_results);
        overall_results.push_back(query_results);
    }

    exportMzTab_(overall_results, mztab_out);
    return;
}



void AccurateMassSearchEngine::exportMzTab_(const QueryResultsTable& overall_results, MzTab& mztab_out)
{
    // iterate the overall results table

    String unit_id("AccMassSearch");
    MzTabSmallMoleculeSectionData sm_data_section;
    MzTabSmallMoleculeSectionRows all_sm_rows;

    Size id_group(1);

    std::map<String, UInt> adduct_stats;                    // adduct --> # occurences
    std::map<String, std::set<Size> > adduct_stats_unique;   // adduct --> # occurences (count for each feature only once)

    for (QueryResultsTable::const_iterator tab_it = overall_results.begin(); tab_it != overall_results.end(); ++tab_it)
    {
        // std::cout << tab_it->first << std::endl;

        for (Size hit_idx = 0; hit_idx < tab_it->size(); ++hit_idx)
        {
            // tab_it->second[hit_idx].outputResults();

            std::vector<String> matching_ids = tab_it->at(hit_idx).getMatchingHMDBids();

            // iterate over multiple IDs, generate a new row for each one

            for (Size id_idx = 0; id_idx < matching_ids.size(); ++id_idx)
            {
                MzTabSmallMoleculeSectionRow mztab_row_record;

                // set the identifier field
                String hid_temp = matching_ids[id_idx];
                MzTabString hmdb_id;
                hmdb_id.set(hid_temp);
                std::vector<MzTabString> hmdb_id_dummy;
                hmdb_id_dummy.push_back(hmdb_id);
                MzTabStringList string_dummy_list;
                string_dummy_list.set(hmdb_id_dummy);

                mztab_row_record.identifier = string_dummy_list;

                // set the chemical formula field
                MzTabString chem_form;
                String form_temp = tab_it->at(hit_idx).getFormulaString();
                chem_form.set(form_temp);

                mztab_row_record.chemical_formula = chem_form;

                // set the smiles field
                String smi_temp = hmdb_properties_mapping_[hid_temp][1];       // extract SMILES from struct mapping file
                MzTabString smi_string;
                smi_string.set(smi_temp);

                mztab_row_record.smiles = smi_string;

                // set the inchi_key field
                String inchi_temp = hmdb_properties_mapping_[hid_temp][2];       // extract INCHIKEY from struct mapping file
                MzTabString inchi_key;
                inchi_key.set(inchi_temp);

                mztab_row_record.inchi_key = inchi_key;

                // set description field (we use it for the common name of the compound)
                String name_temp = hmdb_properties_mapping_[hid_temp][0];
                MzTabString common_name;
                common_name.set(name_temp);

                mztab_row_record.description = common_name;


                // set mass_to_charge field
                DoubleReal mz_temp = tab_it->at(hit_idx).getAdductMass();
                MzTabDouble mass_to_charge;
                mass_to_charge.set(mz_temp);

                mztab_row_record.mass_to_charge = mass_to_charge;


                // set charge field
                Int ch_temp = tab_it->at(hit_idx).getCharge();
                MzTabDouble mcharge;
                mcharge.set(ch_temp);

                mztab_row_record.charge = mcharge;


                // set RT field
                DoubleReal rt_temp(tab_it->at(hit_idx).getObservedRT());
                MzTabDouble rt_temp2;
                rt_temp2.set(rt_temp);
                std::vector<MzTabDouble> rt_temp3;
                rt_temp3.push_back(rt_temp2);
                MzTabDoubleList observed_rt;
                observed_rt.set(rt_temp3);

                mztab_row_record.retention_time = observed_rt;


                // set database field
                String dbname_temp = "HMDB";
                MzTabString dbname;
                dbname.set(dbname_temp);

                mztab_row_record.database = dbname;


                // set database_version field
                String dbver_temp = "3.5";
                MzTabString dbversion;
                dbversion.set(dbver_temp);

                mztab_row_record.database_version = dbversion;


                // set smallmolecule_abundance_sub
                // check if we deal with a feature or consensus feature

                std::vector<DoubleReal> indiv_ints(tab_it->at(hit_idx).getIndividualIntensities());
                std::vector<MzTabDouble> int_temp3;

                if (indiv_ints.size() == 0)
                {
                    DoubleReal int_temp(tab_it->at(hit_idx).getObservedIntensity());
                    MzTabDouble int_temp2;
                    int_temp2.set(int_temp);
                    int_temp3.push_back(int_temp2);
                }
                else
                {
                    for (Size ii = 0; ii < indiv_ints.size(); ++ii)
                    {
                        DoubleReal int_temp(indiv_ints[ii]);
                        MzTabDouble int_temp2;
                        int_temp2.set(int_temp);
                        int_temp3.push_back(int_temp2);
                    }
                }

                mztab_row_record.smallmolecule_abundance_sub = int_temp3;


                // set smallmolecule_abundance_stdev_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
                DoubleReal stdev_temp(0.0);
                MzTabDouble stdev_temp2;
                stdev_temp2.set(stdev_temp);
                std::vector<MzTabDouble> stdev_temp3;

                if (indiv_ints.size() == 0)
                {
                    stdev_temp3.push_back(stdev_temp2);
                }
                else
                {
                    for (Size ii = 0; ii < indiv_ints.size(); ++ii)
                    {
                        stdev_temp3.push_back(stdev_temp2);
                    }
                }

                mztab_row_record.smallmolecule_abundance_stdev_sub = stdev_temp3;


                // set smallmolecule_abundance_std_error_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
                DoubleReal stderr_temp(0.0);
                MzTabDouble stderr_temp2;
                stderr_temp2.set(stderr_temp);
                std::vector<MzTabDouble> stderr_temp3;

                if (indiv_ints.size() == 0)
                {
                    stderr_temp3.push_back(stderr_temp2);
                }
                else
                {
                    for (Size ii = 0; ii < indiv_ints.size(); ++ii)
                    {
                        stderr_temp3.push_back(stderr_temp2);
                    }
                }


                mztab_row_record.smallmolecule_abundance_std_error_sub = stderr_temp3;


                // optional columns:
                std::vector<MzTabOptionalColumnEntry> optionals;

                // ppm error
                MzTabString ppmerr;
                ppmerr.set(String(tab_it->at(hit_idx).getErrorPPM()));
                MzTabOptionalColumnEntry col0;
                col0.first = "opt_ppm_error";
                col0.second = ppmerr;
                optionals.push_back(col0);

                // set found adduct ion
                String addion_temp(tab_it->at(hit_idx).getFoundAdduct());
                MzTabString addion;
                addion.set(addion_temp);
                MzTabOptionalColumnEntry col1;
                col1.first = "opt_adduct_ion";
                col1.second = addion;
                optionals.push_back(col1);
                ++adduct_stats[addion_temp]; // just some stats


                // set isotope similarity score
                DoubleReal sim_score_temp(tab_it->at(hit_idx).getIsotopesSimScore());
                std::stringstream read_in;
                read_in << sim_score_temp;
                String sim_score_temp2(read_in.str());
                MzTabString sim_score;
                sim_score.set(sim_score_temp2);
                MzTabOptionalColumnEntry col2;
                col2.first = "opt_isosim_score";
                col2.second = sim_score;
                optionals.push_back(col2);


                // set id group; rows with the same id group number originated from the same feature
                adduct_stats_unique[addion_temp].insert(id_group); // stats ...
                String id_group_temp(id_group);
                MzTabString id_group_str;
                id_group_str.set(id_group_temp);
                MzTabOptionalColumnEntry col3;
                col3.first = "opt_id_group";
                col3.second = id_group_str;
                optionals.push_back(col3);
                mztab_row_record.opt_ = optionals;
                all_sm_rows.push_back(mztab_row_record);
            }
        }
        ++id_group;
    }

    sm_data_section[unit_id] = all_sm_rows;
    mztab_out.setSmallMoleculeSectionData(sm_data_section);


    // print some adduct stats:
    LOG_INFO << "Adduct stats as 'adduct: #peaks explained (#total db entries)'\n";
    for (std::map<String, UInt>::const_iterator it = adduct_stats.begin(); it!= adduct_stats.end(); ++it)
    {
      LOG_INFO << "  " << it->first << ": " << adduct_stats_unique[it->first].size() << " (" << it->second << ")\n";
    }
    LOG_INFO << std::endl;

}


/// protected methods

void AccurateMassSearchEngine::updateMembers_()
{
    mass_error_value_ = (DoubleReal)param_.getValue("mass_error_value");
    mass_error_unit_ = (String)param_.getValue("mass_error_unit");
    ion_mode_ = (String)param_.getValue("ionization_mode");
    ion_mode_internal_ = ion_mode_; // just copy, since we have not seen any data yet
    
    iso_similarity_ = param_.getValue("isotopic_similarity").toBool();
    
    // use defaults if empty for all .tsv files
    db_mapping_file_ = (String)param_.getValue("db:mapping");
    if (db_mapping_file_.trim().empty()) db_mapping_file_ = (String)defaults_.getValue("db:mapping");
    db_struct_file_ = (String)param_.getValue("db:struct");
    if (db_struct_file_.trim().empty()) db_struct_file_ = (String)defaults_.getValue("db:struct");

    pos_adducts_fname_ = (String)param_.getValue("positive_adducts_file");
    if (pos_adducts_fname_.trim().empty()) pos_adducts_fname_ = (String)defaults_.getValue("positive_adducts_file");
    neg_adducts_fname_ = (String)param_.getValue("negative_adducts_file");
    if (neg_adducts_fname_.trim().empty()) neg_adducts_fname_ = (String)defaults_.getValue("negative_adducts_file");

    // database names might have changed, so parse files again before next query
    is_initialized_ = false;
}

/// private methods

void AccurateMassSearchEngine::parseMappingFile_(const String& db_mapping_file)
{
    mass_mappings_.clear();

    // load map_fname mapping file
    String filename = db_mapping_file;

    // load map_fname mapping file
    if (!File::readable(filename))
    {
      // throws Exception::FileNotFound if not found
      filename = File::find(filename);
    }

    String line;
    std::stringstream str_buf;
    std::istream_iterator<String> eol;

    // LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

    std::ifstream ifs(filename.c_str());
    while (getline(ifs, line))
    {
        str_buf.clear();
        str_buf << line;
        // std::cout << line << std::endl;
        std::istream_iterator<String> istr_it(str_buf);

        Size word_count(0);
        MappingEntry_ entry;

        while (istr_it != eol)
        {
            // LOG_DEBUG << *istr_it << " ";
            if (word_count == 0)
            {
                entry.mass = istr_it->toDouble();
            }
            else if (word_count == 1)
            {
                entry.formula = *istr_it;
                if (entry.mass == 0)
                {  // recompute mass from formula
                   entry.mass = EmpiricalFormula(entry.formula).getMonoWeight();
                   //std::cerr << "mass of " << entry.formula << " is " << entry.mass << "\n";
                }
            }
            else
            {
                entry.massIDs.push_back(*istr_it);
            }

            ++word_count;
            ++istr_it;
        }
        // LOG_DEBUG << std::endl;

        if (entry.massIDs.size() > 0)
        {
            mass_mappings_.push_back(entry);
        }
    }

    std::sort(mass_mappings_.begin(), mass_mappings_.end(), CompareEntryAndMass_());

    LOG_INFO << "Read " << mass_mappings_.size() << " entries from mapping file!" << std::endl;

    return;
}

void AccurateMassSearchEngine::parseStructMappingFile_(const String& db_struct_file)
{
    hmdb_properties_mapping_.clear();

    String filename = db_struct_file;

    // load map_fname mapping file
    if (!File::readable(filename))
    {
      // throws Exception::FileNotFound if not found
      filename = File::find(filename);
    }

    std::ifstream ifs(filename.c_str());
    String line;
    // LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

    std::vector<String> parts;
    while (getline(ifs, line))
    {
        line.split("\t", parts);

        if (parts.size() > 1)
        {

            String hmdb_id_key(parts[0]);
            std::vector<String> props;

            std::copy(parts.begin() + 1, parts.end(), std::back_inserter(props));
            // LOG_DEBUG << std::endl;

            if (props.size() == 3)
            {
                hmdb_properties_mapping_[hmdb_id_key] = props;
            }
            else
            {
                LOG_WARN << "Properties incomplete for " << hmdb_id_key << std::endl;
                for (Size i = 0; i < props.size(); ++i)
                {
                    LOG_WARN << props[i] << std::endl;
                }
            }
        }
    }

    return;
}

void AccurateMassSearchEngine::parseAdductsFile_(const String& filename, StringList& result)
{
  result.clear();

  String fname = filename;
  // search for mapping file
  if (!File::readable(fname))
  {
    // throws Exception::FileNotFound if not found
    fname = File::find(filename);
  }

  std::ifstream ifs(fname.c_str());
  String line;

  // LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

  while (getline(ifs, line))
  {
    line = line.trim();

    if (line != "")
    {
      result.push_back(line);
    }
  }

  LOG_INFO << "Read " << result.size() << " entries from adduct file '" << fname << "'." << std::endl;

  return;
}

void AccurateMassSearchEngine::searchMass_(const DoubleReal& neutral_query_mass, std::vector<Size>& hit_indices)
{
    DoubleReal diff_mz(0.0);
    // check if mass error window is given in ppm or Da
    if (mass_error_unit_ == "ppm")
    {
        diff_mz = (neutral_query_mass / 1e6) * mass_error_value_;
    }
    else
    {
        diff_mz = mass_error_value_;
    }

    //LOG_INFO << "searchMass: neutral_query_mass=" << neutral_query_mass << " diff_mz=" << diff_mz << " ppm allowed:" << mass_error_value_ << std::endl;


    // binary search for formulas which are within diff_mz distance
    if (mass_mappings_.size() < 1)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There are no entries found in mass-to-ids mapping file! Aborting... ", "0");
    }

    std::vector<MappingEntry_>::iterator lower_it = std::lower_bound(mass_mappings_.begin(), mass_mappings_.end(), neutral_query_mass - diff_mz, CompareEntryAndMass_());  // first element equal or larger
    std::vector<MappingEntry_>::iterator upper_it = std::upper_bound(mass_mappings_.begin(), mass_mappings_.end(), neutral_query_mass + diff_mz, CompareEntryAndMass_());  // first element greater than

    //std::cout << *lower_it << " " << *upper_it << "idx: " << lower_it - masskey_table_.begin() << " " << upper_it - masskey_table_.begin() << std::endl;
    Size start_idx = std::distance(mass_mappings_.begin(), lower_it);
    Size end_idx = std::distance(mass_mappings_.begin(), upper_it);

    hit_indices.clear();

    for (Size hit_idx = start_idx; hit_idx < end_idx; ++hit_idx)
    {
        hit_indices.push_back(hit_idx);
        //DoubleReal found_mass(mass_mappings_[hit_idx].mass);
        //DoubleReal found_error_ppm(((neutral_query_mass - found_mass)/neutral_query_mass)*1e6);
        // debug output
        //std::cout << std::setprecision(10) << "found mass: " << found_mass  << " with error: " << found_error_ppm << std::endl;
    }

    return;
}

void AccurateMassSearchEngine::parseAdductString_(const String& addstr, std::vector<String>& components)
{

}

void AccurateMassSearchEngine::computeNeutralMassFromAdduct_(const DoubleReal& adduct_mass, const String& adduct_string, DoubleReal& neutral_mass, Int& charge_value)
{
    // retrieve adduct and charge
    std::vector<String> tmpvec, tmpvec1, tmpvec2;
    String cp_str(adduct_string);
    cp_str.removeWhitespaces();

    bool is_intrinsic = false;

    cp_str.split(";", tmpvec);


    String molform = "", charge_str = "";
    // split term from adduct table into formula and charge, e.g. "M-H" and "1-"
    if (tmpvec.size() == 2)
    {
        // std::cout << "main: " << tmpvec[0] << " ch: " << tmpvec[1] << std::endl;
        molform = tmpvec[0].trim();

        if (molform == "M")
        {
            // std::cout << "intrinsic charge detected! " << adduct_string << " " << adduct_mass << " ";
            is_intrinsic = true;
        }
        charge_str = tmpvec[1].trim();
    }
    else
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not detect molecular ion or charge... maybe the semicolon missing?", cp_str);
    }

    // check if charge string is formatted correctly
    if ((!charge_str.hasSuffix("+")) && (!charge_str.hasSuffix("-")))
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Charge sign +/- in the end of the string is missing! ", charge_str);
    }

    // get charge and sign
    String charge_value_str(charge_str.substr(0, charge_str.size() - 1));
    charge_value = charge_value_str.toInt();
    String sign_char(charge_str.suffix(1));

    //  std::cout << "sign: " << sign_char << " value: " << charge_value << std::endl;

    if (sign_char=="+")
    {
        if (charge_value < 0)
        {
            charge_value *= -1;
        }
    }
    else if (sign_char=="-")
    {
        if (charge_value > 0)
        {
            charge_value *= -1;
        }
    }

    //    std::cout << "charge value: " << charge_value << std::endl;

    tmpvec.clear();

    // split by +
    molform.split("+", tmpvec);

    if (tmpvec.size() >= 2)
    {
        String left_hand(tmpvec[0].removeWhitespaces());

        if (left_hand == "")
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of + operator is missing! Aborting... Offending operator number 1", molform);
        }

        tmpvec1.push_back(left_hand);

        for (Size i = 1; i < tmpvec.size(); ++i)
        {
            String right_hand(tmpvec[i].removeWhitespaces());

            if (right_hand == "")
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of + operator is missing! Aborting... Offending operator number " + String(i + 1), molform);
            }

            tmpvec1.push_back("+");
            tmpvec1.push_back(right_hand);
        }
    }
    else // original formula contains no + sign
    {
        tmpvec1.push_back(molform);
    }

    for (Size i = 0; i < tmpvec1.size(); ++i)
    {
        std::vector<String> splits, newstr;
        tmpvec1[i].split("-", splits);

        if (splits.size() >= 2)
        {
            String left_hand(splits[0].removeWhitespaces());

            if (left_hand == "")
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of - operator is missing! Aborting... Offending operator number 1", molform);
            }

            newstr.push_back(left_hand);

            for (Size j = 1; j < splits.size(); j++)
            {
                String right_hand(splits[j].removeWhitespaces());

                if (right_hand == "")
                {
                    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of + operator is missing! Aborting... Offending operator number " + String(i + 1), molform);
                }

                newstr.push_back("-");
                newstr.push_back(right_hand);
            }

            tmpvec2.insert(tmpvec2.end(), newstr.begin(), newstr.end());
        }
        else
        {
            tmpvec2.push_back(tmpvec1[i]);
        }
    }

    // some further sanity check if adduct formula is correct
    String m_part(tmpvec2[0]);
    // std::cout << m_part.at(m_part.size() - 1) << std::endl;

    if (m_part.compare(m_part.size() - 1, 1, "M") != 0)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "First term of adduct string must contain the molecular entity M: ", m_part);
    }

    String mult_str(m_part.substr(0, m_part.size() - 1));

    DoubleReal mol_multiplier(1.0);

    if (mult_str != "")
    {
        mol_multiplier = mult_str.toDouble();
    }
    // std::cout << m_part << " " << mol_multiplier << std::endl;


    // time to evaluate the adduct string and compute the neutral (query) mass

    // first decharge the observed (=adducted) mass...
    neutral_mass = std::abs(charge_value) * adduct_mass;
    String last_op("");

    // add/subtract each adduct compound...
    for (Size part_idx = 1; part_idx < tmpvec2.size(); ++part_idx)
    {
        if (tmpvec2[part_idx] == "+")
        {
            last_op = "+";
            continue;
        }
        else if (tmpvec2[part_idx] == "-")
        {
            last_op = "-";
            continue;
        }

        // std::cout << "putting " << tmpvec2[part_idx] << " into a formula with mass ";

        // check if formula has got a stoichometry factor in front
        String formula_str(tmpvec2[part_idx]);
        const char first_char = formula_str[0];

        DoubleReal stoichio_factor(1.0);

        if (isdigit(first_char))
        {
            String tmp_factor(first_char);
            stoichio_factor = tmp_factor.toDouble();
            formula_str = formula_str.substr(1, formula_str.size());
        }

        // std::cout << stoichio_factor << "*" << formula_str << " ";
        EmpiricalFormula part_formula(formula_str);
        // std::cout << part_formula.getMonoWeight() << std::endl;

        if (last_op == "+")
        {
            neutral_mass -= stoichio_factor * part_formula.getMonoWeight();
            last_op = "";
        }
        else if (last_op == "-")
        {
            neutral_mass += stoichio_factor * part_formula.getMonoWeight();
            last_op = "";
        }
    }


    // correct for electron masses
    DoubleReal electrons_mass_diff(charge_value * Constants::ELECTRON_MASS_U);

    // std::cout << "electron mass: " << Constants::ELECTRON_MASS_U << " " << Constants::ELECTRON_MASS << " " << electrons_mass_diff << std::endl;
    if (!is_intrinsic)
    {
        neutral_mass += electrons_mass_diff;
    }
    // divide by stoichiometry factor
    neutral_mass /= mol_multiplier;

    // std::cout << " neutral: " << neutral_mass << std::endl;

    return;
}

DoubleReal AccurateMassSearchEngine::computeCosineSim_(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y)
{
    if (x.size() != y.size())
    {
        return 0.0;
    }

    DoubleReal mixed_sum(0.0);
    DoubleReal x_squared_sum(0.0);
    DoubleReal y_squared_sum(0.0);


    for (Size i = 0; i < x.size(); ++i)
    {
        mixed_sum += x[i] * y[i];
        x_squared_sum += x[i] * x[i];
        y_squared_sum += y[i] * y[i];
    }

    DoubleReal denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));

    return (denom > 0.0) ? mixed_sum / denom : 0.0;
}

DoubleReal AccurateMassSearchEngine::computeEuclideanDist_(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y)
{
    if (x.size() != y.size())
    {
        return -100.0;
    }

    DoubleReal sum_of_squares(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
        sum_of_squares += (x[i] - y[i]) * (x[i] - y[i]);
    }

    return std::sqrt(sum_of_squares);
}

DoubleReal AccurateMassSearchEngine::computeIsotopePatternSimilarity_(const Feature& feat, const EmpiricalFormula& form)
{
    Size num_traces = (Size)feat.getMetaValue("num_of_masstraces");
    Size MAX_THEORET_ISOS(5);

    Size common_size = num_traces < MAX_THEORET_ISOS ? num_traces : MAX_THEORET_ISOS;

    IsotopeDistribution iso_dist(form.getIsotopeDistribution(MAX_THEORET_ISOS));

    DoubleReal max_iso_prob(iso_dist.begin()->second);

    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin(); iso_it != (iso_dist.begin() + common_size); ++iso_it)
    {
        // std::cout << "first: " << iso_it->first << " second: " << iso_it->second << std::endl;
        if (iso_it->second > max_iso_prob)
        {
            max_iso_prob = iso_it->second;
        }
    }

    std::vector<DoubleReal> normed_iso_ratios;

    // std::cout << "theoret. iso: ";
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin(); iso_it != (iso_dist.begin() + common_size); ++iso_it)
    {
        DoubleReal temp_ratio((iso_it->second) / max_iso_prob);
        normed_iso_ratios.push_back(temp_ratio);

        // std::cout << temp_ratio << " ";
    }

    // std::cout << "\nact. iso: ";

    DoubleReal max_feat_int((DoubleReal)feat.getMetaValue("masstrace_intensity_0"));

    std::vector<DoubleReal> normed_feat_ratios;

    for (Size int_idx = 0; int_idx < common_size; ++int_idx)
    {
        std::stringstream read_in;
        read_in << int_idx;
        String identifier(read_in.str());
        DoubleReal mt_int = (DoubleReal)feat.getMetaValue("masstrace_intensity_" + identifier);

        normed_feat_ratios.push_back(mt_int);

        if (mt_int > max_feat_int)
        {
            max_feat_int = mt_int;
        }
    }

    // normalize with max isotope intensity
    for (Size int_idx = 0; int_idx < common_size; ++int_idx)
    {
        normed_feat_ratios[int_idx] /= max_feat_int;
    }

    return computeCosineSim_(normed_iso_ratios, normed_feat_ratios);
    // return computeEuclideanDist_(normed_iso_ratios, normed_feat_ratios);
}


} // closing namespace OpenMS
