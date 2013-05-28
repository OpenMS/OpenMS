// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Erhan Kenar $
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
AccurateMassSearchResult::AccurateMassSearchResult()
    :  adduct_mass_(),
      query_mass_(),
      found_mass_(),
      charge_(),
      error_ppm_(),
      matching_index_(),
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
AccurateMassSearchResult::AccurateMassSearchResult(const AccurateMassSearchResult& source)
    :  adduct_mass_(source.adduct_mass_),
      query_mass_(source.query_mass_),
      found_mass_(source.found_mass_),
      charge_(source.charge_),
      error_ppm_(source.error_ppm_),
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
    matching_index_ = rhs.matching_index_;
    found_adduct_ = rhs.found_adduct_;
    empirical_formula_ = rhs.empirical_formula_;
    matching_hmdb_ids_ = rhs.matching_hmdb_ids_;
    isotopes_sim_score_ = rhs.isotopes_sim_score_;

    return *this;
}

void AccurateMassSearchResult::outputResults() const
{
    std::cout << "adduct_mass: " << std::setprecision(8) << adduct_mass_ << std::endl;
    std::cout << "query_mass: " << query_mass_ << std::endl;
    std::cout << "found_mass: " << found_mass_ << std::endl;
    std::cout << "charge: " << charge_ << std::endl;
    std::cout << "error ppm: " << error_ppm_ << std::endl;
    std::cout << "matching idx: " << matching_index_ << std::endl;

    std::cout << "found_adduct_: " << found_adduct_ << std::endl;
    std::cout << "emp. formula: " << empirical_formula_ << std::endl;
    std::cout << "matching HMDB ids:";

    for (Size i = 0; i < matching_hmdb_ids_.size(); ++i)
    {
        std::cout << " " << matching_hmdb_ids_[i];
    }

    std::cout << std::endl;
    std::cout << "isocheck sim score: " << isotopes_sim_score_ << std::endl;
}


AccurateMassSearchEngine::AccurateMassSearchEngine() :
    DefaultParamHandler("AccurateMassSearchEngine"), ProgressLogger()
{
    defaults_.setValue("mass_error_value", 5.0, "Tolerance allowed for accurate mass search.");

    defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
    defaults_.setValidStrings("mass_error_unit", StringList::create(("ppm,Da")));

    defaults_.setValue("ionization_mode", "positive", "Positive or negative ionization mode?");
    defaults_.setValidStrings("ionization_mode", StringList::create(("positive,negative")));

    defaults_.setValue("isotopic_similarity", "true", "Computes a similarity score for each hit (only if the feature exhibits at least two isotopic mass traces).");
    defaults_.setValidStrings("isotopic_similarity", StringList::create(("false,true")));

    defaults_.setValue("report_mode", "all", "Results are reported in one of several modes: Either (all) matching hits, the (top3) scoring hits, or the (best) scoring hit.");
    defaults_.setValidStrings("report_mode", StringList::create(("all,top3,best")));

    std::vector<String> pos_adducts, neg_adducts;

    pos_adducts.push_back("M+3H;3+");
    pos_adducts.push_back("M+2H+Na;3+");
    pos_adducts.push_back("M+H+2Na;3+");
    pos_adducts.push_back("M+3Na;3+");
    pos_adducts.push_back("M+2H;2+");
    pos_adducts.push_back("M+H+NH4;2+");
    pos_adducts.push_back("M+H+Na;2+");
    pos_adducts.push_back("M+H+K;2+");
    pos_adducts.push_back("M+CH3CN+2H;2+"); // acetonitrile
    pos_adducts.push_back("M+2Na;2+");
    pos_adducts.push_back("M+2CH3CN+2H;2+");
    pos_adducts.push_back("M+3CH3CN+2H;2+");
    pos_adducts.push_back("M+H;1+");
    pos_adducts.push_back("M+NH4;1+");
    pos_adducts.push_back("M+Na;1+");
    pos_adducts.push_back("M+CH3OH+H;1+");
    pos_adducts.push_back("M+K;1+");
    pos_adducts.push_back("M+CH3CN+H;1+");
    pos_adducts.push_back("M+2Na-H;1+");
    pos_adducts.push_back("M+C3H8O;1+"); // isopropanol
    pos_adducts.push_back("M+CH3CN+Na;1+");
    pos_adducts.push_back("M+2K-H;1+");
    pos_adducts.push_back("M+C2H6OS;1+"); // DMSO; dimethylsulfoxide
    pos_adducts.push_back("M+2CH3CN+H;1+");
    pos_adducts.push_back("M+C3H8O+Na+H;1+");
    pos_adducts.push_back("2M+H;1+");
    pos_adducts.push_back("2M+NH4;1+");
    pos_adducts.push_back("2M+Na;1+");
    pos_adducts.push_back("2M+3H2O+2H;2+");
    pos_adducts.push_back("2M+K;1+");
    pos_adducts.push_back("2M+CH3CN+H;1+");
    pos_adducts.push_back("2M+CH3CN+Na;1+");

    StringList pos_adduct_list(pos_adducts);

    defaults_.setValue("positive_adducts", pos_adduct_list, "This is the list of potential positive adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts.", StringList::create("advanced"));

    StringList neg_adduct_list(neg_adducts);

    defaults_.setValue("negative_adducts", neg_adduct_list, "This is the list of potential negative adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts.", StringList::create("advanced"));


    defaultsToParam_();

    this->setLogType(CMD);

    // The default constructor loads the default mapping file (chemical formulas -> HMDB IDs)
    parseMappingFile_("");

    // This loads additional properties like common name, smiles, and inchi key for each HMDB id
    parseStructMappingFile_("");
}

AccurateMassSearchEngine::~AccurateMassSearchEngine()
{

}


/// public methods

void AccurateMassSearchEngine::queryByMass(const DoubleReal& adduct_mass, const DoubleReal& adduct_charge, std::vector<AccurateMassSearchResult>& results)
{
    // depending on ionization mode, look for positive oder negative adducts

    if (ion_mode_ == "positive")
    {
        for (Size adduct_idx = 0; adduct_idx < pos_adducts_.size(); ++adduct_idx)
        {
            std::vector<Size> hit_idx;
            DoubleReal query_mass, charge;
            String pos_adduct_name(pos_adducts_[adduct_idx]);

            computeNeutralMassFromAdduct_(adduct_mass, pos_adduct_name, query_mass, charge);

            // std::cout << "looking for " << pos_adducts_[adduct_idx] << std::endl;
            if ((adduct_charge > 0) && (charge != adduct_charge))
            {
                continue;
            }

            // get potential hits as indices in masskey_table
            searchMass_(query_mass, hit_idx);

            // store information from query hits in AccurateMassSearchResult objects
            for (Size i = 0; i < hit_idx.size(); ++i)
            {
                DoubleReal found_mass(*(masskey_table_.begin() + hit_idx[i]));
                DoubleReal found_error_ppm(((query_mass - found_mass)/query_mass)*1000000);
                String found_formula(mass_formula_mapping_[hit_idx[i]]);

                AccurateMassSearchResult ams_result;
                ams_result.setAdductMass(adduct_mass);
                ams_result.setQueryMass(query_mass);
                ams_result.setFoundMass(found_mass);
                ams_result.setCharge(adduct_charge);
                ams_result.setErrorPPM(found_error_ppm);
                ams_result.setMatchingIndex(hit_idx[i]);
                ams_result.setFoundAdduct(pos_adduct_name);
                ams_result.setEmpiricalFormula(found_formula);

                std::vector<String> matching_hmdb_ids;

                for (Size j = 0; j < mass_id_mapping_[hit_idx[i]].size(); ++j)
                {
                    matching_hmdb_ids.push_back(mass_id_mapping_[hit_idx[i]][j]);
                }
                ams_result.setMatchingHMDBids(matching_hmdb_ids);

                results.push_back(ams_result);

                // ams_result.outputResults();
                // std::cout << "****************************************************" << std::endl;
            }

        }
    }
    else
    {
        // same for negative mode

    }
    return ;
}

void AccurateMassSearchEngine::queryByFeature(const Feature& feat, std::vector<AccurateMassSearchResult>& results)
{
    DoubleReal adduct_mass(feat.getMZ());
    DoubleReal adduct_charge(feat.getCharge());

    queryByMass(adduct_mass, adduct_charge, results);
}

void AccurateMassSearchEngine::run(const FeatureMap<> & fmap, MzTab& mztab_out)
{
    //    for (Size i = 0; i < mass_id_mapping_.size(); ++i)
    //    {
    //        std::cout << i << " : " << mass_formula_mapping_[i] << std::endl;
    //    }

    typedef std::map<String, std::vector<AccurateMassSearchResult> > QueryResultsTable;

    // map for storing overall results
    QueryResultsTable overall_results;

    for (Size i = 0; i < fmap.size(); ++i)
    {
        // Feature().getMetaValue(3)
        std::vector<AccurateMassSearchResult> query_results;

        // std::cout << fmap[i].getMetaValue(3) << " mass: " << fmap[i].getMZ() << " num_traces: " << fmap[i].getMetaValue("num_of_masstraces") << " charge: " << fmap[i].getCharge() << std::endl;
        queryByFeature(fmap[i], query_results);

        if (iso_similarity_ && (Size)fmap[i].getMetaValue("num_of_masstraces") > 1 && query_results.size() > 0)
        {
            // compute isotope pattern similarities and determine best matching one
            DoubleReal best_iso_sim(std::numeric_limits<DoubleReal>::max());
            Size best_iso_idx(0);

            for (Size hit_idx = 0; hit_idx < query_results.size(); ++hit_idx)
            {
                String emp_formula(query_results[hit_idx].getFormulaString());
                DoubleReal iso_sim(computeIsotopePatternSimilarity_(fmap[i], emp_formula));
                query_results[hit_idx].setIsotopesSimScore(iso_sim);

                if (iso_sim < best_iso_sim)
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

        String feat_label(fmap[i].getMetaValue(3));
        overall_results[feat_label] = query_results;
    }

    std::cout << "RESULTS: " << std::endl;
    // iterate the overall results table

    MzTabSmallMoleculeSectionData sm_data_section;
    MzTabSmallMoleculeSectionRows all_sm_rows;

    for (QueryResultsTable::const_iterator tab_it = overall_results.begin(); tab_it != overall_results.end(); ++tab_it)
    {
        std::cout << tab_it->first << std::endl;

        for (Size i = 0; i < tab_it->second.size(); ++i)
        {
            tab_it->second[i].outputResults();
            LOG_WARN << "outputted!" << std::endl;

            std::vector<String> matching_ids = tab_it->second[i].getMatchingHMDBids();

            // iterate over multiple IDs, generate a new row for each one

            for (Size i = 0; i < matching_ids.size(); ++i)
            {
                MzTabSmallMoleculeSectionRow mztab_row_record;

                // set the identifier
                MzTabString hmdb_id;
                hmdb_id.set(matching_ids[i]);
                std::vector<MzTabString> hmdb_id_dummy;
                hmdb_id_dummy.push_back(hmdb_id);
                MzTabStringList string_dummy_list;
                string_dummy_list.set(hmdb_id_dummy);

                LOG_WARN << "first " << std::endl;

                mztab_row_record.identifier = string_dummy_list;

                // set the chemical formula
                MzTabString chem_form;
                String str_temp = tab_it->second[i].getFormulaString();

                chem_form.set(str_temp);

                LOG_WARN << tab_it->second[i].getFormulaString() << std::endl;
//                mztab_row_record.chemical_formula = chem_form;
//                LOG_WARN << "second " << std::endl;



                //                MzTabStringList identifier; // The small molecule’s identifier.
                //                MzTabString chemical_formula; // Chemical formula of the identified compound.
                //                MzTabString smiles; // Molecular structure in SMILES format.
                //                MzTabString inchi_key; // InChi Key of the identified compound.
                //                MzTabString description; // Human readable description (i.e. the name)
                //                MzTabDouble mass_to_charge; // Precursor ion’s m/z.
                //                MzTabDouble charge; // Precursor ion’s charge.
                //                MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
                //                MzTabInteger taxid; // NEWT taxonomy for the species.
                //                MzTabString species; // Human readable name of the species
                //                MzTabString database; // Name of the used database.
                //                MzTabString database_version; // String Version of the database (and optionally # of compounds).
                //                MzTabInteger reliability; // (1-3) The identification reliability.
                //                MzTabString uri; // The source entry’s location.
                //                MzTabSpectraRef spectra_ref; // Spectra identifying the small molecule.
                //                MzTabParameterList search_engine; // Search engine(s) identifying the small molecule.
                //                MzTabParameterList search_engine_score; // Search engine(s) identifications score(s).
                //                MzTabModificationList modifications; // Modifications identified on the small molecule.
                //                std::vector<MzTabDouble> smallmolecule_abundance_sub; // Abundance in the subsample;
                //                std::vector<MzTabDouble> smallmolecule_abundance_stdev_sub; // Standard deviation of the abundance.
                //                std::vector<MzTabDouble> smallmolecule_abundance_std_error_sub; // Standard errpr of the abundance.
                //                std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.

                // all_sm_rows.push_back(mztab_row_record);

            }
        }
    }

}

/// protected methods

void AccurateMassSearchEngine::updateMembers_()
{
    mass_error_value_ = (DoubleReal)param_.getValue("mass_error_value");
    mass_error_unit_ = (String)param_.getValue("mass_error_unit");
    ion_mode_ = (String)param_.getValue("ionization_mode");
    iso_similarity_ = param_.getValue("isotopic_similarity").toBool();

    pos_adducts_ = (StringList)param_.getValue("positive_adducts");
    neg_adducts_ = (StringList)param_.getValue("negative_adducts");
}


/// private methods

void AccurateMassSearchEngine::parseMappingFile_(const String& map_fname)
{
    // load map_fname mapping file
    String fname;

    if (map_fname == "")
    {
        fname = File::find("CHEMISTRY/HMDBMappingFile.tsv");
    }
    else
    {
        fname = map_fname;
    }

    std::ifstream ifs(fname.c_str());

    String line;
    std::stringstream str_buf;
    std::istream_iterator<String> eol;

    // LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

    while (getline(ifs, line))
    {
        str_buf.clear();
        str_buf << line;
        std::istream_iterator<String> istr_it(str_buf);

        Size word_count(0);
        DoubleReal mass_key(-1.0);
        String formula_str;
        std::vector<String> hmdb_ids;

        while (istr_it != eol)
        {
            // LOG_DEBUG << *istr_it << " ";
            if (word_count == 0)
            {
                mass_key = istr_it->toDouble();
            }
            else if (word_count == 1)
            {
                formula_str = *istr_it;
            }
            else
            {
                hmdb_ids.push_back(*istr_it);
            }

            ++word_count;
            ++istr_it;
        }
        // LOG_DEBUG << std::endl;

        if (hmdb_ids.size() > 0)
        {
            masskey_table_.push_back(mass_key);
            mass_formula_mapping_.push_back(formula_str);
            mass_id_mapping_.push_back(hmdb_ids);
        }
    }

    if (masskey_table_.size() != mass_id_mapping_.size()
            && masskey_table_.size() != mass_formula_mapping_.size())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Parsing of mass-to-HMDB-IDs mapping failed... Sizes of masskey_table_ and mass_id_mapping_ differ!" + String(masskey_table_.size()), String(mass_id_mapping_.size()));
    }

    LOG_INFO << "masskey_table size: " << masskey_table_.size() << " mass_id_mapping size: " << mass_id_mapping_.size() << " mass_formula_mapping size: " << mass_formula_mapping_.size() << std::endl;

    return ;
}

void AccurateMassSearchEngine::parseStructMappingFile_(const String& map_fname)
{
    // load map_fname mapping file
    String fname;

    if (map_fname == "")
    {
        fname = File::find("CHEMISTRY/HMDB2StructMapping.tsv");
    }
    else
    {
        fname = map_fname;
    }

    std::ifstream ifs(fname.c_str());

    String line;

    // LOG_DEBUG << "parsing " << fname << " file..." << std::endl;

    while (getline(ifs, line))
    {
        std::vector<String> parts;
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

    return ;
}


void AccurateMassSearchEngine::searchMass_(const DoubleReal& neutral_query_mass, std::vector<Size>& hit_indices)
{
    DoubleReal diff_mz(0.0);
    // check if mass error window is given in ppm or Da
    if (mass_error_unit_ == "ppm")
    {
        diff_mz = (neutral_query_mass / 1000000) * mass_error_value_;
    }
    else
    {
        diff_mz = mass_error_value_;
    }

    // LOG_INFO << "searchMass: neutral_query_mass=" << neutral_query_mass << " diff_mz=" << diff_mz << std::endl;


    // binary search for formulas which are within diff_mz distance
    Size n_masskeys = masskey_table_.size();

    if (n_masskeys < 1)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There are no entries found in mass-to-ids mapping file! Aborting... ", String(n_masskeys));
    }

    std::vector<DoubleReal>::const_iterator lower_it = std::lower_bound(masskey_table_.begin(), masskey_table_.end(), neutral_query_mass - diff_mz);
    std::vector<DoubleReal>::const_iterator upper_it = std::upper_bound(masskey_table_.begin(), masskey_table_.end(), neutral_query_mass + diff_mz);

    //std::cout << *lower_it << " " << *upper_it << "idx: " << lower_it - masskey_table_.begin() << " " << upper_it - masskey_table_.begin() << std::endl;
    Size start_idx(lower_it - masskey_table_.begin());
    Size end_idx(upper_it - masskey_table_.begin());

    hit_indices.clear();

    for (Size hit_idx = start_idx; hit_idx < end_idx; ++hit_idx)
    {
        hit_indices.push_back(hit_idx);
        // DoubleReal found_mass(*(masskey_table_.begin() + hit_idx));
        // DoubleReal found_error_ppm(((neutral_query_mass - found_mass)/neutral_query_mass)*1000000);
        // debug output
        // std::cout << std::setprecision(10) << "found mass: " << found_mass  << " with error: " << found_error_ppm << std::endl;
    }

    return ;
}

void AccurateMassSearchEngine::parseAdductString_(const String& addstr, std::vector<String>& components)
{

}

void AccurateMassSearchEngine::computeNeutralMassFromAdduct_(const DoubleReal& adduct_mass, const String& adduct_string, DoubleReal& neutral_mass, DoubleReal& charge_value)
{
    // retrieve adduct and charge
    std::vector<String> tmpvec, tmpvec1, tmpvec2;
    String cp_str(adduct_string);
    cp_str.removeWhitespaces();

    cp_str.split(";", tmpvec);


    String molform = "", charge_str = "";

    if (tmpvec.size() == 2)
    {
        // std::cout << "main: " << tmpvec[0] << " ch: " << tmpvec[1] << std::endl;
        molform = tmpvec[0].trim();
        charge_str = tmpvec[1].trim();
    }
    else
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not detect molecular ion or charge... maybe the semicolon missing?", cp_str);
    }

    // check if charge string is formatted correctly
    if ((charge_str.compare(charge_str.size() - 1, 1, "+") != 0) && (charge_str.compare(charge_str.size() - 1, 1, "-") != 0))
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Charge sign +/- in the end of the string is missing! ", charge_str);
    }

    String charge_value_str(charge_str.substr(0, charge_str.size() - 1));

    charge_value = charge_value_str.toDouble();
    String sign_char(charge_str.substr(charge_str.size() - 1, 1));

    //  std::cout << "sign: " << sign_char << " value: " << charge_value << std::endl;

    if (sign_char.compare(0,1, "+") == 0)
    {
        if (charge_value < 0)
        {
            charge_value *= -1.0;
        }
    }
    else if (sign_char.compare(0,1, "-") == 0)
    {
        if (charge_value > 0)
        {
            charge_value *= -1.0;
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
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of + operator is missing! Aborting... Offending operator number " + String(i+1), molform);
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

        //        for (Size j = 0; j < splits.size(); ++j)
        //        {
        //            std::cout << "splits: " << splits[j] << ",";
        //        }


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
                    throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Left- or righthand side of + operator is missing! Aborting... Offending operator number " + String(i+1), molform);
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

    // first decharge the adduct mass...
    neutral_mass = std::fabs(charge_value) * adduct_mass;
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
        char first_char = formula_str[0];

        DoubleReal stoichio_factor(1.0);

        if (std::isdigit(first_char))
        {
            String tmp_factor(first_char);
            stoichio_factor = tmp_factor.toDouble();
            formula_str = formula_str.substr(1,formula_str.size());
        }

        // std::cout << stoichio_factor << "*" << formula_str << " ";
        EmpiricalFormula part_formula(formula_str);
        // std::cout << part_formula.getMonoWeight() << std::endl;

        if (last_op == "+")
        {
            neutral_mass -= stoichio_factor*part_formula.getMonoWeight();
            last_op = "";
        }
        else if (last_op == "-")
        {
            neutral_mass += stoichio_factor*part_formula.getMonoWeight();
            last_op = "";
        }
    }


    // correct for electron masses
    DoubleReal electrons_mass_diff(charge_value * Constants::ELECTRON_MASS_U);

    // std::cout << "electron mass: " << Constants::ELECTRON_MASS_U << " " << Constants::ELECTRON_MASS << " " << electrons_mass_diff << std::endl;
    neutral_mass += electrons_mass_diff;

    // divide by stoichiometry factor
    neutral_mass /= mol_multiplier;

    return ;
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
        DoubleReal temp_ratio((iso_it->second)/max_iso_prob);
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

    // return computeCosineSim_(normed_iso_ratios, normed_feat_ratios);
    return computeEuclideanDist_(normed_iso_ratios, normed_feat_ratios);
}



} // closing namespace OpenMS
