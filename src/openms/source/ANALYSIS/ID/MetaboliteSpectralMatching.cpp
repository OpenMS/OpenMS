// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>


#include <numeric>
#include <boost/math/special_functions/factorials.hpp>

#include <boost/dynamic_bitset.hpp>

#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

namespace OpenMS
{

SpectralMatch::SpectralMatch() :
  observed_precursor_mass_(),
  observed_precursor_rt_(),
  found_precursor_mass_(),
  found_precursor_charge_(),
  matching_score_(),
  observed_spectrum_idx_(),
  matching_spectrum_idx_(),
  primary_id_(),
  secondary_id_(),
  common_name_(),
  sum_formula_(),
  inchi_string_(),
  smiles_string_(),
  precursor_adduct_()
{

}

/// Default destructor
SpectralMatch::~SpectralMatch()
{

}

/// Copy constructor
SpectralMatch::SpectralMatch(const SpectralMatch& sm) :
  observed_precursor_mass_(sm.observed_precursor_mass_),
  observed_precursor_rt_(sm.observed_precursor_rt_),
  found_precursor_mass_(sm.found_precursor_mass_),
  found_precursor_charge_(sm.found_precursor_charge_),
  matching_score_(sm.matching_score_),
  observed_spectrum_idx_(sm.observed_spectrum_idx_),
  matching_spectrum_idx_(sm.matching_spectrum_idx_),
  primary_id_(sm.primary_id_),
  secondary_id_(sm.secondary_id_),
  common_name_(sm.common_name_),
  sum_formula_(sm.sum_formula_),
  inchi_string_(sm.inchi_string_),
  smiles_string_(sm.smiles_string_),
  precursor_adduct_(sm.precursor_adduct_)
{

}

/// Assignment operator
SpectralMatch& SpectralMatch::operator=(const SpectralMatch& rhs)
{
  if (this == &rhs)
    return *this;

  observed_precursor_mass_ = rhs.observed_precursor_mass_;
  observed_precursor_rt_ = rhs.observed_precursor_rt_;
  found_precursor_mass_ = rhs.found_precursor_mass_;
  found_precursor_charge_ = rhs.found_precursor_charge_;
  matching_score_ = rhs.matching_score_;
  observed_spectrum_idx_ = rhs.observed_spectrum_idx_;
  matching_spectrum_idx_ = rhs.matching_spectrum_idx_;
  primary_id_ = rhs.primary_id_;
  secondary_id_ = rhs.secondary_id_;
  common_name_ = rhs.common_name_;
  sum_formula_ = rhs.sum_formula_;
  inchi_string_ = rhs.inchi_string_;
  smiles_string_ = rhs.smiles_string_;
  precursor_adduct_ = rhs.precursor_adduct_;

  return *this;
}

double SpectralMatch::getObservedPrecursorMass() const
{
  return observed_precursor_mass_;
}

void SpectralMatch::setObservedPrecursorMass(const double& qmass)
{
  observed_precursor_mass_ = qmass;
}

double SpectralMatch::getObservedPrecursorRT() const
{
  return observed_precursor_rt_;
}

void SpectralMatch::setObservedPrecursorRT(const double& prt)
{
  observed_precursor_rt_ = prt;
}

double SpectralMatch::getFoundPrecursorMass() const
{
  return found_precursor_mass_;
}

void SpectralMatch::setFoundPrecursorMass(const double& fmass)
{
  found_precursor_mass_ = fmass;
}

Int SpectralMatch::getFoundPrecursorCharge() const
{
  return found_precursor_charge_;
}

void SpectralMatch::setFoundPrecursorCharge(const Int& pch)
{
  found_precursor_charge_ = pch;
}

double SpectralMatch::getMatchingScore() const
{
  return matching_score_;
}

void SpectralMatch::setMatchingScore(const double& mscore)
{
  matching_score_ = mscore;
}

Size SpectralMatch::getObservedSpectrumIndex() const
{
  return observed_spectrum_idx_;
}

void SpectralMatch::setObservedSpectrumIndex(const Size& obs_spec_idx)
{
  observed_spectrum_idx_ = obs_spec_idx;
}

Size SpectralMatch::getMatchingSpectrumIndex() const
{
  return matching_spectrum_idx_;
}

void SpectralMatch::setMatchingSpectrumIndex(const Size& match_spec_idx)
{
  matching_spectrum_idx_ = match_spec_idx;
}

String SpectralMatch::getPrimaryIdentifier() const
{
  return primary_id_;
}

void SpectralMatch::setPrimaryIdentifier(const String& pid)
{
  primary_id_ = pid;
}

String SpectralMatch::getSecondaryIdentifier() const
{
  return secondary_id_;
}

void SpectralMatch::setSecondaryIdentifier(const String& sid)
{
  secondary_id_ = sid;
}

String SpectralMatch::getCommonName() const
{
  return common_name_;
}

void SpectralMatch::setCommonName(const String& cname)
{
  common_name_ = cname;
}

String SpectralMatch::getSumFormula() const
{
  return sum_formula_;
}

void SpectralMatch::setSumFormula(const String& sf)
{
  sum_formula_ = sf;
}

String SpectralMatch::getInchiString() const
{
  return inchi_string_;
}

void SpectralMatch::setInchiString(const String& istr)
{
  inchi_string_ = istr;
}

String SpectralMatch::getSMILESString() const
{
  return smiles_string_;
}

void SpectralMatch::setSMILESString(const String& sstr)
{
  smiles_string_ = sstr;
}

String SpectralMatch::getPrecursorAdduct() const
{
  return precursor_adduct_;
}

void SpectralMatch::setPrecursorAdduct(const String& padd)
{
  precursor_adduct_ = padd;
}


MetaboliteSpectralMatching::MetaboliteSpectralMatching() :
  DefaultParamHandler("MetaboliteSpectralMatching"), ProgressLogger()
{
  defaults_.setValue("prec_mass_error_value", 100.0, "Error allowed for precursor ion mass.");
  defaults_.setValue("frag_mass_error_value", 500.0, "Error allowed for product ions.");

  defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
  defaults_.setValidStrings("mass_error_unit", ListUtils::create<String>(("ppm,Da")));

  defaults_.setValue("report_mode", "top3", "Which results shall be reported: the top-three scoring ones or the best scoring one?");
  defaults_.setValidStrings("report_mode", ListUtils::create<String>(("top3,best")));

  defaults_.setValue("ionization_mode", "positive", "Positive or negative ionization mode?");
  defaults_.setValidStrings("ionization_mode", ListUtils::create<String>(("positive,negative")));


  defaultsToParam_();

  this->setLogType(CMD);
}

MetaboliteSpectralMatching::~MetaboliteSpectralMatching()
{

}


/// public methods

double MetaboliteSpectralMatching::computeHyperScore(MSSpectrum exp_spectrum, MSSpectrum db_spectrum,
                             const double& fragment_mass_error, const double& mz_lower_bound)
{

  double dot_product(0.0);
  Size matched_ions_count(0);

  // scan for matching peaks between observed and DB stored spectra
  for (MSSpectrum::iterator frag_it = exp_spectrum.MZBegin(mz_lower_bound); frag_it != exp_spectrum.end(); ++frag_it)
  {
    double frag_mz = frag_it->getMZ();

    double mz_offset = fragment_mass_error;

    if (mz_error_unit_ == "ppm")
    {
      mz_offset = frag_mz * 1e-6 * fragment_mass_error;
    }

    MSSpectrum::iterator db_mass_it = db_spectrum.MZBegin(frag_mz - mz_offset);
    MSSpectrum::iterator db_mass_end = db_spectrum.MZEnd(frag_mz + mz_offset);

    std::pair<double, Peak1D> nearest_peak(mz_offset + 1.0, Peak1D());

    // linear search for peak nearest to observed fragment peak
    for (; db_mass_it != db_mass_end; ++db_mass_it)
    {
      double db_mz(db_mass_it->getMZ());
      double abs_mass_diff(std::abs(frag_mz - db_mz));

      if (abs_mass_diff < nearest_peak.first) {
        nearest_peak.first = abs_mass_diff;
        nearest_peak.second = *db_mass_it;
      }
    }

    // update dot product
    if (nearest_peak.second.getIntensity() > 0.0)
    {
      ++matched_ions_count;
      dot_product += frag_it->getIntensity() * nearest_peak.second.getIntensity();
    }
  }

  double matched_ions_term(0.0);

  // return score 0 if too few matched ions
  if (matched_ions_count < 3)
  {
    return matched_ions_term;
  }


  if (matched_ions_count <= boost::math::max_factorial<double>::value)
  {
    matched_ions_term = std::log(boost::math::factorial<double>((double)matched_ions_count));
  }
  else
  {
    matched_ions_term = std::log(boost::math::factorial<double>(boost::math::max_factorial<double>::value));
  }

  double hyperscore(std::log(dot_product) + matched_ions_term);


  if (hyperscore < 0)
  {
    hyperscore = 0;
  }

  return hyperscore;
}

void MetaboliteSpectralMatching::run(PeakMap & msexp, PeakMap & spec_db, MzTab& mztab_out)
{
  std::sort(spec_db.begin(), spec_db.end(), PrecursorMZLess);

  std::vector<double> mz_keys;

  // copy precursor m/z values to vector for searching
  for (Size spec_idx = 0; spec_idx < spec_db.size(); ++spec_idx)
  {
    mz_keys.push_back(spec_db[spec_idx].getPrecursors()[0].getMZ());
  }

  // remove potential noise peaks by selecting the ten most intense peak per 100 Da window
  WindowMower wm;
  Param wm_param;

  wm_param.setValue("windowsize", 20.0);
  wm_param.setValue("movetype", "slide");
  wm_param.setValue("peakcount", 5);
  wm.setParameters(wm_param);

  wm.filterPeakMap(msexp);

  // merge MS2 spectra with same precursor mass
  SpectraMerger spme;
  spme.mergeSpectraPrecursors(msexp);
  wm.filterPeakMap(msexp);


  // container storing results
  std::vector<SpectralMatch> matching_results;


  for (Size spec_idx = 0; spec_idx < msexp.size(); ++spec_idx)
  {
    // std::cout << "merged spectrum no. " << spec_idx << " with #fragment ions: " << msexp[spec_idx].size() << std::endl;

    // iterate over all precursor masses
    for (Size prec_idx = 0; prec_idx < msexp[spec_idx].getPrecursors().size(); ++prec_idx)
    {
      // get precursor m/z
      double precursor_mz(msexp[spec_idx].getPrecursors()[prec_idx].getMZ());

      // std::cout << "precursor no. " << prec_idx << ": mz " << precursor_mz << " ";

      double prec_mz_lowerbound, prec_mz_upperbound;

      if (mz_error_unit_ == "Da")
      {
        prec_mz_lowerbound = precursor_mz - precursor_mz_error_;
        prec_mz_upperbound = precursor_mz + precursor_mz_error_;
      }
      else
      {
        double ppm_offset(precursor_mz * 1e-6 * precursor_mz_error_);
        prec_mz_lowerbound = precursor_mz - ppm_offset;
        prec_mz_upperbound = precursor_mz + ppm_offset;
      }


      // std::cout << "lower mz: " << prec_mz_lowerbound << " ";
      // std::cout << "upper mz: " << prec_mz_upperbound << std::endl;

      std::vector<double>::const_iterator lower_it = std::lower_bound(mz_keys.begin(), mz_keys.end(), prec_mz_lowerbound);
      std::vector<double>::const_iterator upper_it = std::upper_bound(mz_keys.begin(), mz_keys.end(), prec_mz_upperbound);

      Size start_idx(lower_it - mz_keys.begin());
      Size end_idx(upper_it - mz_keys.begin());

      //std::cout << "identifying " << msexp[spec_idx].getMetaValue("Massbank_Accession_ID") << std::endl;

      std::vector<SpectralMatch> partial_results;

      for (Size search_idx = start_idx; search_idx < end_idx; ++search_idx)
      {
        // do spectral matching
        // std::cout << "scanning " << spec_db[search_idx].getPrecursors()[0].getMZ() << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << std::endl;

        // check for charge state of precursor ions: do they match?
        if ( (ion_mode_ == "positive" && spec_db[search_idx].getPrecursors()[0].getCharge() < 0) || (ion_mode_ == "negative" && spec_db[search_idx].getPrecursors()[0].getCharge() > 0))
        {
          continue;
        }

        double hyperscore(computeHyperScore(msexp[spec_idx], spec_db[search_idx], fragment_mz_error_, 0.0));

        // std::cout << " scored with " << hyperScore << std::endl;
        if (hyperscore > 0)
        {
          // std::cout << "  ** detected " << spec_db[search_idx].getMetaValue("Massbank_Accession_ID") << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << " scored with " << hyperscore << std::endl;

          // score result temporarily
          SpectralMatch tmp_match;
          tmp_match.setObservedPrecursorMass(precursor_mz);
          tmp_match.setFoundPrecursorMass(spec_db[search_idx].getPrecursors()[0].getMZ());
          double obs_rt = std::floor(msexp[spec_idx].getRT() * 10)/10.0;
          tmp_match.setObservedPrecursorRT(obs_rt);
          tmp_match.setFoundPrecursorCharge(spec_db[search_idx].getPrecursors()[0].getCharge());
          tmp_match.setMatchingScore(hyperscore);
          tmp_match.setObservedSpectrumIndex(spec_idx);
          tmp_match.setMatchingSpectrumIndex(search_idx);

          tmp_match.setPrimaryIdentifier(spec_db[search_idx].getMetaValue("Massbank_Accession_ID"));
          tmp_match.setSecondaryIdentifier(spec_db[search_idx].getMetaValue("HMDB_ID"));
          tmp_match.setSumFormula(spec_db[search_idx].getMetaValue("Sum_Formula"));
          tmp_match.setCommonName(spec_db[search_idx].getMetaValue("Metabolite_Name"));
          tmp_match.setInchiString(spec_db[search_idx].getMetaValue("Inchi_String"));
          tmp_match.setSMILESString(spec_db[search_idx].getMetaValue("SMILES_String"));
          tmp_match.setPrecursorAdduct(spec_db[search_idx].getMetaValue("Precursor_Ion"));


          partial_results.push_back(tmp_match);

        }
      }

      // sort results by decreasing store
      std::sort(partial_results.begin(), partial_results.end(), SpectralMatchScoreGreater);

      // report mode: top3 or best?
      if (report_mode_ == "top3")
      {
        Size num_results(partial_results.size());

        Size last_result_idx = (num_results >= 3) ? 3 : num_results;

        for (Size result_idx = 0; result_idx < last_result_idx; ++result_idx)
        {
          // std::cout << "score: " << partial_results[result_idx].getMatchingScore() << " " << partial_results[result_idx].getMatchingSpectrumIndex() << std::endl;
          matching_results.push_back(partial_results[result_idx]);
        }
      }

      if (report_mode_ == "best")
      {
        if (partial_results.size() > 0)
        {
          matching_results.push_back(partial_results[0]);
        }
      }

    } // end precursor loop
  } // end spectra loop

  // write final results to MzTab
  exportMzTab_(matching_results, mztab_out);
}

/// protected methods

void MetaboliteSpectralMatching::updateMembers_()
{
  precursor_mz_error_ = (double)param_.getValue("prec_mass_error_value");
  fragment_mz_error_ = (double)param_.getValue("frag_mass_error_value");
  ion_mode_ = (String)param_.getValue("ionization_mode");

  mz_error_unit_ = (String)param_.getValue("mass_error_unit");
  report_mode_ = (String)param_.getValue("report_mode");
}


/// private methods

void MetaboliteSpectralMatching::exportMzTab_(const std::vector<SpectralMatch>& overall_results, MzTab& mztab_out)
{
  // iterate the overall results table
  MzTabSmallMoleculeSectionRows all_sm_rows;

  for (Size id_idx = 0; id_idx < overall_results.size(); ++id_idx)
  {
    SpectralMatch current_id(overall_results[id_idx]);

    MzTabSmallMoleculeSectionRow mztab_row_record;

    // set the identifier field
    String hid_temp = current_id.getPrimaryIdentifier();
    MzTabString prim_id;
    prim_id.set(hid_temp);
    std::vector<MzTabString> id_dummy;
    id_dummy.push_back(prim_id);
    MzTabStringList string_dummy_list;
    string_dummy_list.set(id_dummy);

    mztab_row_record.identifier = string_dummy_list;

    // set the chemical formula field
    MzTabString chem_form;
    String form_temp = current_id.getSumFormula();
    chem_form.set(form_temp);

    mztab_row_record.chemical_formula = chem_form;

    // set the smiles field
    String smi_temp = current_id.getSMILESString();     // extract SMILES from struct mapping file
    MzTabString smi_string;
    smi_string.set(smi_temp);

    mztab_row_record.smiles = smi_string;

    // set the inchi_key field
    String inchi_temp = current_id.getInchiString();    // extract INCHIKEY from struct mapping file
    MzTabString inchi_key;
    inchi_key.set(inchi_temp);

    mztab_row_record.inchi_key = inchi_key;

    // set description field (we use it for the common name of the compound)
    String name_temp = current_id.getCommonName();
    MzTabString common_name;
    common_name.set(name_temp);

    mztab_row_record.description = common_name;


    // set mass_to_charge field (precursor mass here)
    double mz_temp = current_id.getFoundPrecursorMass();
    MzTabDouble mass_to_charge;
    mass_to_charge.set(mz_temp);

    mztab_row_record.exp_mass_to_charge = mass_to_charge;  //TODO: distinguish the experimental precursor mass and spectral library precursor mass (later should probably go into cv_opt_ column)


    // set charge field
    Int ch_temp = current_id.getFoundPrecursorCharge();
    MzTabDouble mcharge;
    mcharge.set(ch_temp);

    mztab_row_record.charge = mcharge;


    // set RT field
    double rt_temp = current_id.getObservedPrecursorRT();
    MzTabDouble rt_temp2;
    rt_temp2.set(rt_temp);
    std::vector<MzTabDouble> rt_temp3;
    rt_temp3.push_back(rt_temp2);
    MzTabDoubleList observed_rt;
    observed_rt.set(rt_temp3);

    mztab_row_record.retention_time = observed_rt;


    // set database field
    String dbname_temp = "MassBank";
    MzTabString dbname;
    dbname.set(dbname_temp);

    mztab_row_record.database = dbname;


    // set database_version field
    String dbver_temp = "Sep 27, 2013";
    MzTabString dbversion;
    dbversion.set(dbver_temp);

    mztab_row_record.database_version = dbversion;


    // set smallmolecule_abundance_sub
    // check if we deal with a feature or consensus feature
    std::vector<MzTabDouble> int_temp3;

    double int_temp(0.0);
    MzTabDouble int_temp2;
    int_temp2.set(int_temp);
    int_temp3.push_back(int_temp2);

    for (Size i = 0; i != int_temp3.size(); ++i)
    {
      mztab_row_record.smallmolecule_abundance_study_variable[i + 1] = int_temp3[i];
    }

    // set smallmolecule_abundance_stdev_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
    double stdev_temp(0.0);
    MzTabDouble stdev_temp2;
    stdev_temp2.set(stdev_temp);
    std::vector<MzTabDouble> stdev_temp3;

    stdev_temp3.push_back(stdev_temp2);

    for (Size i = 0; i != stdev_temp3.size(); ++i)
    {
      mztab_row_record.smallmolecule_abundance_stdev_study_variable[i + 1] = stdev_temp3[i];
    }

    // set smallmolecule_abundance_std_error_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
    double stderr_temp(0.0);
    MzTabDouble stderr_temp2;
    stderr_temp2.set(stderr_temp);
    std::vector<MzTabDouble> stderr_temp3;

    stderr_temp3.push_back(stderr_temp2);

    for (Size i = 0; i != stderr_temp3.size(); ++i)
    {
      mztab_row_record.smallmolecule_abundance_std_error_study_variable[i + 1] = stderr_temp3[i];
    }

    // optional columns:
    std::vector<MzTabOptionalColumnEntry> optionals;

    // ppm error
    double error_ppm(((current_id.getFoundPrecursorMass() - current_id.getObservedPrecursorMass())/current_id.getFoundPrecursorMass())*1e6);
    error_ppm = std::floor(error_ppm*100)/100;

    MzTabString ppmerr;
    ppmerr.set(String(error_ppm));
    MzTabOptionalColumnEntry col0;
    col0.first = "opt_ppm_error";
    col0.second = ppmerr;
    optionals.push_back(col0);

    // set found adduct ion
    String addion_temp = current_id.getPrecursorAdduct();
    MzTabString addion;
    addion.set(addion_temp);
    MzTabOptionalColumnEntry col1;
    col1.first = "opt_adduct_ion";
    col1.second = addion;
    optionals.push_back(col1);

    // set isotope similarity score
    double sim_score_temp = current_id.getMatchingScore();
    std::stringstream read_in;
    read_in << sim_score_temp;
    String sim_score_temp2(read_in.str());
    MzTabString sim_score;
    sim_score.set(sim_score_temp2);
    MzTabOptionalColumnEntry col2;
    col2.first = "opt_match_score";
    col2.second = sim_score;
    optionals.push_back(col2);


    // set secondary ID (here HMDB id)
    String sec_id = current_id.getSecondaryIdentifier();
    MzTabString sec_id_str;
    sec_id_str.set(sec_id);
    MzTabOptionalColumnEntry col3;
    col3.first = "opt_sec_id";
    col3.second = sec_id_str;
    optionals.push_back(col3);

    // set source spectra index
    // TODO: this should use spectra_ref column
    String source_idx = String(current_id.getObservedSpectrumIndex());
    MzTabString source_idx_str;
    source_idx_str.set(source_idx);
    MzTabOptionalColumnEntry col4;
    col4.first = "opt_source_idx";
    col4.second = source_idx_str;
    optionals.push_back(col4);

    mztab_row_record.opt_ = optionals;

    all_sm_rows.push_back(mztab_row_record);
  }

  mztab_out.setSmallMoleculeSectionRows(all_sm_rows);

}

} // closing namespace OpenMS
