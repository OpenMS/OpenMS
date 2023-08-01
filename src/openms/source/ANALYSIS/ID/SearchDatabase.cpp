// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Chris Bielow $
// $Authors: Max Alcer, Heike Einsfeld $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/ID/SearchDatabase.h>

#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>

#include <omp.h>
#include <unordered_set>
#include <vector>

using namespace std;


namespace OpenMS
{

using namespace Constants::UserParam;
using namespace Math;

std::vector<SearchDatabase::Peptide_> SearchDatabase::generate_peptides_(const std::vector<FASTAFile::FASTAEntry>& entries) const
{
  vector<SearchDatabase::Peptide_> all_peptides;

  ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
  ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);

  #pragma omp parallel
  { 
    ProteaseDigestion digestor;
    digestor.setEnzyme(digestor_enzyme_);
    digestor.setMissedCleavages(missed_cleavages_);
    vector<SearchDatabase::Peptide_> all_peptides_pvt;
    #pragma omp for nowait
    for (size_t i = 0; i < entries.size(); i++)
    {
      vector<AASequence> peptides;
      digestor.digest(AASequence::fromString(entries[i].sequence), peptides, peptide_min_length_, peptide_max_length_);
      
      for (const auto& pep : peptides)
      { 
        if (pep.toString().find('X') != string::npos) continue; // filtering peptide with unknown AA, can't calculate MonoWeight

        vector<AASequence> modified_peptides;

        AASequence modified_pep = pep;

        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, modified_pep);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, modified_pep, max_variable_mods_per_peptide_, modified_peptides);

        for (const auto & mod_pep : modified_peptides)
        {
          double mod_seq_mz = mod_pep.getMonoWeight();

          if (!contains(mod_seq_mz, peptide_min_mass_, peptide_max_mass_)) continue;

          all_peptides_pvt.emplace_back(mod_pep, i, mod_seq_mz);
        }

        double seq_mz = pep.getMonoWeight();

        if (!contains(seq_mz, peptide_min_mass_, peptide_max_mass_)) continue;

        all_peptides_pvt.emplace_back(pep, i, seq_mz);
        
      }
    }
    #pragma omp critical (add_critical)
    all_peptides.insert(all_peptides.end(), all_peptides_pvt.begin(), all_peptides_pvt.end());
  }
  return all_peptides;
}

void SearchDatabase::fragment_merge_(int first, int last, const std::vector<int>& chunks, std::vector<SearchDatabase::Fragment_>& input) const
{
  if (last - first > 1)
  {
    int mid = first + (last-first) / 2;
    #pragma omp parallel sections
    {
      #pragma omp section
      SearchDatabase::fragment_merge_(first, mid, chunks, input);
      #pragma omp section
      SearchDatabase::fragment_merge_(mid, last, chunks, input);
    }
    std::inplace_merge(input.begin() + chunks[first], input.begin() + chunks[mid], input.begin() + chunks[last], 
    [&](const  SearchDatabase::Fragment_& l, const  SearchDatabase::Fragment_& r)-> bool
    {return (tie(l.fragment_mz_, all_peptides_[l.peptide_index_].peptide_mz_) < tie(r.fragment_mz_, all_peptides_[r.peptide_index_].peptide_mz_));});
  }
}

std::vector<SearchDatabase::Fragment_> SearchDatabase::generate_fragments_() const
{
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  std::vector<Fragment_> all_frags;
  std::vector<int> chunk_start = {0}; //for fragment_merge need start and end of presorted chunks
    
  for (size_t i = 0; i < all_peptides_.size(); i++)
  { 
    tsg.getSpectrum(b_y_ions, all_peptides_[i].sequence_, 1, 1);      
    for (const auto& frag : b_y_ions)
    { 
      if (!contains(frag.getMZ(), fragment_min_mz_, fragment_max_mz_)) continue;
      all_frags.emplace_back(i, frag);        
    }
    chunk_start.emplace_back(all_frags.size());
    b_y_ions.clear(true);
  }

  SearchDatabase::fragment_merge_(0, chunk_start.size()-1, chunk_start, all_frags);
  
  return all_frags;
}

SearchDatabase::SearchDatabase(const std::vector<FASTAFile::FASTAEntry>& entries) : DefaultParamHandler("SearchDatabase")
{ 
  vector<String> all_enzymes;  
  ProteaseDB::getInstance()->getAllNames(all_enzymes);
  vector<String> all_mods;
  ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
  vector<string> tolerance_units{UNIT_DA, UNIT_PPM};
  defaults_.setValue("digestor_enzyme", "Trypsin", "Enzyme for digestion");
  defaults_.setValidStrings("digestor_enzyme", ListUtils::create<std::string>(all_enzymes));
  defaults_.setValue("missed_cleavages", 0, "Missed cleavages for digestion");
  defaults_.setValue("peptide_min_mass", 500, "Minimal peptide mass for database");
  defaults_.setValue("peptide_max_mass", 5000, "Maximal peptide mass for database");
  defaults_.setValue("peptide_min_length", 5, "Minimal peptide length for database");
  defaults_.setValue("peptide_max_length", 50, "Maximal peptide length for database");
  defaults_.setValue("fragment_min_mz", 150, "Minimal fragment mz for database");
  defaults_.setValue("fragment_max_mz", 2000, "Maximal fragment mz for database");
  defaults_.setValue("precursor_mz_tolerance", 2.0, "Tolerance for precursor-m/z in search");
  defaults_.setValue("fragment_mz_tolerance", 0.05, "Tolerance for fragment-m/z in search");
  defaults_.setValue("precursor_mz_tolerance_unit", UNIT_DA, "Unit of tolerance for precursor-m/z");
  defaults_.setValidStrings("precursor_mz_tolerance_unit", tolerance_units);
  defaults_.setValue("fragment_mz_tolerance_unit", UNIT_DA, "Unit of tolerance for fragment-m/z");
  defaults_.setValidStrings("fragment_mz_tolerance_unit", tolerance_units);
  defaults_.setValue("modifications_fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
  defaults_.setValidStrings("modifications_fixed", ListUtils::create<std::string>(all_mods));
  defaults_.setValue("modifications_variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
  defaults_.setValidStrings("modifications_variable", ListUtils::create<std::string>(all_mods));
  defaults_.setValue("max_variable_mods_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");  

  defaultsToParam_();

  if (entries.size() == 0) return;

  all_peptides_ = generate_peptides_(entries);

  all_fragments_ = generate_fragments_();
  
  bucketsize_ = size_t(sqrt(all_fragments_.size())); //calculating bucketsize for balanced tree (performance)

  for (size_t i = 0; i < all_fragments_.size(); i += bucketsize_)
  {
    bucket_frags_mz_.emplace_back(all_fragments_[i].fragment_mz_);
  }
  // sorting the fragments of each bucket by the precursor-m/z
  #pragma omp parallel for  
  for (size_t i = 0; i < all_fragments_.size(); i += bucketsize_)
  { 
    auto bucket_begin = all_fragments_.begin() + i;
    auto condition = distance(all_fragments_.begin(), bucket_begin + bucketsize_) >= int(all_fragments_.size());

    sort(bucket_begin, condition ? all_fragments_.end() : bucket_begin+bucketsize_, 
    [&](const SearchDatabase::Fragment_& l, const SearchDatabase::Fragment_& r) -> bool 
    {return (all_peptides_[l.peptide_index_].peptide_mz_ < all_peptides_[r.peptide_index_].peptide_mz_);});
  }
}

void SearchDatabase::updateMembers_()
{
  digestor_enzyme_ = param_.getValue("digestor_enzyme").toString();
  missed_cleavages_ = param_.getValue("missed_cleavages");
  peptide_min_mass_ = param_.getValue("peptide_min_mass");
  peptide_max_mass_ = param_.getValue("peptide_max_mass");
  peptide_min_length_ = param_.getValue("peptide_min_length");
  peptide_max_length_ = param_.getValue("peptide_max_length");
  fragment_min_mz_ = param_.getValue("fragment_min_mz");
  fragment_max_mz_ = param_.getValue("fragment_max_mz");
  precursor_mz_tolerance_ = param_.getValue("precursor_mz_tolerance");
  fragment_mz_tolerance_ = param_.getValue("fragment_mz_tolerance");
  precursor_mz_tolerance_unit_ = param_.getValue("precursor_mz_tolerance_unit").toString();
  fragment_mz_tolerance_unit_ = param_.getValue("fragment_mz_tolerance_unit").toString();
  modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_fixed"));
  modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_variable"));
  max_variable_mods_per_peptide_ = param_.getValue("max_variable_mods_per_peptide");
}

void SearchDatabase::search(MSSpectrum& spectrum, std::vector<SearchDatabase::Candidate>& candidates) const
{ 
  candidates.clear();
  if (spectrum.empty()) return;
  
  unordered_set<size_t> index_hash; // saving every candidate only once

  std::vector<Precursor> precursor = spectrum.getPrecursors();

  if (precursor.size() != 1) return;

  double prec_mz = precursor[0].getUnchargedMass();

  spectrum.sortByIntensity(true);
  size_t count_found = 0;
  size_t count_rounds = 0;

  auto prec_unit_da = precursor_mz_tolerance_unit_ == UNIT_DA;
  auto frag_unit_da = fragment_mz_tolerance_unit_ == UNIT_DA;

  for (const auto& peak : spectrum)
  {
    if (count_found == 0 && count_rounds == 5) return; // if peak with the 5 highest intensity don't match, skip whole spectrum

    double new_frag_mz_tolerance = frag_unit_da ? fragment_mz_tolerance_ : ppmToMassAbs(fragment_mz_tolerance_, peak.getMZ());

    auto itr_lower = lower_bound(bucket_frags_mz_.begin(), bucket_frags_mz_.end(), peak.getMZ() -  new_frag_mz_tolerance);
    
    size_t index_lower = distance(bucket_frags_mz_.begin(), itr_lower);

    if (index_lower != 0) index_lower--; // lower bound returns one too high (because searching the minimum of each bucket)

    while (bucket_frags_mz_[index_lower] <= (peak.getMZ() + new_frag_mz_tolerance) && index_lower < bucket_frags_mz_.size())
    { 
      auto bucket_start = all_fragments_.begin() + (index_lower * bucketsize_);
      auto bucket_end = all_fragments_.begin() + (index_lower + 1) * bucketsize_;

      auto condition = distance(all_fragments_.begin(), bucket_end) >= int(all_fragments_.size());

      double new_prec_mz_tolerance = prec_unit_da ? precursor_mz_tolerance_ : ppmToMassAbs(precursor_mz_tolerance_, prec_mz);

      auto itr_lower_inner = lower_bound(bucket_start, condition ? all_fragments_.end() : bucket_end, prec_mz - new_prec_mz_tolerance,
      [&](const SearchDatabase::Fragment_& r, double l) {
      return (all_peptides_[r.peptide_index_].peptide_mz_ < l);
      });

      while (all_peptides_[(*itr_lower_inner).peptide_index_].peptide_mz_ <= (prec_mz + new_prec_mz_tolerance) && itr_lower_inner != bucket_end && itr_lower_inner != all_fragments_.end())
      { 
        double theoretical_peak_tolerance = frag_unit_da ? fragment_mz_tolerance_ : ppmToMassAbs(fragment_mz_tolerance_, (*itr_lower_inner).fragment_mz_);
        // filter by matching fragment-m/z
        if (!contains(peak.getMZ(), (*itr_lower_inner).fragment_mz_ - theoretical_peak_tolerance, (*itr_lower_inner).fragment_mz_ + theoretical_peak_tolerance)) 
        {
          ++itr_lower_inner;
          continue;
        }
        index_hash.insert((*itr_lower_inner).peptide_index_);
        ++itr_lower_inner;
        count_found++;
      }
      index_lower++;
    }    
    count_rounds++;
  }
  
  for (size_t j : index_hash)
  {
    candidates.emplace_back(&all_peptides_[j].sequence_, all_peptides_[j].protein_index_);
  }
}

void SearchDatabase::search(MSExperiment& experiment, std::vector<SearchDatabase::CandidatesWithIndex>& candidates) const
{
  candidates.clear();
  if (experiment.size() == 0) return;  
  
  #pragma omp parallel
  { 
    #pragma omp for 
    for (size_t i = 0; i < experiment.size(); i++)
    {
      if (experiment[i].empty()) continue;

      vector<Candidate> temp_cand;
      
      SearchDatabase::search(experiment[i], temp_cand);

      #pragma omp critical (add_critical)
      {
        candidates.emplace_back(temp_cand, i);
      }
    }
  }
}

} // end namespace OpenMS





