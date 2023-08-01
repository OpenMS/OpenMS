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


#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <vector>

namespace OpenMS
{
  //class MSExperiment;

  /** @brief Creating a two-level tree from a given FASTAFile by generating the theoretical peptides and their b-y-ions and sorting
    * them in the outer tree by fragment-m/z and the inner tree by precursor-m/z. Search function to find candidates to experimental MS2 sepctra
    * (implementation inspired by https://lazear.github.io/sage/)
    * @ingroup ID
   */
class OPENMS_DLLAPI FragmentIndex : public DefaultParamHandler
{
  public:

  /** @brief Storing found candiates from search 
   */
  struct Candidate
  {
    size_t peptide_start; ///< Start position of peptide in protein
    size_t peptide_end; ///< End position of peptide in protein
    size_t protein_index; ///< Index to Entry of std::vector<FASTAFile::FASTAEntry>
    bool is_modified_; ///< true if candidate is modified
    size_t modification_index_; ///< index to modified version of this peptide
  };
  
  /** @brief Storing vector of found candiates from search with the Index of the MSSpectrum in the MSExperiment
   */
  struct CandidatesWithIndex
  {
    std::vector<FragmentIndex::Candidate> candidates; ///< Vector of found candidates
    size_t spectrum_index; ///< Index to the MSSpectrum in the MSExperiment
    CandidatesWithIndex(std::vector<FragmentIndex::Candidate>&& cand, size_t spectrum_i): candidates(std::move(cand)), spectrum_index(spectrum_i){}
  };


  /** @brief Sets default parameters     
  */
  FragmentIndex();

  /** @brief Builds up the Search Datastructure     
  * @param entries Input vector of FASTAFile-Entries to base SearchDatastructure on
  */
  void build(const std::vector<FASTAFile::FASTAEntry>& entries);

  /** @brief Searches Peaks of every MSSpectrum of the MSExperiment in the Database       
     * @param experiment Input MSExperiment containing of MS2 Spectra (MS1 Spectra has an empty output)
     * @param candidates Output vector of found candidates with the index of MSSpectrum in experiment
    */
  void search(MSExperiment& experiment, std::vector<FragmentIndex::CandidatesWithIndex>& candidates) const;

  /** @brief Searches Peaks of the MSSpectrum in the Database       
     * @param spectrum Input MS2 Spectrum (MS1 Spectra has an empty output)
     * @param candidates Output vector of found candidatest
    */
  void search(MSSpectrum& spectrum, std::vector<FragmentIndex::Candidate>& candidates) const;

  protected:
  
  void updateMembers_() override;
  
  //Saving b_y_ion after creating theoretical spectrum with index to precursor-peptide
  struct Fragment_
  {
    size_t peptide_index_;
    double fragment_mz_;
    Fragment_(size_t prec, const Peak1D& frag):peptide_index_(prec), fragment_mz_(frag.getMZ()){}
  };

  //Saving first and last index of precursor in FASTA-Entry, index to protein and the MonoWeight after digesting
  struct Peptide_
  {
    size_t peptide_begin_;
    size_t peptide_end_;
    size_t protein_index_;
    double peptide_mz_;
    bool is_modified_;
    size_t modification_index_;
  };
  
  /// in-silico digest protein database
  std::vector<FragmentIndex::Peptide_> generate_peptides_(const std::vector<FASTAFile::FASTAEntry>& entries) const;
  /// Merges presorted Chunks of Peptide-Fragments inplace
  void fragment_merge_(int first, int last, const std::vector<int>& chunks, std::vector<FragmentIndex::Fragment_>& input) const;
  ///generates sorted vector with all theoretical Fragments for all theoretical Peptides
  std::vector<FragmentIndex::Fragment_> generate_fragments_(const std::vector<FASTAFile::FASTAEntry>& entries) const;

  bool is_build_;
  std::string digestor_enzyme_;
  size_t missed_cleavages_; ///< number of missed cleavages
  double peptide_min_mass_;
  double peptide_max_mass_;
  size_t peptide_min_length_;
  size_t peptide_max_length_;
  double fragment_min_mz_;
  double fragment_max_mz_;
  size_t fragment_min_charge_;
  size_t fragment_max_charge_;
  size_t bucketsize_; ///< number of fragments per outer node
  double precursor_mz_tolerance_;
  std::string precursor_mz_tolerance_unit_;
  double fragment_mz_tolerance_;
  std::string fragment_mz_tolerance_unit_;
  StringList modifications_fixed_;
  StringList modifications_variable_;
  size_t max_variable_mods_per_peptide_;
  size_t max_missed_peaks_;
  std::vector<Peptide_> all_peptides_;
  std::vector<double> bucket_frags_mz_; ///< Minimum fragment-m/z for each other node
  std::vector<Fragment_> all_fragments_;
};

}