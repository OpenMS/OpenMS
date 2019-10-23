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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

// turn on additional debug output
// #define DEBUG_OPXLHELPER


using namespace std;

namespace OpenMS
{
  // Enumerate all pairs of peptides from the searched database and calculate their masses (inlcuding mono-links and loop-links)
  vector<OPXLDataStructs::XLPrecursor> OPXLHelper::enumerateCrossLinksAndMasses(const vector<OPXLDataStructs::AASeqWithMass>& peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, const vector< double >& spectrum_precursors, vector< int >& precursor_correction_positions, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
  {
    // initialize empty vector for the results
    vector<OPXLDataStructs::XLPrecursor> mass_to_candidates;

    double min_precursor = spectrum_precursors[0];
    double max_precursor = spectrum_precursors[spectrum_precursors.size()-1];

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize p1 = 0; p1 < static_cast<SignedSize>(peptides.size()); ++p1)
    {
      // get the amino acid sequence of this peptide as a character string
      String seq_first = peptides[p1].peptide_seq.toUnmodifiedString();

      // generate mono-links: one cross-linker with one peptide attached to one side
      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
      {
        // Monoisotopic weight of the peptide + cross-linker
        double cross_linked_pair_mass = peptides[p1].peptide_mass + cross_link_mass_mono_link[i];

        // Make sure it is clear only one peptide is considered here. Use an out-of-range value for the second peptide.
        // to check: if(precursor.beta_index < peptides.size()) returns "false" for a mono-link
        OPXLDataStructs::XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = peptides.size() + 1; // an out-of-range index to represent an empty index

        // call function to compare with spectrum precursor masses
        // will only add this candidate, if the mass is within the given tolerance to any precursor in the spectra data
        // after the first monolink is added, stop enumerating masses (if other candidates fit within the same precursor, they will have exactly the same fragment matching)
        if (filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_correction_positions, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor))
        {
          break;
        }
      }

       // test if this peptide could have loop-links: one cross-link with both sides attached to the same peptide
       // TODO check for distance between the two linked residues
      bool first_res = false; // is there a residue the first side of the linker can attach to?
      bool second_res = false; // is there a residue the second side of the linker can attach to?
      for (Size k = 0; k < seq_first.size()-1; ++k)
      {
        for (Size i = 0; i < cross_link_residue1.size(); ++i)
        {
          if (cross_link_residue1[i].size() == 1 && seq_first.substr(k, 1) == cross_link_residue1[i])
          {
            first_res = true;
          }
        }
        for (Size i = 0; i < cross_link_residue2.size(); ++i)
        {
          if (cross_link_residue2[i].size() == 1 && seq_first.substr(k, 1) == cross_link_residue2[i])
          {
            second_res = true;
          }
        }
      }

      // If both sides of a cross-linker can link to this peptide, generate the loop-link
      if (first_res && second_res)
      {
        // Monoisotopic weight of the peptide + cross-linker
        double cross_linked_pair_mass = peptides[p1].peptide_mass + cross_link_mass;

        // also only one peptide
        OPXLDataStructs::XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = peptides.size() + 1; // an out-of-range index to represent an empty index

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_correction_positions, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }

      // check for minimal mass of second peptide, jump farther than current peptide if possible
      double allowed_error = 0;
      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = min_precursor * precursor_mass_tolerance * 1e-6;
      }
      else // Dalton
      {
        allowed_error = precursor_mass_tolerance;
      }
      double min_second_peptide_mass = min_precursor - cross_link_mass - peptides[p1].peptide_mass - allowed_error;

      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = max_precursor * precursor_mass_tolerance * 1e-6;
      }
      double max_second_peptide_mass = max_precursor - cross_link_mass - peptides[p1].peptide_mass + allowed_error;

      // Generate cross-links: one cross-linker linking two separate peptides, the most important case
      // Loop over all p2 peptide candidates, that come after p1 in the list
      for (Size p2 = p1; p2 < peptides.size(); ++p2)
      {
        // skip peptides, that are too small in any case
        if (peptides[p2].peptide_mass < min_second_peptide_mass)
        {
          continue;
        }
        else if (peptides[p2].peptide_mass > max_second_peptide_mass)
        {
          break;
        }

        // Monoisotopic weight of the first peptide + the second peptide + cross-linker
        double cross_linked_pair_mass = peptides[p1].peptide_mass + peptides[p2].peptide_mass + cross_link_mass;

        // this time both peptides have valid indices
        OPXLDataStructs::XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = p2;

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_correction_positions, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }
    } // end of parallelized for-loop
    return mass_to_candidates;
  }

  bool OPXLHelper::filter_and_add_candidate(vector<OPXLDataStructs::XLPrecursor>& mass_to_candidates, const vector< double >& spectrum_precursors, vector< int >& precursor_correction_positions, bool precursor_mass_tolerance_unit_ppm, double precursor_mass_tolerance, OPXLDataStructs::XLPrecursor precursor)
  {
    vector< double >::const_iterator low_it;
    vector< double >::const_iterator up_it;

    // compute absolute tolerance from relative, if necessary
    double allowed_error = 0;
    if (precursor_mass_tolerance_unit_ppm) // ppm
    {
      allowed_error = precursor.precursor_mass * precursor_mass_tolerance * 1e-6;
    }
    else // Dalton
    {
      allowed_error = precursor_mass_tolerance;
    }

    // find precursor with m/z >= low end of range
    low_it = lower_bound(spectrum_precursors.begin(), spectrum_precursors.end(), precursor.precursor_mass - allowed_error);
    // find precursor with m/z > (not equal to) high end of range
    up_it =  upper_bound(spectrum_precursors.begin(), spectrum_precursors.end(), precursor.precursor_mass + allowed_error);
    // if these two are equal, there is no precursor within the range

    if (low_it != up_it) // if they are not equal, there are matching precursors in the data
    {
#pragma omp critical (mass_to_candidates_access)
      {
        mass_to_candidates.push_back(precursor);
        // take the position of the highest matching precursor mass in the vector (prioritize smallest correction)
        precursor_correction_positions.push_back(std::distance(spectrum_precursors.begin(), std::prev(up_it, 1)));
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  std::vector<OPXLDataStructs::AASeqWithMass> OPXLHelper::digestDatabase(
    vector<FASTAFile::FASTAEntry> fasta_db,
    EnzymaticDigestion digestor,
    Size min_peptide_length,
    StringList cross_link_residue1,
    StringList cross_link_residue2,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
    Size max_variable_mods_per_peptide)
  {
    multimap<StringView, AASequence> processed_peptides;
    vector<OPXLDataStructs::AASeqWithMass> peptide_masses;

    bool n_term_linker = false;
    bool c_term_linker = false;
    for (const String& res : cross_link_residue1)
    {
      if (res == "N-term")
      {
        n_term_linker = true;
      }
      else if (res == "C-term")
      {
        c_term_linker = true;
      }
    }
    for (const String& res : cross_link_residue2)
    {
      if (res == "N-term")
      {
        n_term_linker = true;
      }
      else if (res == "C-term")
      {
        c_term_linker = true;
      }
    }

    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {
      // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)
      vector<StringView> current_digest;
      digestor.digestUnmodified(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        // skip peptides with invalid AAs // TODO is this necessary?
        if (cit->getString().has('B') || cit->getString().has('O') || cit->getString().has('U') || cit->getString().has('X') || cit->getString().has('Z')) continue;

        OPXLDataStructs::PeptidePosition position = OPXLDataStructs::INTERNAL;
        if (fasta_db[fasta_index].sequence.hasPrefix(cit->getString()))
        {
          position = OPXLDataStructs::N_TERM;
        }
        else if (fasta_db[fasta_index].sequence.hasSuffix(cit->getString()))
        {
          position = OPXLDataStructs::C_TERM;
        }

        // skip if no cross-linked residue
        bool skip = true;
        if (n_term_linker && position == OPXLDataStructs::N_TERM)
        {
          skip = false;
        }
        else if (c_term_linker && position == OPXLDataStructs::C_TERM)
        {
          skip = false;
        }
        else
        {
          for (const String& res : cross_link_residue1)
          {
            if (res.size() == 1 && (cit->getString().find(res) < cit->getString().size()-1))
            {
              skip = false;
            }
          }
          for (const String& res : cross_link_residue2)
          {
            if (res.size() == 1 && (cit->getString().find(res) < cit->getString().size()-1))
            {
              skip = false;
            }
          }
        }
        if (skip) continue;

        bool already_processed = false;

        if (processed_peptides.find(*cit) != processed_peptides.end())
        {
          // peptide (and all modified variants) already processed so skip it
          already_processed = true;
        }

        if (already_processed) continue;

        vector<AASequence> all_modified_peptides;

        // generate all modified variants of a peptide
        AASequence aas = AASequence::fromString(cit->getString());
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < static_cast<SignedSize>(all_modified_peptides.size()); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          OPXLDataStructs::AASeqWithMass pep_mass;
          pep_mass.peptide_mass = candidate.getMonoWeight();
          pep_mass.peptide_seq = candidate;
          pep_mass.position = position;

          processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
          peptide_masses.push_back(pep_mass);
        }
      }
    }
    sort(peptide_masses.begin(), peptide_masses.end(), OPXLDataStructs::AASeqWithMassComparator());
    return peptide_masses;
  }

  vector <OPXLDataStructs::ProteinProteinCrossLink> OPXLHelper::buildCandidates(const std::vector< OPXLDataStructs::XLPrecursor > & candidates,
                                                                                const std::vector< int > & precursor_corrections,
                                                                                const std::vector< int >& precursor_correction_positions,
                                                                                const std::vector<OPXLDataStructs::AASeqWithMass> & peptide_masses,
                                                                                const StringList & cross_link_residue1,
                                                                                const StringList & cross_link_residue2,
                                                                                double cross_link_mass,
                                                                                const DoubleList & cross_link_mass_mono_link,
                                                                                const std::vector< double >& spectrum_precursor_vector,
                                                                                const std::vector< double >& allowed_error_vector,
                                                                                String cross_link_name)
  {
    bool n_term_linker = false;
    bool c_term_linker = false;
    for (const String& res : cross_link_residue1)
    {
      if (res == "N-term")
      {
        n_term_linker = true;
      }
      else if (res == "C-term")
      {
        c_term_linker = true;
      }
    }
    for (const String& res : cross_link_residue2)
    {
      if (res == "N-term")
      {
        n_term_linker = true;
      }
      else if (res == "C-term")
      {
        c_term_linker = true;
      }
    }

    vector <OPXLDataStructs::ProteinProteinCrossLink> cross_link_candidates;

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      OPXLDataStructs::XLPrecursor candidate = candidates[i];
      vector <SignedSize> link_pos_first;
      vector <SignedSize> link_pos_second;
      const AASequence* peptide_first = &(peptide_masses[candidate.alpha_index].peptide_seq);
      OPXLDataStructs::PeptidePosition peptide_pos_first = peptide_masses[candidate.alpha_index].position;
      const AASequence* peptide_second = nullptr;
      OPXLDataStructs::PeptidePosition peptide_pos_second = OPXLDataStructs::INTERNAL;
      if (candidate.beta_index < peptide_masses.size())
      {
        peptide_second = &(peptide_masses[candidate.beta_index].peptide_seq);
        peptide_pos_second = peptide_masses[candidate.beta_index].position;
      }
      String seq_first = peptide_first->toUnmodifiedString();
      String seq_second;
      if (peptide_second) { seq_second = peptide_second->toUnmodifiedString(); }

      // mono-links and loop-links with different masses can be generated for the same precursor mass, but only one of them can be valid each time.
      // Find out which is the case. But it should not happen often enough to slow down the tool significantly.
      bool is_loop = false;
      for (Size f = 0; f < allowed_error_vector.size(); ++f)
      {
        if (abs(spectrum_precursor_vector[f] - (peptide_first->getMonoWeight() + cross_link_mass)) <= allowed_error_vector[f])
        {
          is_loop = true;
        }
      }

      for (Size k = 0; k < seq_first.size()-1; ++k)
      {
        for (Size x = 0; x < cross_link_residue1.size(); ++x)
        {
          if (seq_first.substr(k, 1) == cross_link_residue1[x]) link_pos_first.push_back(k);
        }
      }
      if (candidate.beta_index < peptide_masses.size())
      {
        for (Size k = 0; k < seq_second.size()-1; ++k)
        {
          for (Size x = 0; x < cross_link_residue2.size(); ++x)
          {
            if (seq_second.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
          }
        }
      }
      else
      {
        // Second position defining a mono-link and the second positions on the same peptide for loop links (only one of these two is valid for any specific precursor)
        if (!is_loop)
        {
          link_pos_second.push_back(-1);
        }
        else
        {
          for (Size k = 0; k < seq_first.size()-1; ++k)
          {
            for (Size x = 0; x < cross_link_residue2.size(); ++x)
            {
              if (seq_first.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
            }
          }
        }
      }

      // Determine larger peptide (alpha) by sequence length, use mass as tie breaker
      bool alpha_first = true;

      if (seq_second.size() > seq_first.size())
      {
        alpha_first = false;
      }
      else if ( (seq_second.size() == seq_first.size()) && peptide_second && (peptide_second->getMonoWeight() > peptide_first->getMonoWeight()) )
      {
        alpha_first = false;
      }

      // TODO remodel this, there should be a simpler way, e.g. the peptides were sorted so "second" is always heavier?
      // generate cross_links for all valid combinations
      for (Size x = 0; x < link_pos_first.size(); ++x)
      {
        for (Size y = 0; y < link_pos_second.size(); ++y)
        {
          OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate;
          cross_link_candidate.precursor_correction = precursor_corrections[i];
          cross_link_candidate.cross_linker_name = cross_link_name;
          // if loop link, and the positions are the same, then it is linking the same residue with itself,  skip this combination, also pos1 > pos2 would be the same link as pos1 < pos2
          if (((seq_second.size() == 0) && (link_pos_first[x] >= link_pos_second[y])) && (link_pos_second[y] != -1))
          {
            continue;
          }
          if (alpha_first)
          {
            cross_link_candidate.alpha = peptide_first;
            cross_link_candidate.beta = peptide_second;
            cross_link_candidate.cross_link_position.first = link_pos_first[x];
            cross_link_candidate.cross_link_position.second = link_pos_second[y];
            cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
            cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
          }
          else
          {
            cross_link_candidate.alpha = peptide_second;
            cross_link_candidate.beta = peptide_first;
            cross_link_candidate.cross_link_position.first = link_pos_second[y];
            cross_link_candidate.cross_link_position.second = link_pos_first[x];
            cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
            cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
          }
          // Cross-linker mass is only one of the mono-link masses, if there is no second position (second == -1), otherwise the normal linker mass
          if (link_pos_second[y] != -1)
          {
            cross_link_candidate.cross_linker_mass = cross_link_mass;

#pragma omp critical (cross_link_candidates_access)
            cross_link_candidates.push_back(cross_link_candidate);
          }
          else
          {
            for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
            {
              // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
              bool is_correct_monolink = false;
              if (abs(spectrum_precursor_vector[precursor_correction_positions[i]] - (peptide_first->getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error_vector[precursor_correction_positions[i]])
              {
                is_correct_monolink = true;
              }
              if (is_correct_monolink)
              {
                cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];

#pragma omp critical (cross_link_candidates_access)
                cross_link_candidates.push_back(cross_link_candidate);
              }
            }
          }
        }
      }

      if (peptide_pos_second != OPXLDataStructs::INTERNAL)
      {
        ResidueModification::TermSpecificity second_spec = ResidueModification::N_TERM;
        Size mod_pos = 0;
        bool compatible = false;
        // If the peptide is at the N-terminus of the protein and the first residue is otherwise not linkable, add a terminal cross-linker
        if (n_term_linker && (peptide_pos_second == OPXLDataStructs::N_TERM) && find(link_pos_second.begin(), link_pos_second.end(), 0) == link_pos_second.end())
        {
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_second == OPXLDataStructs::C_TERM) && peptide_second)
        {
          second_spec = ResidueModification::C_TERM;
          mod_pos = peptide_second->size() - 1;
          compatible = true;
        }
        if (compatible)
        {
          for (Size x = 0; x < link_pos_first.size(); ++x)
          {
            OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate;
            if (alpha_first)
            {
              cross_link_candidate.alpha = peptide_first;
              cross_link_candidate.beta = peptide_second;
              cross_link_candidate.cross_link_position.first = link_pos_first[x];
              cross_link_candidate.cross_link_position.second = mod_pos;
              cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
              cross_link_candidate.term_spec_beta = second_spec;
            }
            else
            {
              cross_link_candidate.alpha = peptide_second;
              cross_link_candidate.beta = peptide_first;
              cross_link_candidate.cross_link_position.first = mod_pos;
              cross_link_candidate.cross_link_position.second = link_pos_first[x];
              cross_link_candidate.term_spec_alpha = second_spec;
              cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
            }
            // If second peptide has a term specificity, there must be a second peptide, so we don't have to consider mono or loop-links
            cross_link_candidate.cross_linker_mass = cross_link_mass;
            cross_link_candidate.cross_linker_name = cross_link_name;
            cross_link_candidate.precursor_correction = precursor_corrections[i];

#pragma omp critical (cross_link_candidates_access)
            cross_link_candidates.push_back(cross_link_candidate);
          }
        }
      }

      if (peptide_pos_first != OPXLDataStructs::INTERNAL)
      {
        ResidueModification::TermSpecificity first_spec = ResidueModification::N_TERM;
        Size mod_pos = 0;
        bool compatible = false;
        // If the peptide is at the N-terminus of the protein and the first residue is otherwise not linkable, add a terminal cross-linker
        if (n_term_linker && (peptide_pos_first == OPXLDataStructs::N_TERM) && find(link_pos_first.begin(), link_pos_first.end(), 0) == link_pos_first.end())
        {
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_first == OPXLDataStructs::C_TERM))
        {
          first_spec = ResidueModification::C_TERM;
          mod_pos = peptide_first->size() - 1;
          compatible = true;
        }
        if (compatible)
        {
          for (Size x = 0; x < link_pos_second.size(); ++x)
          {
            OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate;
            cross_link_candidate.cross_linker_name = cross_link_name;
            cross_link_candidate.precursor_correction = precursor_corrections[i];
            if (alpha_first)
            {
              cross_link_candidate.alpha = peptide_first;
              cross_link_candidate.beta = peptide_second;
              cross_link_candidate.cross_link_position.first = mod_pos;
              cross_link_candidate.cross_link_position.second = link_pos_second[x];
              cross_link_candidate.term_spec_alpha = first_spec;
              cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;;
            }
            else
            {
              cross_link_candidate.alpha = peptide_second;
              cross_link_candidate.beta = peptide_first;
              cross_link_candidate.cross_link_position.first = link_pos_second[x];
              cross_link_candidate.cross_link_position.second = mod_pos;
              cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;;
              cross_link_candidate.term_spec_beta = first_spec;
            }
            // Cross-linker mass is only one of the mono-link masses, if there is no second position (second == -1), otherwise the normal linker mass
            if (link_pos_second[x] != -1)
            {
              cross_link_candidate.cross_linker_mass = cross_link_mass;

#pragma omp critical (cross_link_candidates_access)
              cross_link_candidates.push_back(cross_link_candidate);
            }
            else
            {
              for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
              {
                // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
                bool is_correct_monolink = false;
                if (abs(spectrum_precursor_vector[precursor_correction_positions[i]] - (peptide_first->getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error_vector[precursor_correction_positions[i]])
                {
                  is_correct_monolink = true;
                }
                if (is_correct_monolink)
                {
                  cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];

#pragma omp critical (cross_link_candidates_access)
                  cross_link_candidates.push_back(cross_link_candidate);
                }
              }
            }
          }
        }
      }
    } // end of parallelized for-loop
    return cross_link_candidates;
  }

  void OPXLHelper::buildFragmentAnnotations(std::vector<PeptideHit::PeakAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum)
  {
    if (theoretical_spectrum.empty() || experiment_spectrum.empty())
    {
      return;
    }
    PeakSpectrum::IntegerDataArray charges = theoretical_spectrum.getIntegerDataArrays()[0];
    PeakSpectrum::StringDataArray names = theoretical_spectrum.getStringDataArrays()[0];
    for (Size k = 0; k < matching.size(); ++k)
    {
      PeptideHit::PeakAnnotation frag_anno;
      frag_anno.mz = experiment_spectrum[matching[k].second].getMZ();
      frag_anno.intensity = experiment_spectrum[matching[k].second].getIntensity();

      frag_anno.charge = charges[matching[k].first];
      frag_anno.annotation = names[matching[k].first];
      frag_annotations.push_back(frag_anno);
    }
  }

  void OPXLHelper::buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy)
  {
    for (Size i = 0; i < top_csms_spectrum.size(); ++i)
    {
      PeptideIdentification peptide_id;

      const PeakSpectrum& spectrum_light = spectra[scan_index];
      double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();

      String xltype = "cross-link";
      SignedSize alpha_pos = top_csms_spectrum[i].cross_link.cross_link_position.first;
      SignedSize beta_pos = top_csms_spectrum[i].cross_link.cross_link_position.second;

      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::MONO)
      {
        xltype = "mono-link";
      }
      else if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::LOOP)
      {
        xltype = "loop-link";
      }

      PeptideHit ph_alpha, ph_beta;
      // Set monolink as a modification or add MetaValue for cross-link identity and mass
      AASequence seq_alpha = *top_csms_spectrum[i].cross_link.alpha;
      ResidueModification::TermSpecificity alpha_term_spec = top_csms_spectrum[i].cross_link.term_spec_alpha;
      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::MONO)
      {
        vector< String > mods;
        const String residue = seq_alpha[alpha_pos].getOneLetterCode();

#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
        OPENMS_LOG_DEBUG << "Searching mono-link for " << residue << " | " << alpha_pos << endl;
#endif
        ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001, residue, ResidueModification::ANYWHERE);

#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
        OPENMS_LOG_DEBUG << "number of modifications fitting the diff mass: " << mods.size() << endl;
#endif

        bool mod_set = false;
        if (mods.size() > 0) // If several mods have the same diff mass, try to resolve ambiguity by cross-linker name (e.g. DSS and BS3 are different reagents, but have the same result after the reaction)
        {
          for (Size s = 0; s < mods.size(); ++s)
          {
            if (mods[s].hasSubstring(top_csms_spectrum[i].cross_link.cross_linker_name))
            {
#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
              OPENMS_LOG_DEBUG << "applied modification: " << mods[s] << endl;
#endif
              seq_alpha.setModification(alpha_pos, mods[s]);
              mod_set = true;
              break;
            }
          }
        }
        else if (mods.size() == 0 && (alpha_pos == 0 || alpha_pos == static_cast<int>(seq_alpha.size())-1))
        {
#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
          OPENMS_LOG_DEBUG << "No residue specific mono-link found, searching for terminal mods..." << endl;
#endif
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001, "", alpha_term_spec);
          if (mods.size() > 0)
          {
            Size mod_index = 0;
            for (Size s = 0; s < mods.size(); ++s)
            {
              if (mods[s].hasSubstring(top_csms_spectrum[i].cross_link.cross_linker_name))
              {
                mod_index = s;
              }
            }
            if (alpha_term_spec == ResidueModification::N_TERM)
            {
#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
              OPENMS_LOG_DEBUG << "Setting N-term mono-link: " << mods[mod_index] << endl;
#endif
              seq_alpha.setNTerminalModification(mods[mod_index]);
            }
            else
            {
#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
              OPENMS_LOG_DEBUG << "Setting C-term mono-link: " << mods[mod_index] << endl;
#endif
              seq_alpha.setCTerminalModification(mods[mod_index]);
            }
            mod_set = true;
          }
        }

        if ( (mods.size() > 0) && (!mod_set) ) // If resolving by name did not work, use any with matching diff mass
        {
          seq_alpha.setModification(alpha_pos, mods[0]);
          mod_set = true;
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, mods[0]);
        }
        if (!mod_set) // If no equivalent mono-link exists in the UNIMOD or XLMOD databases, use the given name to construct a placeholder
        {
          String mod_name = String("unknown mono-link " + top_csms_spectrum[i].cross_link.cross_linker_name + " mass " + String(top_csms_spectrum[i].cross_link.cross_linker_mass));
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, mod_name);
        }
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, top_csms_spectrum[i].cross_link.cross_linker_name);
      }
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS, DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));

      String alpha_term = "ANYWHERE";
      if (alpha_term_spec == ResidueModification::N_TERM)
      {
        alpha_term = "N_TERM";
      }
      else if (alpha_term_spec == ResidueModification::C_TERM)
      {
        alpha_term = "C_TERM";
      }

      ResidueModification::TermSpecificity beta_term_spec = top_csms_spectrum[i].cross_link.term_spec_beta;
      String beta_term = "ANYWHERE";
      if (beta_term_spec == ResidueModification::N_TERM)
      {
        beta_term = "N_TERM";
      }
      else if (beta_term_spec == ResidueModification::C_TERM)
      {
        beta_term = "C_TERM";
      }

      vector<PeptideHit> phs;
      if (beta_pos >= 0)
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue(beta_pos));
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, DataValue("-"));
      }

      ph_alpha.setSequence(seq_alpha);
      ph_alpha.setCharge(precursor_charge);
      ph_alpha.setScore(top_csms_spectrum[i].score);
      ph_alpha.setRank(DataValue(i+1));

      ph_alpha.setMetaValue("xl_chain", "MS:1002509");  // donor (longer, heavier, alphabetically earlier)
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, DataValue(alpha_pos));
      ph_alpha.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectra[scan_index].getNativeID());
      ph_alpha.setMetaValue("spectrum_index", scan_index);
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, xltype);
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK, DataValue(i + 1));
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA, alpha_term);
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, beta_term);
      ph_alpha.setMetaValue(Constants::UserParam::PRECURSOR_ASSIGNMENT_CORRECTION, top_csms_spectrum[i].precursor_correction);

      if (scan_index_heavy != scan_index)
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT, spectra[scan_index_heavy].getRT());
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ, spectra[scan_index_heavy].getPrecursors()[0].getMZ());
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF, spectra[scan_index_heavy].getNativeID());
        ph_alpha.setMetaValue("spectrum_index_heavy", scan_index_heavy);
      }
      ph_alpha.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, top_csms_spectrum[i].precursor_error_ppm);

      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_SCORE, top_csms_spectrum[i].score); // important for Percolator feature set because the PeptideHit score might be overwritten by a q-value
      ph_alpha.setMetaValue("OpenPepXL:xquest_score", top_csms_spectrum[i].xquest_score);
      ph_alpha.setMetaValue("OpenPepXL:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
      ph_alpha.setMetaValue("OpenPepXL:xcorr common", top_csms_spectrum[i].xcorrc_max);
      ph_alpha.setMetaValue("OpenPepXL:match-odds", top_csms_spectrum[i].match_odds);
      ph_alpha.setMetaValue("OpenPepXL:intsum", top_csms_spectrum[i].int_sum);
      ph_alpha.setMetaValue("OpenPepXL:intsum_alpha", top_csms_spectrum[i].intsum_alpha);
      ph_alpha.setMetaValue("OpenPepXL:intsum_beta", top_csms_spectrum[i].intsum_beta);
      ph_alpha.setMetaValue("OpenPepXL:total_current", top_csms_spectrum[i].total_current);

      ph_alpha.setMetaValue("OpenPepXL:wTIC", top_csms_spectrum[i].wTIC);
      ph_alpha.setMetaValue("OpenPepXL:TIC", top_csms_spectrum[i].percTIC);
      ph_alpha.setMetaValue("OpenPepXL:prescore", top_csms_spectrum[i].pre_score);

      ph_alpha.setMetaValue("OpenPepXL:log_occupancy", top_csms_spectrum[i].log_occupancy);
      ph_alpha.setMetaValue("OpenPepXL:log_occupancy_alpha", top_csms_spectrum[i].log_occupancy_alpha);
      ph_alpha.setMetaValue("OpenPepXL:log_occupancy_beta", top_csms_spectrum[i].log_occupancy_beta);

      ph_alpha.setMetaValue("matched_xlink_alpha",top_csms_spectrum[i].matched_xlink_alpha);
      ph_alpha.setMetaValue("matched_xlink_beta",top_csms_spectrum[i].matched_xlink_beta);
      ph_alpha.setMetaValue("matched_linear_alpha",top_csms_spectrum[i].matched_linear_alpha);
      ph_alpha.setMetaValue("matched_linear_beta",top_csms_spectrum[i].matched_linear_beta);

      ph_alpha.setMetaValue("ppm_error_abs_sum_linear_alpha", top_csms_spectrum[i].ppm_error_abs_sum_linear_alpha);
      ph_alpha.setMetaValue("ppm_error_abs_sum_linear_beta", top_csms_spectrum[i].ppm_error_abs_sum_linear_beta);
      ph_alpha.setMetaValue("ppm_error_abs_sum_xlinks_alpha", top_csms_spectrum[i].ppm_error_abs_sum_xlinks_alpha);
      ph_alpha.setMetaValue("ppm_error_abs_sum_xlinks_beta", top_csms_spectrum[i].ppm_error_abs_sum_xlinks_beta);

      ph_alpha.setMetaValue("ppm_error_abs_sum_linear", top_csms_spectrum[i].ppm_error_abs_sum_linear);
      ph_alpha.setMetaValue("ppm_error_abs_sum_xlinks", top_csms_spectrum[i].ppm_error_abs_sum_xlinks);
      ph_alpha.setMetaValue("ppm_error_abs_sum_alpha", top_csms_spectrum[i].ppm_error_abs_sum_alpha);
      ph_alpha.setMetaValue("ppm_error_abs_sum_beta", top_csms_spectrum[i].ppm_error_abs_sum_beta);
      ph_alpha.setMetaValue("ppm_error_abs_sum", top_csms_spectrum[i].ppm_error_abs_sum);

      ph_alpha.setMetaValue("precursor_total_intensity", top_csms_spectrum[i].precursor_total_intensity);
      ph_alpha.setMetaValue("precursor_target_intensity", top_csms_spectrum[i].precursor_target_intensity);
      ph_alpha.setMetaValue("precursor_signal_proportion", top_csms_spectrum[i].precursor_signal_proportion);
      ph_alpha.setMetaValue("precursor_target_peak_count", top_csms_spectrum[i].precursor_target_peak_count);
      ph_alpha.setMetaValue("precursor_residual_peak_count", top_csms_spectrum[i].precursor_residual_peak_count);

      ph_alpha.setMetaValue("selected", "false");

      ph_alpha.setPeakAnnotations(top_csms_spectrum[i].frag_annotations);

#ifdef DEBUG_OPXLHELPER
#pragma omp critical (LOG_DEBUG_access)
      OPENMS_LOG_DEBUG << "Annotations of size " << ph_alpha.getPeakAnnotations().size() << endl;
#endif

      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::CROSS)
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, (*top_csms_spectrum[i].cross_link.beta).toString());

        ph_beta.setSequence(*top_csms_spectrum[i].cross_link.beta);
        ph_beta.setCharge(precursor_charge);
        ph_beta.setScore(top_csms_spectrum[i].score);

        ph_beta.setMetaValue("xl_chain", "MS:1002510"); // receiver
        ph_beta.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectra[scan_index].getNativeID());
        ph_beta.setMetaValue("spectrum_index", scan_index);

        phs.push_back(ph_alpha);
        phs.push_back(ph_beta);
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, "-");
        phs.push_back(ph_alpha);
      }

      peptide_id.setRT(spectrum_light.getRT());
      peptide_id.setMZ(precursor_mz);
      String specIDs;
      if (scan_index_heavy != scan_index)
      {
        specIDs = spectra[scan_index].getNativeID() + "," + spectra[scan_index_heavy].getNativeID();
      }
      else
      {
        specIDs = spectra[scan_index].getNativeID();
      }

      peptide_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, specIDs);
      peptide_id.setHits(phs);
      peptide_id.setScoreType(Constants::UserParam::OPENPEPXL_SCORE);

// This critical section is called this way, because access to all_top_csms also happens in OpenPepXLAlgorithm and OpenPepXLLFAlgorithm.
// Access to peptide_ids is also critical, but it is only accessed here during parallel processing.
#ifdef _OPENMP
#pragma omp critical (all_top_csms_access)
#endif
      {
        peptide_ids.push_back(peptide_id);
        all_top_csms[all_top_csms_current_index][i].peptide_id_index = peptide_ids.size()-1;
      }
    }
  }

  void OPXLHelper::addProteinPositionMetaValues(std::vector< PeptideIdentification > & peptide_ids)
  {
    for (PeptideIdentification& id : peptide_ids)
    {
      if (id.getHits().empty()) continue;

      PeptideHit& ph_alpha = id.getHits()[0];
      String prot1_pos;

      // cross-link position in Protein (alpha)
      const std::vector<PeptideEvidence> pevs = ph_alpha.getPeptideEvidences();
      for (std::vector<PeptideEvidence>::const_iterator pev = pevs.begin(); pev != pevs.end(); ++pev)
      {
        // start counting at 1: pev->getStart() and xl_pos are both starting at 0,  with + 1 the N-term residue is number 1
        Int prot_link_pos = pev->getStart() + String(ph_alpha.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)).toInt() + 1;
        prot1_pos = prot1_pos + "," + prot_link_pos;
      }
      // remove leading "," of first position
      if (!prot1_pos.empty())
      {
        prot1_pos = prot1_pos.suffix(prot1_pos.size()-1);
      }
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT, prot1_pos);

      // cross-link position in Protein (beta)
      if (id.getHits().size() == 2)
      {
        PeptideHit& ph_beta = id.getHits()[1];
        String prot2_pos;
        String prot2_accessions;

        const std::vector<PeptideEvidence> pevs_beta = ph_beta.getPeptideEvidences();
        for (std::vector<PeptideEvidence>::const_iterator pev = pevs_beta.begin(); pev != pevs_beta.end(); ++pev)
        {
          // start counting at 1: pev->getStart() and xl_pos are both starting at 0,  with + 1 the N-term residue is number 1
          Int prot_link_pos = pev->getStart() + String(ph_alpha.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)).toInt() + 1;
          prot2_pos = prot2_pos + "," + prot_link_pos;
          prot2_accessions = prot2_accessions + "," + pev->getProteinAccession();
        }
        // remove leading "," of first position
        if (!prot2_pos.empty())
        {
          prot2_pos = prot2_pos.suffix(prot2_pos.size()-1);
        }
        ph_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1_PROT, prot1_pos);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT, prot2_pos);
        ph_beta.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT, prot2_pos);
      }
      else
      {
        // second cross-link position in Protein (loop-links)
        if (ph_alpha.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2) != "-")
        {
          String prot2_pos;
          for (std::vector<PeptideEvidence>::const_iterator pev = pevs.begin(); pev != pevs.end(); ++pev)
          {
            // start counting at 1: pev->getStart() and xl_pos are both starting at 0,  with + 1 the N-term residue is number 1
            Int prot_link_pos = pev->getStart() + String(ph_alpha.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)).toInt() + 1;
            prot2_pos = prot2_pos + "," + prot_link_pos;
          }
          // remove leading "," of first position
          if (!prot2_pos.empty())
          {
            prot2_pos = prot2_pos.suffix(prot2_pos.size()-1);
          }
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT, prot2_pos);
        }
        else
        {
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2_PROT, "-");
        }
      }
    }
  }

  void OPXLHelper::addBetaAccessions(std::vector< PeptideIdentification > & peptide_ids)
  {
    for (PeptideIdentification& id : peptide_ids)
    {
      if (id.getHits().empty()) continue;

      PeptideHit& ph_alpha = id.getHits()[0];

      if (id.getHits().size() == 2)
      {
        PeptideHit& ph_beta = id.getHits()[1];
        String prot2_accessions;

        const std::vector<PeptideEvidence> pevs_beta = ph_beta.getPeptideEvidences();
        for (std::vector<PeptideEvidence>::const_iterator pev = pevs_beta.begin(); pev != pevs_beta.end(); ++pev)
        {
          prot2_accessions = prot2_accessions + "," + pev->getProteinAccession();
        }

        if (!prot2_accessions.empty())
        {
          prot2_accessions = prot2_accessions.suffix(prot2_accessions.size()-1);
        }
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, prot2_accessions);
        ph_beta.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, prot2_accessions);
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_ACCESSIONS, "-");
      }
    }
  }

  void OPXLHelper::addXLTargetDecoyMV(std::vector< PeptideIdentification > & peptide_ids)
  {
    for (PeptideIdentification& id : peptide_ids)
    {
      if (id.getHits().empty()) continue;

      PeptideHit& ph_alpha = id.getHits()[0];
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_ALPHA, ph_alpha.getMetaValue(Constants::UserParam::TARGET_DECOY));

      // cross-link position in Protein (beta)
      if (id.getHits().size() == 2)
      {
        PeptideHit& ph_beta = id.getHits()[1];
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA, ph_beta.getMetaValue(Constants::UserParam::TARGET_DECOY));

        // if at least one of the two accession lists only contains decoys, the cross-link will be treated as a decoy
        if ( (!String(ph_alpha.getMetaValue(Constants::UserParam::TARGET_DECOY)).hasSubstring("target")) ||
             (!String(ph_beta.getMetaValue(Constants::UserParam::TARGET_DECOY)).hasSubstring("target")) )
        {
          ph_alpha.setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");
        }
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_TARGET_DECOY_BETA, "-");
      }
    }
  }

  void OPXLHelper::removeBetaPeptideHits(std::vector< PeptideIdentification > & peptide_ids)
  {
    // add PeptideEvidence data from Beta peptides as MetaValues to the Alpha peptide and remove the Beta PeptideHit
    for (PeptideIdentification& id : peptide_ids)
    {
      if (id.getHits().size() == 2 && id.getHits()[1].getMetaValue("xl_chain") == "MS:1002510")
      {

        const std::vector<PeptideEvidence>& peptide_evidences = id.getHits()[1].getPeptideEvidences();
        String pre, post, start, end;
        for (const PeptideEvidence& pe : peptide_evidences)
        {
          if (!pre.empty())
          {
            pre += ",";
            post += ",";
            start += ",";
            end += ",";
          }
          pre += pe.getAABefore();
          post += pe.getAAAfter();
          start += pe.getStart();
          end += pe.getEnd();
        }

        std::vector<PeptideHit> hits;
        PeptideHit hit = id.getHits()[0];
        hit.removeMetaValue("xl_chain");
        hit.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_PRE, pre);
        hit.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_POST, post);
        hit.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_START, start);
        hit.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_PEPEV_END, end);
        hits.push_back(hit);
        id.setHits(hits);
      }
    }

    // Collect PeptideHits to the same spectrum under one PeptideIdentification
    map<String, PeptideIdentification> new_peptide_ids;
    for (PeptideIdentification& id : peptide_ids)
    {
      if (id.getHits().size() > 0)
      {
        PeptideHit& hit = id.getHits()[0];
        PeptideIdentification new_id;
        String current_spectrum = id.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE);
        if (new_peptide_ids.find(current_spectrum) != new_peptide_ids.end())
        {
          new_id = (*new_peptide_ids.find(current_spectrum)).second;
        }
        else
        {
          new_id.setRT(id.getRT());
          new_id.setMZ(id.getMZ());
          new_id.setScoreType(Constants::UserParam::OPENPEPXL_SCORE);
          new_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, current_spectrum);
        }
        hit.removeMetaValue("xl_chain");
        new_id.insertHit(hit);
        new_peptide_ids[current_spectrum] = new_id;

      }
    }
    vector<PeptideIdentification> new_peptide_ids_vector;
    for (pair<String, PeptideIdentification> mit : new_peptide_ids)
    {
      new_peptide_ids_vector.push_back(mit.second);
    }
    peptide_ids = new_peptide_ids_vector;
  }

  void OPXLHelper::addPercolatorFeatureList(ProteinIdentification& prot_id)
  {
    // add features for percolator
    /* default features added in PercolatorAdapter:
    * SpecId, ScanNr, ExpMass, CalcMass, mass,
    * peplen, charge#min..#max, enzN, enzC, enzInt, dm, absdm
    */
    StringList feature_set;
    feature_set
      << Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM
      << Constants::UserParam::OPENPEPXL_SCORE
      << Constants::UserParam::PRECURSOR_ASSIGNMENT_CORRECTION
      << "OpenPepXL:xquest_score"
      << "OpenPepXL:xcorr xlink"
      << "OpenPepXL:xcorr common"
      << "OpenPepXL:match-odds"
      << "OpenPepXL:intsum"
      << "OpenPepXL:wTIC"
      << "OpenPepXL:TIC"
      << "OpenPepXL:prescore"
      << "OpenPepXL:log_occupancy"
      << "OpenPepXL:log_occupancy_alpha"
      << "OpenPepXL:log_occupancy_beta"
      << "matched_xlink_alpha"
      << "matched_xlink_beta"
      << "matched_linear_alpha"
      << "matched_linear_beta"
      << "ppm_error_abs_sum_linear_alpha"
      << "ppm_error_abs_sum_linear_beta"
      << "ppm_error_abs_sum_xlinks_alpha"
      << "ppm_error_abs_sum_xlinks_beta"
      << "ppm_error_abs_sum_linear"
      << "ppm_error_abs_sum_xlinks"
      << "ppm_error_abs_sum_alpha"
      << "ppm_error_abs_sum_beta"
      << "ppm_error_abs_sum"
      << "precursor_total_intensity"
      << "precursor_target_intensity"
      << "precursor_signal_proportion"
      << "precursor_target_peak_count"
      << "precursor_residual_peak_count";

    ProteinIdentification::SearchParameters search_params = prot_id.getSearchParameters();
    search_params.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    search_params.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    prot_id.setSearchParameters(search_params);
  }

  void OPXLHelper::computeDeltaScores(std::vector< PeptideIdentification > & peptide_ids)
  {
    for (PeptideIdentification& pep_id : peptide_ids)
    {
      pep_id.sort();
      vector<PeptideHit>& phs = pep_id.getHits();

      // at least two PeptideHits needed for an actual delta score
      if (phs.size() > 1)
      {
        for (Size rank = 0; rank < phs.size()-1; ++rank)
        {
          double delta_score = phs[rank+1].getScore() / phs[rank].getScore();
          phs[rank].setMetaValue(Constants::UserParam::DELTA_SCORE, delta_score);
        }
      }
      // set delta score to 0 for the last ranked PeptideHit, or if only one Peptide hit is available
      if (phs.size() > 0)
      {
        phs[phs.size()-1].setMetaValue(Constants::UserParam::DELTA_SCORE, 0.0);
      }
    }
  }

  std::vector< PeptideIdentification > OPXLHelper::combineTopRanksFromPairs(std::vector< PeptideIdentification > & peptide_ids, Size number_top_hits)
  {
    std::vector< PeptideIdentification > new_peptide_ids;
    std::vector< PeptideIdentification > current_spectrum_peptide_ids;
    std::set< String > spectrum_indices;

    for (PeptideIdentification& id : peptide_ids)
    {
      if (!id.getHits().empty())
      {
        spectrum_indices.insert(id.getHits()[0].getMetaValue("spectrum_index"));
      }
    }

    for (String index : spectrum_indices)
    {
      for (PeptideIdentification& id : peptide_ids)
      {
        if (!id.getHits().empty())
        {
          if (String(id.getHits()[0].getMetaValue("spectrum_index")) == index)
          {
            current_spectrum_peptide_ids.push_back(id);
          }
        }
      }

      std::sort(current_spectrum_peptide_ids.rbegin(), current_spectrum_peptide_ids.rend(), OPXLHelper::PeptideIDScoreComparator());
      current_spectrum_peptide_ids.erase( std::unique( current_spectrum_peptide_ids.begin(), current_spectrum_peptide_ids.end() ), current_spectrum_peptide_ids.end() );

      if (current_spectrum_peptide_ids.size() > number_top_hits)
      {
        current_spectrum_peptide_ids.resize(number_top_hits);
      }

      Size rank_count(1);
      for (PeptideIdentification& current_id : current_spectrum_peptide_ids)
      {
        if (!current_id.getHits().empty())
        {
          current_id.getHits()[0].setMetaValue("xl_rank", rank_count);
          if (current_id.getHits().size() == 2)
          {
            current_id.getHits()[1].setMetaValue("xl_rank", rank_count);
          }
        }
        rank_count++;
      }
      new_peptide_ids.insert(new_peptide_ids.end(), current_spectrum_peptide_ids.begin(), current_spectrum_peptide_ids.end());
      current_spectrum_peptide_ids.clear();
    }
    return new_peptide_ids;
  }

  std::vector <OPXLDataStructs::ProteinProteinCrossLink> OPXLHelper::collectPrecursorCandidates(const IntList& precursor_correction_steps, double precursor_mass, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, const vector<OPXLDataStructs::AASeqWithMass>& filtered_peptide_masses, double cross_link_mass, DoubleList cross_link_mass_mono_link, StringList cross_link_residue1, StringList cross_link_residue2, String cross_link_name)
  {
    // determine candidates
    std::vector< OPXLDataStructs::XLPrecursor > candidates;

    std::vector< double > spectrum_precursor_vector;
    std::vector< double > allowed_error_vector;

    for (int correction_mass : precursor_correction_steps)
    {
      double allowed_error = 0;

      double corrected_precursor_mass = precursor_mass - (static_cast<double>(correction_mass) * Constants::C13C12_MASSDIFF_U);

      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = corrected_precursor_mass * precursor_mass_tolerance * 1e-6;
      }
      else // Dalton
      {
        allowed_error = precursor_mass_tolerance;
      }

      spectrum_precursor_vector.push_back(corrected_precursor_mass);
      allowed_error_vector.push_back(allowed_error);

    } // end correction mass loop

    std::vector< int > precursor_correction_positions;
    candidates = OPXLHelper::enumerateCrossLinksAndMasses(filtered_peptide_masses, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, spectrum_precursor_vector, precursor_correction_positions, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

    vector< int > precursor_corrections;
    for (Size pc = 0; pc < precursor_correction_positions.size(); ++pc)
    {
      precursor_corrections.push_back(precursor_correction_steps[precursor_correction_positions[pc]]);
    }
    vector <OPXLDataStructs::ProteinProteinCrossLink> cross_link_candidates = OPXLHelper::buildCandidates(candidates, precursor_corrections, precursor_correction_positions, filtered_peptide_masses, cross_link_residue1, cross_link_residue2, cross_link_mass, cross_link_mass_mono_link, spectrum_precursor_vector, allowed_error_vector, cross_link_name);
    return cross_link_candidates;
  }

  double OPXLHelper::computePrecursorError(OPXLDataStructs::CrossLinkSpectrumMatch csm, double precursor_mz, int precursor_charge)
  {
    // Error calculation
    double weight = csm.cross_link.alpha->getMonoWeight();
    if (csm.cross_link.getType() == OPXLDataStructs::CROSS)
    {
      weight += csm.cross_link.beta->getMonoWeight() + csm.cross_link.cross_linker_mass;
    }
    else
    {
      weight += csm.cross_link.cross_linker_mass;
    }
    double precursor_mass = (precursor_mz * static_cast<double>(precursor_charge)) - (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)
                              - (static_cast<double>(csm.precursor_correction) * Constants::C13C12_MASSDIFF_U);
    double error = precursor_mass - weight;
    double rel_error = (error / precursor_mass) / 1e-6;
    return rel_error;
  }

  void OPXLHelper::isoPeakMeans(OPXLDataStructs::CrossLinkSpectrumMatch& csm, DataArrays::IntegerDataArray& num_iso_peaks_array, std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha, std::vector< std::pair< Size, Size > >& matched_spec_linear_beta, std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta)
  {
    csm.num_iso_peaks_mean = Math::mean(num_iso_peaks_array.begin(), num_iso_peaks_array.end());

    vector< double > iso_peaks_linear_alpha;
    vector< double > iso_peaks_linear_beta;
    vector< double > iso_peaks_xlinks_alpha;
    vector< double > iso_peaks_xlinks_beta;

    if (!matched_spec_linear_alpha.empty())
    {
      for (const auto& match : matched_spec_linear_alpha)
      {
        iso_peaks_linear_alpha.push_back(num_iso_peaks_array[match.second]);
      }
      csm.num_iso_peaks_mean_linear_alpha = Math::mean(iso_peaks_linear_alpha.begin(), iso_peaks_linear_alpha.end());
    }

    if (!matched_spec_linear_beta.empty())
    {
      for (const auto& match : matched_spec_linear_beta)
      {
        iso_peaks_linear_beta.push_back(num_iso_peaks_array[match.second]);
      }
      csm.num_iso_peaks_mean_linear_beta = Math::mean(iso_peaks_linear_beta.begin(), iso_peaks_linear_beta.end());
    }

    if (!matched_spec_xlinks_alpha.empty())
    {
      for (const auto& match : matched_spec_xlinks_alpha)
      {
        iso_peaks_xlinks_alpha.push_back(num_iso_peaks_array[match.second]);
      }
      csm.num_iso_peaks_mean_xlinks_alpha = Math::mean(iso_peaks_xlinks_alpha.begin(), iso_peaks_xlinks_alpha.end());
    }

    if (!matched_spec_xlinks_beta.empty())
    {
      for (const auto& match : matched_spec_xlinks_beta)
      {
        iso_peaks_xlinks_beta.push_back(num_iso_peaks_array[match.second]);
      }
      csm.num_iso_peaks_mean_xlinks_beta = Math::mean(iso_peaks_xlinks_beta.begin(), iso_peaks_xlinks_beta.end());
    }
  }

}
