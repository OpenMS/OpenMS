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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
//#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/FORMAT/FASTAFile.h>


using namespace std;

namespace OpenMS
{
  // check whether the candidate pair is within the given tolerance to at least one precursor mass in the spectra data
  void filter_and_add_candidate (vector<OPXLDataStructs::XLPrecursor>& mass_to_candidates, vector< double >& spectrum_precursors, bool precursor_mass_tolerance_unit_ppm, double precursor_mass_tolerance, OPXLDataStructs::XLPrecursor precursor)
  {
    bool found_matching_precursors = false;

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
      found_matching_precursors = true;
    }


    // if precursors were found in the above for-loop, add candidate to results vector
    if (found_matching_precursors)
    {
// don't access this vector from two processing threads at the same time
#ifdef _OPENMP
#pragma omp critical
#endif
      mass_to_candidates.push_back(precursor);
    }
  }


  // Enumerate all pairs of peptides from the searched database and calculate their masses (inlcuding mono-links and loop-links)
  vector<OPXLDataStructs::XLPrecursor> OPXLHelper::enumerateCrossLinksAndMasses(const vector<OPXLDataStructs::AASeqWithMass>&  peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
  {
    // initialize empty vector for the results
    vector<OPXLDataStructs::XLPrecursor> mass_to_candidates;
    // initialize progress counter
    Size countA = 0;

    double min_precursor = spectrum_precursors[0];
    double max_precursor = spectrum_precursors[spectrum_precursors.size()-1];

// Multithreading options: schedule: static, dynamic, guided
// use OpenMP to run this for-loop on multiple CPU cores
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize p1 = 0; p1 < static_cast<SignedSize>(peptides.size()); ++p1)
    {
      // get the amino acid sequence of this peptide as a character string
      String seq_first = peptides[p1].peptide_seq.toUnmodifiedString();

      // every 500 peptides print current progress to console
      countA += 1;
      if (countA % 500 == 0)
      {
        cout << "Enumerating pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << " | current size in mb: " << mass_to_candidates.size() * sizeof(OPXLDataStructs::XLPrecursor) / 1024 / 1024 << endl;
      }

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
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }


       // test if this peptide could have loop-links: one cross-link with both sides attached to the same peptide
       // TODO check for distance between the two linked residues
      bool first_res = false; // is there a residue the first side of the linker can attach to?
      bool second_res = false; // is there a residue the second side of the linker can attach to?
      for (Size k = 0; k < seq_first.size()-1; ++k)
      {
        for (Size i = 0; i < cross_link_residue1.size(); ++i)
        {
          if (seq_first.substr(k, 1) == cross_link_residue1[i])
          {
            first_res = true;
          }
        }
        for (Size i = 0; i < cross_link_residue2.size(); ++i)
        {
          if (seq_first.substr(k, 1) == cross_link_residue2[i])
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
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
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
        } else if (peptides[p2].peptide_mass > max_second_peptide_mass)
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
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }
    }
    cout << "Enumerated pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << " | current size in mb: " << mass_to_candidates.size() * sizeof(OPXLDataStructs::XLPrecursor) / 1024 / 1024 << endl;
    return mass_to_candidates;
  }

  vector<ResidueModification> OPXLHelper::getModificationsFromStringList(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  std::vector<OPXLDataStructs::AASeqWithMass> OPXLHelper::digestDatabase(vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins, Size count_peptides, bool n_term_linker, bool c_term_linker)
  {
    multimap<StringView, AASequence> processed_peptides;
    vector<OPXLDataStructs::AASeqWithMass> peptide_masses;
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {
//#ifdef _OPENMP
//#pragma omp atomic
//#endif
      ++count_proteins;

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
        } else if (fasta_db[fasta_index].sequence.hasSuffix(cit->getString()))
        {
          position = OPXLDataStructs::C_TERM;
        }

        // skip if no cross-linked residue
        bool skip = true;
        for (Size k = 0; k < cross_link_residue1.size(); k++)
        {
          if (cit->getString().find(cross_link_residue1[k]) < cit->getString().size()-1)
          {
            skip = false;
          }
          if (n_term_linker && position == OPXLDataStructs::N_TERM)
          {
            skip = false;
          }
          if (c_term_linker && position == OPXLDataStructs::C_TERM)
          {
            skip = false;
          }
        }
        for (Size k = 0; k < cross_link_residue2.size(); k++)
        {
          if (cit->getString().find(cross_link_residue2[k]) < cit->getString().size()-1)
          {
            skip = false;
          }
          if (n_term_linker && position == OPXLDataStructs::N_TERM)
          {
            skip = false;
          }
          if (c_term_linker && position == OPXLDataStructs::C_TERM)
          {
            skip = false;
          }
        }
        if (skip) continue;

        bool already_processed = false;
//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
        {
          if (processed_peptides.find(*cit) != processed_peptides.end())
          {
            // peptide (and all modified variants) already processed so skip it
            already_processed = true;
          }
        }

        if (already_processed)
        {
          continue;
        }
//        if (cit->getString().find('K') >= cit->getString().size()-1)
//        {
//          continue;
//        }



//#ifdef _OPENMP
//#pragma omp atomic
//#endif
        ++count_peptides;

        vector<AASequence> all_modified_peptides;

        // generate all modified variants of a peptide
        // Note: no critial section is needed despite ResidueDB not beeing thread sage.
        //       It is only written to on introduction of novel modified residues. These residues have been already added above (single thread context).
        {
          AASequence aas = AASequence::fromString(cit->getString());
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < static_cast<SignedSize>(all_modified_peptides.size()); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          OPXLDataStructs::AASeqWithMass pep_mass;
          pep_mass.peptide_mass = candidate.getMonoWeight();
          pep_mass.peptide_seq = candidate;
          pep_mass.position = position;

//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
          {
            processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
            peptide_masses.push_back(pep_mass);
          }
        }
      }
    }
    sort(peptide_masses.begin(), peptide_masses.end(), OPXLDataStructs::AASeqWithMassComparator());
    return peptide_masses;
  }

  vector <OPXLDataStructs::ProteinProteinCrossLink> OPXLHelper::buildCandidates(const std::vector< OPXLDataStructs::XLPrecursor > & candidates, const std::vector<OPXLDataStructs::AASeqWithMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker)
  {
    vector <OPXLDataStructs::ProteinProteinCrossLink> cross_link_candidates;
    for (Size i = 0; i < candidates.size(); ++i)
    {
      OPXLDataStructs::XLPrecursor candidate = candidates[i];
      vector <SignedSize> link_pos_first;
      vector <SignedSize> link_pos_second;
      AASequence peptide_first = peptide_masses[candidate.alpha_index].peptide_seq;
      OPXLDataStructs::PeptidePosition peptide_pos_first = peptide_masses[candidate.alpha_index].position;
      AASequence peptide_second;
      OPXLDataStructs::PeptidePosition peptide_pos_second = OPXLDataStructs::INTERNAL;
      if (candidate.beta_index < peptide_masses.size())
      {
        peptide_second = peptide_masses[candidate.beta_index].peptide_seq;
        peptide_pos_second = peptide_masses[candidate.beta_index].position;
      }
      String seq_first = peptide_first.toUnmodifiedString();
      String seq_second =  peptide_second.toUnmodifiedString();

      // TODO mono-links and loop-links with different masses can be generated for the same precursor mass, but only one of them can be valid each time.
      // Find out which is the case. But it should not happen often enough to slow down the tool significantly.
      bool is_loop = abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass)) <= allowed_error;

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
      } else
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
      } else if (seq_second.size() == seq_first.size() && peptide_second.getMonoWeight() > peptide_first.getMonoWeight())
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
            cross_link_candidates.push_back(cross_link_candidate);
          }
          else
          {
            for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
            {
              // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
              if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
              {
                cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];;
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
        if (n_term_linker && (peptide_pos_second == OPXLDataStructs::N_TERM))
        {
//          second_spec = ResidueModification::N_TERM;
//          mod_pos = 0;
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_second == OPXLDataStructs::C_TERM))
        {
          second_spec = ResidueModification::C_TERM;
          mod_pos = peptide_second.size()-1;
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
            cross_link_candidates.push_back(cross_link_candidate);

          }
        }
      }

      if (peptide_pos_first != OPXLDataStructs::INTERNAL)
      {
        ResidueModification::TermSpecificity first_spec = ResidueModification::N_TERM;
        Size mod_pos = 0;
        bool compatible = false;
        if (n_term_linker && (peptide_pos_first == OPXLDataStructs::N_TERM))
        {
//          first_spec = ResidueModification::N_TERM;
//          mod_pos = 0;
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_first == OPXLDataStructs::C_TERM))
        {
          first_spec = ResidueModification::C_TERM;
          mod_pos = peptide_first.size()-1;
          compatible = true;
        }
        if (compatible)
        {
          for (Size x = 0; x < link_pos_second.size(); ++x)
          {
            OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate;
            cross_link_candidate.cross_linker_name = cross_link_name;
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
              cross_link_candidates.push_back(cross_link_candidate);
            }
            else
            {
              for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
              {
                // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
                if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
                {
                  cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];
                  cross_link_candidates.push_back(cross_link_candidate);
                }
              }
            }
          }
        }
      }
    }
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
      AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
      ResidueModification::TermSpecificity alpha_term_spec = top_csms_spectrum[i].cross_link.term_spec_alpha;
      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::MONO)
      {
        //AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
        vector< String > mods;
        const String residue = seq_alpha[alpha_pos].getOneLetterCode();
        LOG_DEBUG << "Searching mono-link for " << residue << " | " << alpha_pos << endl;
        ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001, residue, ResidueModification::ANYWHERE);
        LOG_DEBUG << "number of modifications fitting the diff mass: " << mods.size() << endl;
        bool mod_set = false;
        if (mods.size() > 0) // If several mods have the same diff mass, try to resolve ambiguity by cross-linker name (e.g. DSS and BS3 are different reagents, but have the same result after the reaction)
        {
          for (Size s = 0; s < mods.size(); ++s)
          {
            if (mods[s].hasSubstring(top_csms_spectrum[i].cross_link.cross_linker_name))
            {
              LOG_DEBUG << "applied modification: " << mods[s] << endl;
              seq_alpha.setModification(alpha_pos, mods[s]);
              mod_set = true;
              break;
            }
          }
        }
        else if (mods.size() == 0 && (alpha_pos == 0 || alpha_pos == static_cast<int>(seq_alpha.size())-1))
        {
          LOG_DEBUG << "No residue specific mono-link found, searching for terminal mods..." << endl;
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
//            cout << "Terminal Mod Test; mod: " << mod_index <<  " | " << mods[mod_index] << " | term_spec: " << alpha_term_spec << endl;
            if (alpha_term_spec == ResidueModification::N_TERM)
            {
              LOG_DEBUG << "Setting N-term mono-link: " << mods[mod_index] << endl;
              seq_alpha.setNTerminalModification(mods[mod_index]);
            }
            else
            {
              LOG_DEBUG << "Setting C-term mono-link: " << mods[mod_index] << endl;
              seq_alpha.setCTerminalModification(mods[mod_index]);
            }
            mod_set = true;
          }
        }

        if ( (mods.size() > 0) && (!mod_set) ) // If resolving by name did not work, use any with matching diff mass
        {
          seq_alpha.setModification(alpha_pos, mods[0]);
          mod_set = true;
        }
        if (!mod_set) // If no equivalent mono-link exists in the UNIMOD or XLMOD databases, use the given name to construct a placeholder
        {
          String mod_name = String("unknown mono-link " + top_csms_spectrum[i].cross_link.cross_linker_name + " mass " + String(top_csms_spectrum[i].cross_link.cross_linker_mass));
          //seq_alpha.setModification(alpha_pos, mod_name);
          LOG_DEBUG << "unknown mono-link" << endl;
          ph_alpha.setMetaValue("xl_mod", mod_name);
          ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
        }
      }
      else
      {
        ph_alpha.setMetaValue("xl_mod", top_csms_spectrum[i].cross_link.cross_linker_name);
        ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
      }


      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::LOOP)
      {
        ph_alpha.setMetaValue("xl_pos2", DataValue(beta_pos));
      }

      // Error calculation
      double weight = seq_alpha.getMonoWeight() + top_csms_spectrum[i].cross_link.cross_linker_mass;
      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::CROSS)
      {
        weight += top_csms_spectrum[i].cross_link.beta.getMonoWeight();
      }
      double theo_mz = (weight + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / static_cast<double>(precursor_charge);
      double error = precursor_mz - theo_mz;
      double rel_error = (error / theo_mz) / 1e-6;

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

      ph_alpha.setSequence(seq_alpha);
      ph_alpha.setCharge(precursor_charge);
      ph_alpha.setScore(top_csms_spectrum[i].score);
      ph_alpha.setRank(DataValue(i+1));
      ph_alpha.setMetaValue("xl_chain", "MS:1002509");  // donor (longer, heavier, alphabetically earlier)
      ph_alpha.setMetaValue("xl_pos", DataValue(alpha_pos));
      ph_alpha.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
      ph_alpha.setMetaValue("xl_type", xltype);
      ph_alpha.setMetaValue("xl_rank", DataValue(i + 1));
      ph_alpha.setMetaValue("xl_term_spec", alpha_term);

      if (scan_index_heavy != scan_index)
      {
        ph_alpha.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
        ph_alpha.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
        ph_alpha.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
      }
      ph_alpha.setMetaValue(Constants::PRECURSOR_ERROR_PPM_USERPARAM, rel_error);

      ph_alpha.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
      ph_alpha.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
      ph_alpha.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
      ph_alpha.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
      ph_alpha.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

      ph_alpha.setMetaValue("OpenXQuest:log_occupancy", top_csms_spectrum[i].log_occupancy);
      ph_alpha.setMetaValue("OpenXQuest:log_occupancy_alpha", top_csms_spectrum[i].log_occupancy_alpha);
      ph_alpha.setMetaValue("OpenXQuest:log_occupancy_beta", top_csms_spectrum[i].log_occupancy_beta);
      ph_alpha.setMetaValue("OpenXQuest:log_occupancy_full_spec", top_csms_spectrum[i].log_occupancy_full_spec);

      ph_alpha.setMetaValue("OpenPepXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
      ph_alpha.setMetaValue("OpenPepXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
      ph_alpha.setMetaValue("OpenPepXL:HyperAlpha", top_csms_spectrum[i].HyperAlpha);
      ph_alpha.setMetaValue("OpenPepXL:HyperBeta", top_csms_spectrum[i].HyperBeta);
      ph_alpha.setMetaValue("OpenPepXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

      ph_alpha.setMetaValue("OpenPepXL:PScoreCommon",top_csms_spectrum[i].PScoreCommon);
      ph_alpha.setMetaValue("OpenPepXL:PScoreXlink",top_csms_spectrum[i].PScoreXlink);
      ph_alpha.setMetaValue("OpenPepXL:PScoreAlpha",top_csms_spectrum[i].PScoreAlpha);
      ph_alpha.setMetaValue("OpenPepXL:PScoreBeta",top_csms_spectrum[i].PScoreBeta);
      ph_alpha.setMetaValue("OpenPepXL:PScoreBoth",top_csms_spectrum[i].PScoreBoth);

      ph_alpha.setMetaValue("selected", "false");

      ph_alpha.setPeakAnnotations(top_csms_spectrum[i].frag_annotations);
      LOG_DEBUG << "Annotations of size " << ph_alpha.getPeakAnnotations().size() << endl;
      phs.push_back(ph_alpha);

      if (top_csms_spectrum[i].cross_link.getType() == OPXLDataStructs::CROSS)
      {
        ph_beta.setSequence(top_csms_spectrum[i].cross_link.beta);
        ph_beta.setCharge(precursor_charge);
        ph_beta.setScore(top_csms_spectrum[i].score);
        ph_beta.setRank(DataValue(i+1));
        ph_beta.setMetaValue("xl_chain", "MS:1002510"); // receiver
        ph_beta.setMetaValue("xl_pos", DataValue(beta_pos));
        ph_beta.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
        ph_beta.setMetaValue("xl_term_spec", beta_term);

        if (scan_index_heavy != scan_index)
        {
          ph_beta.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
          ph_beta.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
          ph_beta.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
        }
        ph_beta.setMetaValue(Constants::PRECURSOR_ERROR_PPM_USERPARAM, rel_error);

        ph_beta.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
        ph_beta.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
        ph_beta.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
        ph_beta.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
        ph_beta.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

        ph_beta.setMetaValue("OpenXQuest:log_occupancy", top_csms_spectrum[i].log_occupancy);
        ph_beta.setMetaValue("OpenXQuest:log_occupancy_alpha", top_csms_spectrum[i].log_occupancy_alpha);
        ph_beta.setMetaValue("OpenXQuest:log_occupancy_beta", top_csms_spectrum[i].log_occupancy_beta);
        ph_beta.setMetaValue("OpenXQuest:log_occupancy_full_spec", top_csms_spectrum[i].log_occupancy_full_spec);

        ph_beta.setMetaValue("OpenPepXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
        ph_beta.setMetaValue("OpenPepXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
        ph_beta.setMetaValue("OpenPepXL:HyperAlpha",top_csms_spectrum[i].HyperAlpha);
        ph_beta.setMetaValue("OpenPepXL:HyperBeta",top_csms_spectrum[i].HyperBeta);
        ph_beta.setMetaValue("OpenPepXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

        ph_beta.setMetaValue("OpenPepXL:PScoreCommon",top_csms_spectrum[i].PScoreCommon);
        ph_beta.setMetaValue("OpenPepXL:PScoreXlink",top_csms_spectrum[i].PScoreXlink);
        ph_beta.setMetaValue("OpenPepXL:PScoreAlpha",top_csms_spectrum[i].PScoreAlpha);
        ph_beta.setMetaValue("OpenPepXL:PScoreBeta",top_csms_spectrum[i].PScoreBeta);
        ph_beta.setMetaValue("OpenPepXL:PScoreBoth",top_csms_spectrum[i].PScoreBoth);

        ph_beta.setMetaValue("selected", "false");

        phs.push_back(ph_beta);
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

      peptide_id.setMetaValue("spectrum_reference", specIDs);
//      peptide_id.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
//      peptide_id.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
//      peptide_id.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
//      peptide_id.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
//      peptide_id.setMetaValue("xl_type", xltype); // TODO: needs CV term
//      peptide_id.setMetaValue("xl_rank", DataValue(i + 1));

      peptide_id.setHits(phs);
      peptide_id.setScoreType("OpenXQuest:combined score");

#ifdef _OPENMP
#pragma omp critical (peptides_ids_access)
#endif
      {
        peptide_ids.push_back(peptide_id);
        all_top_csms[all_top_csms_current_index][i].peptide_id_index = peptide_ids.size()-1;
      }
    }
  }

}
