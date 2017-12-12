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

#ifndef OPENMS_ANALYSIS_XLMS_OPXLHELPER_H
#define OPENMS_ANALYSIS_XLMS_OPXLHELPER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <numeric>

namespace OpenMS
{
  /**
   * @brief The OPXLHelper class contains functions needed by OpenPepXL and OpenPepXLLF to reduce duplicated code
   */
  class OPENMS_DLLAPI OPXLHelper
  {
    public:
    /**
       * @brief Enumerates precursor masses for all candidates in an XL-MS search
       * @param peptides The peptides with precomputed masses from the digestDatabase function
       * @param cross_link_mass_light mass of the cross-linker, only the light one if a labeled linker is used
       * @param cross_link_mass_mono_link A list of possible masses for the cross-link, if it is attached to a peptide on one side
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param spectrum_precursors A vector of all MS2 precursor masses of the searched spectra. Used to filter out candidates.
       * @param precursor_mass_tolerance The precursor mass tolerance
       * @param precursor_mass_tolerance_unit_ppm The unit of the precursor mass tolerance ("Da" or "ppm")
       * @return A vector of XLPrecursors containing all possible candidate cross-links
       */
      static std::vector<OPXLDataStructs::XLPrecursor> enumerateCrossLinksAndMasses(const std::vector<OPXLDataStructs::AASeqWithMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, std::vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);

      /**
       * @brief A helper function, that turns a StringList with modification names into a vector of ResidueModifications
       * @param modNames The list of modification names
       * @return A vector of modifications
       */
      static std::vector<ResidueModification> getModificationsFromStringList(StringList modNames);

      /**
       * @brief Digests a database with the given EnzymaticDigestion settings and precomputes masses for all peptides

          Also keeps track of the peptides at protein terminals and builds peptide candidates with all possible modification patterns
          according to the parameters.

       * @param fasta_db The protein database
       * @param digestor The object containing the digestion settings, e.g. the enzyme
       * @param min_peptide_length The minimal peptide length for the digestion
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param fixed_modifications The list of fixed modifications
       * @param variable_modifications The list of variable modifications
       * @param max_variable_mods_per_peptide The maximal number of variable modifications per peptide
       * @param count_proteins A variable to keep track of the number of proteins in the database. Should be an externally declared variable and = 0 when calling this function.
       * @param count_peptides A variable to keep track of the number of peptides after digestion. Should be an externally declared variable and = 0 when calling this function.
       * @param n_term_linker True, if the cross-linker can react with the N-terminal of a protein
       * @param c_term_linker True, if the cross-linker can react with the C-terminal of a protein
       * @return A vector of AASeqWithMass containing the peptides, their masses and information about terminal peptides
       */
      static std::vector<OPXLDataStructs::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins = 0, Size count_peptides = 0, bool n_term_linker = false, bool c_term_linker = false);

      /**
       * @brief Builds specific cross-link candidates with all possible combinations of linked positions from peptide pairs. Used to build candidates for the precursor mass window of a single MS2 spectrum.
       * @param candidates The XLPrecursors containing indices of two peptides
       * @param peptide_masses The digested peptide database, that the indices in the XLPrecursors refer to
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param cross_link_mass mass of the cross-linker, only the light one if a labeled linker is used
       * @param cross_link_mass_mono_link A list of possible masses for the cross-link, if it is attached to a peptide on one side
       * @param precursor_mass The precursor mass of the experimental spectrum (used to filter out certain candidates, e.g. mono- and loop-links have a different mass)
       * @param allowed_error The maximal precursor mass error in Da
       * @param cross_link_name The name of the cross-linker
       * @param n_term_linker True, if the cross-linker can react with the N-terminal of a protein
       * @param c_term_linker True, if the cross-linker can react with the C-terminal of a protein
       * @return A vector of ProteinProteinCrossLink candidates containing all necessary information to generate theoretical spectra
       */
      static std::vector <OPXLDataStructs::ProteinProteinCrossLink> buildCandidates(const std::vector< OPXLDataStructs::XLPrecursor > & candidates, const std::vector<OPXLDataStructs::AASeqWithMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker);

      /**
       * @brief Fills up the given FragmentAnnotation vector with annotations from a theoretical spectrum4

          This function takes an alignment of a theoretical spectrum with meta information and an experimental spectrum
          and builds annotations taking the MZ and intensity values from the experimental spectrum and the ion names and charges from the theoretical spectrum
          to annotate matched experimental peaks.
       * @param frag_annotations The vector to fill. Does not have to be empty, as annotations from several alignments can just be added on to the same vector.
       * @param matching The alignment between the two spectra
       * @param theoretical_spectrum The theoretical spectrum with meta information
       * @param experiment_spectrum The experimental specrum
       */
      static void buildFragmentAnnotations(std::vector<PeptideHit::PeakAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum);

      /**
       * @brief Builds PeptideIdentifications and PeptideHits
       * @param peptide_ids The vector of PeptideIdentifications for the whole experiment. The created PepIds will be pushed on this one.
       * @param top_csms_spectrum All CrossLinkSpectrumMatches from the current spectrum to be written out
       * @param all_top_csms A vector of all CrossLinkSpectrumMatches of the experiment, that is also extended in this function
       * @param all_top_csms_current_index The index of the current spectrum in all_top_csms (some spectra have no matches, so this is not equal to the spectrum index)
       * @param spectra The searched spectra as a PeakMap
       * @param scan_index The index of the current spectrum
       * @param scan_index_heavy The index of the heavy spectrum in a spectrum pair with labeled linkers
       */
      static void buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy);
  };
}

#endif // OPXLHELPER_H
