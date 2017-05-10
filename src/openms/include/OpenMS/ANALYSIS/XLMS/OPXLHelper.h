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
#include <OpenMS/ANALYSIS/XLMS/OpenProXLUtils.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <numeric>

namespace OpenMS
{
  class OPENMS_DLLAPI OPXLHelper
  {
    public:
      static std::vector<OpenProXLUtils::XLPrecursor> enumerateCrossLinksAndMasses_(const std::vector<OpenProXLUtils::AASeqWithMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, std::vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);

      static std::vector<ResidueModification> getModificationsFromStringList(StringList modNames);

      static std::vector<OpenProXLUtils::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins = 0, Size count_peptides = 0, bool n_term_linker = false, bool c_term_linker = false);

      static std::vector <ProteinProteinCrossLink> buildCandidates(const std::vector< OpenProXLUtils::XLPrecursor > & candidates, const std::vector<OpenProXLUtils::AASeqWithMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker);

      static void buildFragmentAnnotations(std::vector<PeptideHit::FragmentAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum);

      static void buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy = NULL);
  };
}

#endif // OPXLHELPER_H
