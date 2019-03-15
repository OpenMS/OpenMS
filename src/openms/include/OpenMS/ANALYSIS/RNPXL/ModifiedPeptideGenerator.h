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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  class OPENMS_DLLAPI ModifiedPeptideGenerator
  {
   /*
    * @brief Modifications can be generated and applied to AASequences. 
    */

  public
:
    // Applies fixed modifications to a single peptide
    static void applyFixedModifications(const std::vector<ResidueModification>::const_iterator& fixed_mods_begin, const std::vector<ResidueModification>::const_iterator& fixed_mods_end, AASequence& peptide);

    // Applies variable modifications to a single peptide. If keep_original is set the original (e.g. unmodified version) is also returned
    static void applyVariableModifications(const std::vector<ResidueModification>::const_iterator& var_mods_begin, const std::vector<ResidueModification>::const_iterator& var_mods_end, const AASequence& peptide, Size max_variable_mods_per_peptide, std::vector<AASequence>& all_modified_peptides, bool keep_original=true);

  protected:
    // Recursively generate all combinatoric placements at compatible sites
    static void recurseAndGenerateVariableModifiedPeptides_(const std::vector<int>& subset_indices, const std::map<int, std::vector<ResidueModification> >& map_compatibility, int depth, const AASequence& current_peptide, std::vector<AASequence>& modified_peptides);

    // Fast implementation of modification placement. No combinatoric placement is needed in this case - just every site is modified once by each compatible modification. Already modified residues are skipped
    static void applyAtMostOneVariableModification_(const std::vector<ResidueModification>::const_iterator& var_mods_begin, const std::vector<ResidueModification>::const_iterator& var_mods_end, const AASequence& peptide, std::vector<AASequence>& all_modified_peptides, bool keep_original=true);

  };
}


