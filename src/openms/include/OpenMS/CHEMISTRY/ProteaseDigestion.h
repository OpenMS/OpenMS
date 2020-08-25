// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
     @brief Class for the enzymatic digestion of proteins

     Digestion can be performed using simple regular expressions,
     e.g. [KR] | [^P]
     for trypsin. Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
     due to enzyme malfunction/access restrictions. If @em n missed cleavages are allowed, all possible resulting
     peptides (cleaved and uncleaved) with up to @em n missed cleavages are returned.
     Thus @b no random selection of just @em n specific missed cleavage sites is performed.

     An alternative model is also available in EnzymaticDigestionLogModel.

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProteaseDigestion: public EnzymaticDigestion
  {
  public:
    using EnzymaticDigestion::setEnzyme;

    /// Sets the enzyme for the digestion (by name)
    void setEnzyme(const String& name);

    /** 
       @brief: Performs the enzymatic digestion of a protein.

       @param protein Sequence to digest
       @param output Digestion products (peptides)
       @param min_length Minimal length of reported products
       @param max_length Maximal length of reported products (0 = no restriction)
       @return Number of discarded digestion products (which are not matching length restrictions)

    */
    Size digest(const AASequence& protein, std::vector<AASequence>& output, Size min_length = 1, Size max_length = 0) const;

    /// Returns the number of peptides a digestion of @p protein would yield under the current enzyme and missed cleavage settings.
    Size peptideCount(const AASequence& protein);

    /**
      @brief Variant of EnzymaticDigestion::isValidProduct() with support for n-term protein cleavage and random D|P cleavage

      Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the flags provided here.

      @param protein Protein sequence
      @param pep_pos Starting index of potential peptide
      @param pep_length Length of potential peptide
      @param ignore_missed_cleavages Do not compare MC's of potential peptide to the maximum allowed MC's
      @param allow_nterm_protein_cleavage Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 and protein starts with 'M'
      @param allow_random_asp_pro_cleavage Allow cleavage at D|P sites to count as n/c-terminal.
      @return True if peptide has correct n/c terminals (according to enzyme, specificity and above flags)

    */
    bool isValidProduct(const String& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages = true, bool allow_nterm_protein_cleavage = false, bool allow_random_asp_pro_cleavage = false) const;

    /// forwards to isValidProduct using protein.toUnmodifiedString()
    bool isValidProduct(const AASequence& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages = true, bool allow_nterm_protein_cleavage = false, bool allow_random_asp_pro_cleavage = false) const;

  };

} // namespace OpenMS

