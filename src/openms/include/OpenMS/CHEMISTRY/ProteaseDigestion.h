// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <vector>

namespace OpenMS
{
  /**
     @brief Class for the enzymatic digestion of proteins represented as AASequence or String

     Digestion can be performed using simple regular expressions,
     e.g. [KR] | [^P]
     for trypsin. Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
     due to enzyme malfunction/access restrictions. If @em n missed cleavages are allowed, all possible resulting
     peptides (cleaved and uncleaved) with up to @em n missed cleavages are returned.
     Thus @b no random selection of just @em n specific missed cleavage sites is performed.

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProteaseDigestion: public EnzymaticDigestion
  {
  public:
    using EnzymaticDigestion::setEnzyme;

    /// Sets the enzyme for the digestion (by name)
    void setEnzyme(const String& name);

    /** 
       @brief Performs the enzymatic digestion of a protein represented as AASequence

       @param protein Sequence to digest
       @param output Digestion products (peptides)
       @param min_length Minimal length of reported products
       @param max_length Maximal length of reported products (0 = no restriction)
       @return Number of discarded digestion products (which are not matching length restrictions)

    */
    Size digest(const AASequence& protein, std::vector<AASequence>& output, Size min_length = 1, Size max_length = 0) const;

    /** 
       @brief Performs the enzymatic digestion of a protein represented as AASequence

       @param protein Sequence to digest
       @param output Digestion products (start and end indices of peptides)
       @param min_length Minimal length of reported products
       @param max_length Maximal length of reported products (0 = no restriction)
       @return Number of discarded digestion products (which are not matching length restrictions)

    */
    Size digest(const AASequence& protein, std::vector<std::pair<size_t,size_t>>& output, Size min_length = 1, Size max_length = 0) const;

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

