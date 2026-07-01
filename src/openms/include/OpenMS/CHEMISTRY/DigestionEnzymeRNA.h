// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

namespace OpenMS
{
  /**
     @brief Representation of a digestion enzyme for RNA (RNase)

     The cutting sites of these enzymes are defined using two different mechanisms:
     First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
     This is the same mechanism that is used for proteases (@see @ref ProteaseDigestion).
     However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
     Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
     If both expressions match, then there is a cutting site between the two ribonucleotides.

     There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
     A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI DigestionEnzymeRNA: public DigestionEnzyme
  {
  public:
    /// sets the "cuts after ..." regular expression
    void setCutsAfterRegEx(const std::string& value);

    /// returns the "cuts after ..." regular expression
    std::string getCutsAfterRegEx() const;

    /// sets the "cuts before ..." regular expression
    void setCutsBeforeRegEx(const std::string& value);

    /// returns the "cuts before ..." regular expression
    std::string getCutsBeforeRegEx() const;

    /// sets the 3' gain (as a nucleotide modification code)
    void setThreePrimeGain(const std::string& value);

    /// returns the 3' gain (as a nucleotide modification code)
    std::string getThreePrimeGain() const;

    /// sets the 5' gain (as a nucleotide modification code)
    void setFivePrimeGain(const std::string& value);

    /// returns the 5' gain (as a nucleotide modification code)
    std::string getFivePrimeGain() const;

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    bool setValueFromFile(const std::string& key, const std::string& value) override;

  protected:
    std::string three_prime_gain_;
    std::string five_prime_gain_;
    std::string cuts_after_regex_;
    std::string cuts_before_regex_;
  };

  typedef DigestionEnzymeRNA RNase;
}

