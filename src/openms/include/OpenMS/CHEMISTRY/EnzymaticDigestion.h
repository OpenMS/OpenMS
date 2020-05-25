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
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

#include <boost/regex.hpp>
#include <string>
#include <vector>

#include <functional> // for std::function

namespace OpenMS
{
  /**
     @brief Class for the enzymatic digestion of sequences

     Digestion can be performed using simple regular expressions,
     e.g. [KR] | [^P]
     for trypsin. Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
     due to enzyme malfunction/access restrictions. If @em n missed cleavages are given, all possible resulting
     peptides (cleaved and uncleaved) with up to @em n missed cleavages are returned.
     Thus @b no random selection of just @em n specific missed cleavage sites is performed.

     @see ProteaseDigestion for functionality specific to protein digestion.

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI EnzymaticDigestion
  {
public:
    /// when querying for valid digestion products, this determines if the specificity of the two peptide ends is considered important
    enum Specificity
    {
      SPEC_NONE = 0, ///< no requirements on start / end
      SPEC_SEMI = 1, ///< semi specific, i.e., one of the two cleavage sites must fulfill requirements
      SPEC_FULL = 2, ///< fully enzyme specific, e.g., tryptic (ends with KR, AA-before is KR), or peptide is at protein terminal ends
      SPEC_UNKNOWN = 3,
      SPEC_NOCTERM = 8, ///< no requirements on CTerm (currently not supported in the class)
      SPEC_NONTERM = 9, ///< no requirements on NTerm (currently not supported in the class)
      SIZE_OF_SPECIFICITY = 10
    };
    /// Names of the Specificity
    static const std::string NamesOfSpecificity[SIZE_OF_SPECIFICITY];

    /// Name for no cleavage
    static const std::string NoCleavage;

    /// Name for unspecific cleavage
    static const std::string UnspecificCleavage;

    /// Default constructor
    EnzymaticDigestion();

    /// Destructor
    virtual ~EnzymaticDigestion();

    /// Returns the number of missed cleavages for the digestion
    Size getMissedCleavages() const;

    /// Sets the number of missed cleavages for the digestion (default is 0). This setting is ignored when log model is used.
    void setMissedCleavages(Size missed_cleavages);

    /// Returns the enzyme for the digestion
    String getEnzymeName() const;

    /// Sets the enzyme for the digestion
    virtual void setEnzyme(const DigestionEnzyme* enzyme);

    /// Returns the specificity for the digestion
    Specificity getSpecificity() const;

    /// Sets the specificity for the digestion (default is SPEC_FULL).
    void setSpecificity(Specificity spec);

    /// convert spec string name to enum
    /// returns SPEC_UNKNOWN if @p name is not valid
    static Specificity getSpecificityByName(const String& name);

    /**
       @brief Performs the enzymatic digestion of an unmodified sequence.

       By returning only references into the original string this is very fast.

       @param sequence Sequence to digest
       @param output Digestion products
       @param min_length Minimal length of reported products
       @param max_length Maximal length of reported products (0 = no restriction)
       @return Number of discarded digestion products (which are not matching length restrictions)
       */
    Size digestUnmodified(const StringView& sequence, std::vector<StringView>& output, Size min_length = 1, Size max_length = 0) const;

    /**
     @brief Performs the enzymatic digestion of an unmodified sequence.

     By returning only positions into the original string this is very fast and compared to the StringView output
     version of this function it is independent of the original sequence. Can be used for matching products to
     determine e.g. missing ones. @todo could be set of pairs.

     @param sequence Sequence to digest
     @param output Digestion products as vector of pairs of start and end positions
     @param min_length Minimal length of reported products
     @param max_length Maximal length of reported products (0 = no restriction)
     @return Number of discarded digestion products (which are not matching length restrictions)
     */
    Size digestUnmodified(const StringView& sequence, std::vector<std::pair<Size,Size>>& output, Size min_length = 1, Size max_length = 0) const;

    /**
    @brief Is the peptide fragment starting at position @p pos with length @p length within the sequence @p sequence generated by the current enzyme?

    Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here.

    @param protein Protein sequence
    @param pep_pos Starting index of potential peptide
    @param pep_length Length of potential peptide
    @param ignore_missed_cleavages Do not compare MC's of potential peptide to the maximum allowed MC's
    @return True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
    */
    bool isValidProduct(const String& sequence, int pos, int length, bool ignore_missed_cleavages = true) const;

    /**
       @brief Filter based on the number of missed cleavages.

       @param sequence Unmodified (!) amino acid sequence to check.
       @param filter A predicate that takes as parameter the number of missed cleavages in the sequence and returns true if the sequence should be filtered out.
       @return Whether the sequence should be filtered out.
     */
    bool filterByMissedCleavages(const String& sequence, const std::function<bool(const Int)>& filter) const;

protected:

    /**
      @brief supports functionality for ProteaseDigestion as well (which is deeply weaved into the function)
             To avoid code duplication, this is stored here and called by wrappers.
             Do not duplicate the code, just for the sake of semantics (unless we can come up with a clean separation)
             Note: the overhead of allow_nterm_protein_cleavage and allow_random_asp_pro_cleavage is marginal; the main runtime is spend during tokenize_()
    */
    bool isValidProduct_(const String& sequence,
                         int pos,
                         int length,
                         bool ignore_missed_cleavages,
                         bool allow_nterm_protein_cleavage,
                         bool allow_random_asp_pro_cleavage) const;
    /**
      @brief Digests the sequence using the enzyme's regular expression

      The resulting split positions include @p start as first position, but not end.
      If start is negative, it is reset to zero.
      If end is negative or beyond @p sequence's size(), it is set to size().
      All returned positions are relative to the full @p sequence.
      
      Returned positions include @p start and any positions between start and end matching the regex.      
            
      @param sequence ...
      @param start Start digestion after this point
      @param end Past-the-end index into @p sequence
      @return Cleavage positions (this includes @p start, but not @p end)
     */
    std::vector<int> tokenize_(const String& sequence, int start = 0, int end = -1) const;

    /**
       @brief Helper function for digestUnmodified()

       This function implements digestUnmodified() starting from the result of tokenize_().
       The separation enables derived classes to modify the result of tokenize_() during the in-silico digestion.

       @return number of digestion products NOT matching the length restrictions
    */
    Size digestAfterTokenize_(const std::vector<int>& fragment_positions, const StringView& sequence, std::vector<StringView>& output, Size min_length = 0, Size max_length = -1) const;
    Size digestAfterTokenize_(const std::vector<int>& fragment_positions, const StringView& sequence, std::vector<std::pair<Size,Size>>& output, Size min_length = 0, Size max_length = -1) const;

    /**
       @brief Counts the number of missed cleavages in a sequence fragment

       @param cleavage_positions Positions of cleavage in protein as obtained from tokenize_()
       @param seq_start Index into sequence
       @param seq_end Past-the-end index into sequence
       @return number of missed cleavages of peptide
    */
    Size countMissedCleavages_(const std::vector<int>& cleavage_positions, Size seq_start, Size seq_end) const;

    /// Number of missed cleavages
    Size missed_cleavages_;

    /// Used enzyme
    const DigestionEnzyme* enzyme_;
    /// Regex for tokenizing (huge speedup by making this a member instead of stack object in tokenize_())
    boost::regex re_;

    /// specificity of enzyme
    Specificity specificity_;
  };

} // namespace OpenMS


