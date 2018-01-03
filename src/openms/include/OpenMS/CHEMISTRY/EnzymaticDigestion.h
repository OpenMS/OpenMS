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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H
#define OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

#include <string>
#include <vector>
#include <functional>

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
      SPEC_FULL, ///< fully enzyme specific, e.g., tryptic (ends with KR, AA-before is KR), or peptide is at protein terminal ends
      SPEC_SEMI, ///< semi specific, i.e., one of the two cleavage sites must fulfill requirements
      SPEC_NONE, ///< no requirements on start / end
      SIZE_OF_SPECIFICITY
    };
    /// Names of the Specificity
    static const std::string NamesOfSpecificity[SIZE_OF_SPECIFICITY];

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
    void setEnzyme(const DigestionEnzyme* enzyme);

    /// Returns the specificity for the digestion
    Specificity getSpecificity() const;

    /// Sets the specificity for the digestion (default is SPEC_FULL).
    void setSpecificity(Specificity spec);

    /// convert spec string name to enum
    /// returns SIZE_OF_SPECIFICITY if @p name is not valid
    static Specificity getSpecificityByName(const String& name);

    /**
       @brief Performs the enzymatic digestion of an unmodified sequence.

       By returning only references into the original string this is very fast.
       @p max_length restricts the maximum length of reported peptides (0 = no restriction)
    */
    void digestUnmodified(const StringView& sequence, std::vector<StringView>& output, Size min_length = 1, Size max_length = 0) const;

    /// Returns true if the fragment starting at position @p pos with length @p length within the sequence @p sequence could be generated by the current model
    bool isValidProduct(const String& sequence, Size pos, Size length, bool ignore_missed_cleavages = true) const;

    /**
       @brief Filter based on the number of missed cleavages.

       @param sequence Unmodified (!) amino acid sequence to check.
       @param filter A predicate that takes as parameter the number of missed cleavages in the sequence and returns true if the sequence should be filtered out.
       @return Whether the sequence should be filtered out.
     */
    bool filterByMissedCleavages(const String& sequence, std::function<bool(const Int)> filter) const;

protected:
    /// Returns the naive cleavage site positions without specificity (including '0' as first position, but not size() as last)
    std::vector<Size> tokenize_(const String& sequence) const;

    /**
       @brief Helper function for digestUnmodified()

       This function implements digestUnmodified() starting from the result of tokenize_().
       The separation enables derived classes to modify the result of tokenize_() during the in-silico digestion.
    */
    void digestAfterTokenize_(const std::vector<Size>& fragment_positions, const StringView& sequence, std::vector<StringView>& output, Size min_length = 0, Size max_length = -1) const;

    /**
       @brief Counts the number of missed cleavages in a sequence fragment

       @param cleavage_positions Positons of cleavage in protein as obtained from tokenize_()
       @param seq_start Index into sequence
       @param seq_end Past-the-end index into sequence
       @return number of missed cleavages of peptide
    */
    inline Size countMissedCleavages_(const std::vector<Size>& cleavage_positions, Size seq_start, Size seq_end) const;

    /// Number of missed cleavages
    Size missed_cleavages_;

    /// Used enzyme
    const DigestionEnzyme* enzyme_;

    /// specificity of enzyme
    Specificity specificity_;
  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H

