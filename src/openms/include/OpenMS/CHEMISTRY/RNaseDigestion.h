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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <boost/regex.hpp>

namespace OpenMS
{
  /**
     @brief Class for the enzymatic digestion of RNAs

     @see @ref DigestionEnzymeRNA

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI RNaseDigestion: public EnzymaticDigestion
  {
  public:
    /// Sets the enzyme for the digestion
    void setEnzyme(const DigestionEnzyme* enzyme) override;

    /// Sets the enzyme for the digestion (by name)
    void setEnzyme(const String& name);

    /**
       @brief Performs the enzymatic digestion of a (potentially modified) RNA

       Only fragments of appropriate length (between @p min_length and @p max_length) are returned.
    */
    void digest(const NASequence& rna, std::vector<NASequence>& output,
                Size min_length = 0, Size max_length = 0) const;

    /**
       @brief Performs the enzymatic digestion of all RNA parent molecules in @p IdentificationData

       Digestion products are stored as IdentifiedOligos with corresponding MoleculeParentMatch annotations.
       Only fragments of appropriate length (between @p min_length and @p max_length) are included.
    */
    void digest(IdentificationData& id_data, Size min_length = 0,
                Size max_length = 0) const;

  protected:
    const Ribonucleotide* five_prime_gain_; ///< 5' mod added by the enzyme
    const Ribonucleotide* three_prime_gain_; ///< 3' mod added by the enzyme
    boost::regex cuts_after_regex_; ///< reg. exp. for enzyme cutting pattern
    boost::regex cuts_before_regex_; ///< reg. exp. for enzyme cutting pattern

    /**
       @brief Returns the positions of digestion products in the RNA as pairs: (start, length)
     */
    std::vector<std::pair<Size, Size>> getFragmentPositions_(
      const NASequence& rna, Size min_length, Size max_length)
      const;
  };

} // namespace OpenMS

