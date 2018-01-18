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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_TAGGING_H
#define OPENMS_METADATA_TAGGING_H

#include <OpenMS/METADATA/Modification.h>

namespace OpenMS
{
  /**
      @brief Meta information about tagging of a sample e.g. ICAT labeling.

      Holds information about the mass difference between light and heavy tag.
      All other relevant information is provided by Modification.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Tagging :
    public Modification
  {
public:
    /// Isotope variants (light and heavy)
    enum IsotopeVariant {LIGHT, HEAVY, SIZE_OF_ISOTOPEVARIANT};
    /// Names of isotope variants
    static const std::string NamesOfIsotopeVariant[SIZE_OF_ISOTOPEVARIANT];

    /// default constructor
    Tagging();
    /// copy constructor
    Tagging(const Tagging &);
    /// destructor
    ~Tagging() override;

    /// assignment operator
    Tagging & operator=(const Tagging &);

    /**
        @brief Equality operator

    Although this operator takes a reference to a SampleTreatment as argument
    it tests for the equality of Tagging instances!
  */
    bool operator==(const SampleTreatment & rhs) const override;

    /// clone method. See SampleTreatment
    SampleTreatment * clone() const override;

    /// returns the mass difference between light and heavy variant (default is 0.0)
    double getMassShift() const;
    /// sets the mass difference between light and heavy variant
    void setMassShift(double mass_shift);

    /// returns the isotope variant of the tag (default is LIGHT)
    const IsotopeVariant & getVariant() const;
    /// sets the isotope variant of the tag
    void setVariant(const IsotopeVariant & variant);

protected:
    double mass_shift_;
    IsotopeVariant variant_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_TAGGING_H
