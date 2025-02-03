// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

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

    /// Default constructor
    Tagging();
    /// Copy constructor
    Tagging(const Tagging &) = default;
    /// Move constructor
    Tagging(Tagging&&) = default;
    /// Destructor
    ~Tagging() override;

    /// Assignment operator
    Tagging & operator=(const Tagging &) = default;
    /// Move assignment operator
    Tagging& operator=(Tagging&&) & = default;

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

