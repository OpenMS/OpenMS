// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Representation of an isotope

      This class extends Element to store additional isotope-specific properties
      including half-life and decay mode. Isotopes are specific nuclear configurations
      of an element characterized by their mass number (protons + neutrons).
  */
  class OPENMS_DLLAPI Isotope : public Element
  {
public:

    /// Enum for decay modes
    enum class DecayMode
    {
      STABLE,       ///< Stable isotope (no decay)
      UNKNOWN,      ///< Decay mode unknown
      ALPHA,        ///< Alpha decay
      BETA_PLUS,    ///< Beta plus decay (positron emission)
      BETA_MINUS,   ///< Beta minus decay (electron emission)
      PROTON,       ///< Proton emission
      SIZE_OF_DECAYMODE
    };

    /** @name Constructor and Destructors
    */
    //@{
    /// default constructor
    Isotope();

    /// copy constructor
    Isotope(const Isotope& isotope);

    /// detailed constructor
    Isotope(const std::string& name,
            const std::string& symbol,
            unsigned int atomic_number,
            unsigned int mass_number,
            double atomic_mass,
            double abundance,
            double half_life = -1.0,
            DecayMode decay_mode = DecayMode::STABLE);

    /// destructor
    ~Isotope() override;
    //@}

    /** @name Accessors
    */
    //@{
    /// sets the mass number (protons + neutrons)
    void setMassNumber(unsigned int mass_number);

    /// returns the mass number (protons + neutrons)
    unsigned int getMassNumber() const;

    /// returns the number of neutrons (mass_number - atomic_number)
    unsigned int getNeutrons() const;

    /// sets the natural abundance (0.0 to 1.0)
    void setAbundance(double abundance);

    /// returns the natural abundance (0.0 to 1.0)
    double getAbundance() const;

    /// sets the half-life in seconds (negative value means stable)
    void setHalfLife(double half_life);

    /// returns the half-life in seconds (negative value means stable)
    double getHalfLife() const;

    /// sets the decay mode
    void setDecayMode(DecayMode decay_mode);

    /// returns the decay mode
    DecayMode getDecayMode() const;

    /// returns whether the isotope is stable (half_life < 0)
    bool isStable() const;

    /// returns whether this is an isotope (always true for Isotope class)
    virtual bool isIsotope() const;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    Isotope& operator=(const Isotope& isotope);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const Isotope& isotope) const;

    /// inequality operator
    bool operator!=(const Isotope& isotope) const;

    /// less operator (compares by mass number)
    bool operator<(const Isotope& isotope) const;
    //@}

    /// writes the isotope to an output stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Isotope& isotope);

    /// Convert DecayMode enum to string
    static std::string decayModeToString(DecayMode mode);

    /// Convert string to DecayMode enum
    static DecayMode stringToDecayMode(const std::string& mode_str);

protected:

    /// mass number (protons + neutrons)
    unsigned int mass_number_;

    /// natural abundance (0.0 to 1.0)
    double abundance_;

    /// half-life in seconds (negative means stable)
    double half_life_;

    /// decay mode
    DecayMode decay_mode_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream&, const Isotope&);

} // namespace OpenMS
