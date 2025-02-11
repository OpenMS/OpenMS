// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <set>

namespace OpenMS
{
  /**
      @brief Precursor meta information.

      This class contains precursor information:

        - isolation window
        - activation
        - selected ion (m/z, intensity, charge, possible charge states)
        - ion mobility drift time

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Precursor :
    public CVTermList,
    public Peak1D
  {

public:
    /// Constructor
    Precursor() = default;
    /// Copy constructor
    Precursor(const Precursor&) = default;

    // note: we implement the move constructor ourselves due to a bug in MSVS
    // 2015/2017 which cannot produce a default move constructor for classes
    // that contain STL containers (other than vector).

    /// Move constructor
    Precursor(Precursor&&) noexcept;
    /// Destructor
    ~Precursor() override = default;

    /// Assignment operator
    Precursor& operator=(const Precursor&) = default;
    /// Move assignment operator
    Precursor& operator=(Precursor&&) & = default;

    /// Method of activation
    enum ActivationMethod
    {
      CID,                      ///< Collision-induced dissociation (MS:1000133) (also CAD; parent term, but unless otherwise stated often used as synonym for trap-type CID)
      PSD,                      ///< Post-source decay
      PD,                       ///< Plasma desorption
      SID,                      ///< Surface-induced dissociation
      BIRD,                     ///< Blackbody infrared radiative dissociation
      ECD,                      ///< Electron capture dissociation (MS:1000250)
      IMD,                      ///< Infrared multiphoton dissociation
      SORI,                     ///< Sustained off-resonance irradiation
      HCID,                     ///< High-energy collision-induced dissociation
      LCID,                     ///< Low-energy collision-induced dissociation
      PHD,                      ///< Photodissociation
      ETD,                      ///< Electron transfer dissociation
      ETciD,                    ///< Electron transfer and collision-induced dissociation (MS:1003182)
      EThcD,                    ///< Electron transfer and higher-energy collision dissociation (MS:1002631) 
      PQD,                      ///< Pulsed q dissociation (MS:1000599)
      TRAP,                     ///< trap-type collision-induced dissociation (MS:1002472)
      HCD,                      ///< beam-type collision-induced dissociation (MS:1000422)
      INSOURCE,                 ///< in-source collision-induced dissociation (MS:1001880)
      LIFT,                     ///< Bruker proprietary method (MS:1002000)
      SIZE_OF_ACTIVATIONMETHOD
    };
    /// Names of activation methods
    static const std::string NamesOfActivationMethod[SIZE_OF_ACTIVATIONMETHOD];
    static const std::string NamesOfActivationMethodShort[SIZE_OF_ACTIVATIONMETHOD];

    /// Equality operator
    bool operator==(const Precursor & rhs) const;
    /// Equality operator
    bool operator!=(const Precursor & rhs) const;
  
    /// returns a const reference to the activation methods
    const std::set<ActivationMethod>& getActivationMethods() const;
    /// returns a mutable reference to the activation methods
    std::set<ActivationMethod>& getActivationMethods();
    /// convenience function, returning string representation of getActivationMethods()
    StringList getActivationMethodsAsString() const;    
    /// sets the activation methods
    void setActivationMethods(const std::set<ActivationMethod> & activation_methods);

    /// returns the activation energy (in electronvolt)
    double getActivationEnergy() const;
    /// sets the activation energy (in electronvolt)
    void setActivationEnergy(double activation_energy);

    /**
     * @brief Returns the lower offset from the target m/z
     *
     * @note This is an offset relative to the target m/z. The start of the
     * mass isolation window should thus be computed as:
     *
     *   p.getMZ() - p.getIsolationWindowLowerOffset()
     *
     * @return the lower offset from the target m/z
     */
    double getIsolationWindowLowerOffset() const;
    /// sets the lower offset from the target m/z
    void setIsolationWindowLowerOffset(double bound);

    /**
     * @brief Returns the upper offset from the target m/z
     *
     * @note This is an offset relative to the target m/z. The end of the mass
     * isolation window should thus be computed as:
     *
     *   p.getMZ() + p.getIsolationWindowUpperOffset()
     *
     * @return the upper offset from the target m/z
     */
    double getIsolationWindowUpperOffset() const;
    /// sets the upper offset from the target m/z
    void setIsolationWindowUpperOffset(double bound);

    /**
      @brief Returns the ion mobility drift time in milliseconds (-1 means it is not set)

      @note It is possible for the spectrum to not have a Precursor but still
      have a drift time, please check getDriftTime of MSSpectrum first and only
      use this function if you need find-grained access to individual precursors.
    */
    double getDriftTime() const;
    /// sets the ion mobility drift time in milliseconds
    void setDriftTime(double drift_time);

    /**
      @brief Returns the ion mobility drift time unit
    */
    DriftTimeUnit getDriftTimeUnit() const;

    /**
      @brief Sets the ion mobility drift time unit
    */
    void setDriftTimeUnit(DriftTimeUnit dt);


    /**
     * @brief Returns the lower offset from the target ion mobility in milliseconds
     *
     * @note This is an offset relative to the target ion mobility. The start
     * of the ion mobility isolation window should thus be computed as:
     *
     *   p.getDriftTime() + p.getDriftTimeWindowLowerOffset()
     *
     * @return the lower offset from the target ion mobility
    */
    double getDriftTimeWindowLowerOffset() const;
    /// sets the lower offset from the target ion mobility
    void setDriftTimeWindowLowerOffset(double drift_time);

    /**
     * @brief Returns the upper offset from the target ion mobility in milliseconds
     *
     * @note This is an offset relative to the target ion mobility. The end
     * of the ion mobility isolation window should thus be computed as:
     *
     *   p.getDriftTime() + p.getDriftTimeWindowUpperOffset()
     *
     * @return the upper offset from the target ion mobility
    */
    double getDriftTimeWindowUpperOffset() const;
    /// sets the upper offset from the target ion mobility
    void setDriftTimeWindowUpperOffset(double drift_time);

    /// Non-mutable access to the charge
    Int getCharge() const;
    /// Mutable access to the charge
    void setCharge(Int charge);

    /// Mutable access to possible charge states
    std::vector<Int>& getPossibleChargeStates();
    /// Non-mutable access to possible charge states
    const std::vector<Int>& getPossibleChargeStates() const;
    /// Sets the possible charge states
    void setPossibleChargeStates(const std::vector<Int> & possible_charge_states);

    /// Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0, our best guess is doubly charged
    inline double getUnchargedMass() const
    {
      int c = charge_;
      (c == 0) ? c = 2 : c = charge_;
      return getMZ() * c - c * Constants::PROTON_MASS_U;
    }

protected:

    std::set<ActivationMethod> activation_methods_;
    double activation_energy_{};
    double window_low_{};
    double window_up_{};
    double drift_time_{-1};
    double drift_window_low_{};
    double drift_window_up_{};
    DriftTimeUnit drift_time_unit_{DriftTimeUnit::NONE};
    Int charge_{};
    std::vector<Int> possible_charge_states_;
  };
} // namespace OpenMS

