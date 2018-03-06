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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PRECURSOR_H
#define OPENMS_METADATA_PRECURSOR_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/Constants.h>

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

    /// Method of activation
    enum ActivationMethod
    {
      CID,                      ///< Collision-induced dissociation
      PSD,                      ///< Post-source decay
      PD,                       ///< Plasma desorption
      SID,                      ///< Surface-induced dissociation
      BIRD,                             ///< Blackbody infrared radiative dissociation
      ECD,                              ///< Electron capture dissociation
      IMD,                              ///< Infrared multiphoton dissociation
      SORI,                             ///< Sustained off-resonance irradiation
      HCID,                             ///< High-energy collision-induced dissociation
      LCID,                             ///< Low-energy collision-induced dissociation
      PHD,                              ///< Photodissociation
      ETD,                              ///< Electron transfer dissociation
      PQD,                              ///< Pulsed q dissociation
      SIZE_OF_ACTIVATIONMETHOD
    };

    /// Names of activation methods
    static const std::string NamesOfActivationMethod[SIZE_OF_ACTIVATIONMETHOD];
    static const std::string NamesOfActivationMethodShort[SIZE_OF_ACTIVATIONMETHOD];

    /// Constructor
    Precursor();
    /// Copy constructor
    Precursor(const Precursor & source);
    /// Destructor
    ~Precursor() override;

    /// Assignment operator
    Precursor & operator=(const Precursor & source);

    /// Equality operator
    bool operator==(const Precursor & rhs) const;
    /// Equality operator
    bool operator!=(const Precursor & rhs) const;

    /// returns a const reference to the activation methods
    const std::set<ActivationMethod> & getActivationMethods() const;
    /// returns a mutable reference to the activation methods
    std::set<ActivationMethod> & getActivationMethods();
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

    /// Non-mutable access to the charge
    Int getCharge() const;
    /// Mutable access to the charge
    void setCharge(Int charge);

    ///Mutable access to possible charge states
    std::vector<Int> & getPossibleChargeStates();
    ///Non-mutable access to possible charge states
    const std::vector<Int> & getPossibleChargeStates() const;
    ///Sets the possible charge states
    void setPossibleChargeStates(const std::vector<Int> & possible_charge_states);

    /// Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
    inline double getUnchargedMass() const
    {
      int c = charge_;
      (c == 0) ? c = 2 : c = charge_;
      return getMZ() * c - c * Constants::PROTON_MASS_U;
    }

protected:

    std::set<ActivationMethod> activation_methods_;
    double activation_energy_;
    double window_low_;
    double window_up_;
    double drift_time_;
    Int charge_;
    std::vector<Int> possible_charge_states_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_PRECURSOR_H
