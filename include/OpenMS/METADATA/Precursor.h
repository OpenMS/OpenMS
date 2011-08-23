// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

		@ingroup Metadata
	*/
  class OPENMS_DLLAPI Precursor
    : public CVTermList,
			public Peak1D
  {

    public:

    /// Method of activation
      enum ActivationMethod
      {
      	CID,					///< Collision-induced dissociation
      	PSD,					///< Post-source decay
      	PD,						///< Plasma desorption
      	SID,					///< Surface-induced dissociation
				BIRD,					///< Blackbody infrared radiative dissociation
				ECD,					///< Electron capture dissociation
				IMD,					///< Infrared multiphoton dissociation
				SORI,					///< Sustained off-resonance irradiation
				HCID,					///< High-energy collision-induced dissociation
				LCID,					///< Low-energy collision-induced dissociation
				PHD,					///< Photodissociation
				ETD,					///< Electron transfer dissociation
				PQD,					///< Pulsed q dissociation
      	SIZE_OF_ACTIVATIONMETHOD
      };
			/// Names of activation methods
			static const std::string NamesOfActivationMethod[SIZE_OF_ACTIVATIONMETHOD];

      /// Constructor
      Precursor();
      /// Copy constructor
      Precursor(const Precursor& source);
      /// Destructor
      virtual ~Precursor();

      /// Assignment operator
      Precursor& operator= (const Precursor& source);

      /// Equality operator
      bool operator== (const Precursor& rhs) const;
      /// Equality operator
      bool operator!= (const Precursor& rhs) const;

			/// returns a const reference to the activation methods
      const std::set<ActivationMethod>& getActivationMethods() const;
			/// returns a mutable reference to the activation methods
      std::set<ActivationMethod>& getActivationMethods();
      /// sets the activation methods
      void setActivationMethods(const std::set<ActivationMethod>& activation_methods);

			/// returns the activation energy (in electronvolt)
      DoubleReal getActivationEnergy() const;
      /// sets the activation energy (in electronvolt)
      void setActivationEnergy(DoubleReal activation_energy);

      /// returns the lower offset from the target m/z
      DoubleReal getIsolationWindowLowerOffset() const;
      /// sets the lower offset from the target m/z
      void setIsolationWindowLowerOffset(DoubleReal bound);

      /// returns the upper offset from the target m/z
      DoubleReal getIsolationWindowUpperOffset() const;
      /// sets the upper offset from the target m/z
      void setIsolationWindowUpperOffset(DoubleReal bound);

      /// Non-mutable access to the charge
      Int getCharge() const;
      /// Mutable access to the charge
      void setCharge( Int charge );

      ///Mutable access to possible charge states
      std::vector<Int>& getPossibleChargeStates();
      ///Non-mutable access to possible charge states
      const std::vector<Int>& getPossibleChargeStates() const;
      ///Sets the possible charge states
      void setPossibleChargeStates(const std::vector<Int>& possible_charge_states);

	/// Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
	inline DoubleReal getUnchargedMass() const
	{
		int c = charge_;
		(c==0)?c=2:c=charge_;
		return (getMZ()*c - c*Constants::PROTON_MASS_U);
	}

    protected:

      std::set<ActivationMethod> activation_methods_;
      DoubleReal activation_energy_;
      DoubleReal window_low_;
      DoubleReal window_up_;
      Int charge_;
      std::vector<Int> possible_charge_states_;
    };
} // namespace OpenMS

#endif // OPENMS_METADATA_PRECURSOR_H
