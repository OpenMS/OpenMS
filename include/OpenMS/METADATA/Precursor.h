// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PRECURSOR_H
#define OPENMS_METADATA_PRECURSOR_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
namespace OpenMS 
{
	/**
		@brief Precursor meta information.
		
		This class stores precursor meta information, that is not already 
		covered by DSpectrum::getPrecursorPeak().
		
		@ingroup Metadata
	*/  
  class Precursor: public MetaInfoInterface
  {
    public:
    	/// Method of activation
      enum ActivationMethod
      {
      	ACTMETHNULL,	///< Unknown activation method
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

      /// Energy unit
      enum EnergyUnits
      {
      	UNITSNULL,	///< Unknown energy unit
      	EV,					///< Electron volt
      	PERCENT,		///< Percent
      	SIZE_OF_ENERGYUNITS
      };
			/// Names of energy units
			static const std::string NamesOfEnergyUnits[SIZE_OF_ENERGYUNITS];
      
      /// Constructor
      Precursor();
      /// Copy constructor
      Precursor(const Precursor& source);
      /// Destructor
      ~Precursor();
      
      /// Assignment operator
      Precursor& operator= (const Precursor& source);

      /// Equality operator
      bool operator== (const Precursor& rhs) const;
      /// Equality operator
      bool operator!= (const Precursor& rhs) const;
			
			/// returns the activation method
      ActivationMethod getActivationMethod() const;
      /// sets the activation method
      void setActivationMethod(ActivationMethod activation_method);
			
			/// returns the activation energy
      float getActivationEnergy() const;
      /// sets the activation energy
      void setActivationEnergy(float activation_energy);
			
			/// return the actication energy unit
      EnergyUnits getActivationEnergyUnit() const;
      /// sets the activation energy unit
      void setActivationEnergyUnit(EnergyUnits activation_energy_unit);
      
      /// returns the window size
      float getWindowSize() const;
      /// sets the window size
      void setWindowSize(float size);
 
    protected:
      ActivationMethod activation_method_;
      float activation_energy_;
      EnergyUnits activation_energy_unit_;
      float window_size_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_PRECURSOR_H
