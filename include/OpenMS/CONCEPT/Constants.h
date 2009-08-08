// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_CONSTANTS_H
#define OPENMS_CONCEPT_CONSTANTS_H

#include <OpenMS/config.h>

/**	
	@brief Main %OpenMS namespace.
		
	In this namespace all the main %OpenMS classes are located.
*/
namespace OpenMS 
{
	/**	
		@brief Mathematical and physical constants namespace.
			
		This namespace contains definitions for some basic mathematical and physical constants.
		All constants are double precision.
		<BR>
		There are basically two ways of accessing these constants:
		<UL>
			<LI> specify all namespaces:
			<BR>
			<tt>double my_pi = OpenMS::Constants::PI</tt>
			<BR>
			<LI>shortcut via the <tt>using directive</tt>:
			<BR>
			<tt>using namespace OpenMS::Constants;
			<BR>
			double my_pi = PI;</tt>
		</UL>
		
		@ingroup Concept
	*/
	namespace Constants 
	{
		/**	@name	Mathematical constants.
		*/
		//@{

		/// PI
		extern OPENMS_DLLAPI const double  PI;

		/// Euler's number - base of the natural logarithm
		extern OPENMS_DLLAPI const double  E;

		/**	Internal theshold for equality comparisons.
				Default value is 1e-6.
		*/
		extern double EPSILON;
		//@}
			
		/**	@name Chemical/physical constants.
		*/
		//@{
		
		/**	Elementary charge.
			  In units of C (\f$1.60217738 \cdot 10^{-19} C\f$).
		*/
		extern OPENMS_DLLAPI const double	ELEMENTARY_CHARGE;  	 // C     
	
		/// Elementary charge (alias)
		extern OPENMS_DLLAPI const double	e0;

		/** Electron mass.
				In units of kg (\f$9.1093897 \cdot 10^{-31}\f$ kg).
		*/
		extern OPENMS_DLLAPI const double	ELECTRON_MASS   	;   	 // kg

		/** Electron mass
		 		In units (\f$1,822.88850204(77)^{-1}\f$u).
		*/
		extern OPENMS_DLLAPI const double ELECTRON_MASS_U   ;      // u

		/** Proton mass.
				In units of kg (\f$1.6726230 \cdot 10^{-27}\f$ kg).
		*/
		extern OPENMS_DLLAPI const double	PROTON_MASS     	;   	 // kg

		/** Proton mass.
				In units (\f$1.00727646677(10)\f$u)
		*/
		extern OPENMS_DLLAPI const double PROTON_MASS_U     ;      // u

		/** Neutron mass.
				In units of kg (\f$1.6749286 \cdot 10^{-27}\f$ kg).
		*/
		extern OPENMS_DLLAPI const double	NEUTRON_MASS    	;   	 // kg

		/** Neutron mass.
		 		In units (\f$1.0086649156(6)\f$u)
		*/
		extern OPENMS_DLLAPI const double NEUTRON_MASS_U    ;      // u

		/** Avogadro constant.
				In units of \f$mol^{-1}\f$ (\f$6.0221367 \cdot 10^{23} mol^{-1}\f$).
		*/
		extern OPENMS_DLLAPI const double	AVOGADRO;

		/** Avogadro constant (alias)
		*/
		extern OPENMS_DLLAPI const double	NA;

		/** Avogadro constant (alias)
		*/
		extern OPENMS_DLLAPI const double	MOL;

		/** Boltzmann constant.
				In units of J/K (\f$1.380657 \cdot 10^{-23}\f$ J/K).
		*/
		extern OPENMS_DLLAPI const double	BOLTZMANN;

		/** Boltzmann constant (alias)
		*/
		extern OPENMS_DLLAPI const double	k;
		
		/** Planck constant.
				In units of Js (\f$6.6260754 \cdot 10^{-34}\f$ Js).
		*/
		extern OPENMS_DLLAPI const double	PLANCK;

		/** Planck constant (alias)
		*/
		extern OPENMS_DLLAPI const double	h;

		/** Gas constant (= NA * k)	
		*/
		extern OPENMS_DLLAPI const double	GAS_CONSTANT;

		/** Gas constant (alias)
		*/
		extern OPENMS_DLLAPI const double R;

		/** Faraday constant (= NA * e0)
		*/
		extern OPENMS_DLLAPI const double	FARADAY;

		/** Faraday constant (alias)
		*/
		extern OPENMS_DLLAPI const double	F;

		/** Bohr radius.
				In units m (\f$5.29177249 \cdot 10^{-11}\f$ m).
		*/
		extern OPENMS_DLLAPI const double	BOHR_RADIUS;

		/** Bohr radius (alias)
		*/
		extern OPENMS_DLLAPI const double	a0;

		//  the following values from: 
		//  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

		/** Vacuum permittivity.
				In units of \f$C^2J^{-1}m^{-1}\f$ (\f$8.85419 \cdot 10^{-12} C^2J^{-1}m^{-1}\f$).
		*/
		extern OPENMS_DLLAPI const double	VACUUM_PERMITTIVITY;

		/** Vacuum permeability.
				In units of \f$Js^2C^{-2}m^{-1}\f$ (\f$4\pi \cdot 10^{-7} Js^2C^{-2}m^{-1}\f$).
		*/
		extern OPENMS_DLLAPI const double	VACUUM_PERMEABILITY;

		/** Speed of light.
				In units of m/s (\f$2.99792458 \cdot 10^8 ms^{-1}\f$).
		*/
		extern OPENMS_DLLAPI const double	SPEED_OF_LIGHT;

		/** Speed of Light (alias)
		*/
		extern OPENMS_DLLAPI const double	c;

		/** Gravitational constant.
				In units of \f$Nm^2kg^{-2}\f$ (\f$6.67259 \cdot 10^{-11} Nm^2kg^{-2}\f$).
		*/
		extern OPENMS_DLLAPI const double	GRAVITATIONAL_CONSTANT;

		/** Fine structure constant.
				Without unit (\f$7.29735 \cdot 10^{-3}\f$).
		*/
		extern OPENMS_DLLAPI const double	FINE_STRUCTURE_CONSTANT;
		//@}

		/**	@name	Conversion factors
		*/
		//@{		
			
		/** Degree per rad.
				57.2957795130823209
		*/
		extern OPENMS_DLLAPI const double	DEG_PER_RAD;

		/** Rad per degree.
				0.0174532925199432957
		*/
		extern OPENMS_DLLAPI const double	RAD_PER_DEG;

		/** mm per inch.
				25.4
		*/
		extern OPENMS_DLLAPI const double	MM_PER_INCH 			;

		/** m per foot.
				3.048
		*/
		extern OPENMS_DLLAPI const double	M_PER_FOOT  			;

		/** Joules per calorie.
				4.184
		*/
		extern OPENMS_DLLAPI const double	JOULE_PER_CAL;

		/** Calories per Joule.
				1/JOULE_PER_CAL
		*/
		extern OPENMS_DLLAPI const double	CAL_PER_JOULE;

		//@}
	}
}

#endif // OPENMS_CONCEPT_CONSTANTS_H
