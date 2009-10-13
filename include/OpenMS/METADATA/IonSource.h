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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IONSOURCE_H
#define OPENMS_METADATA_IONSOURCE_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of a ion source (part of a MS Instrument)
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI IonSource
  	: public MetaInfoInterface
  {
    public:
    	/// inlet type
	    enum InletType
	    {
	    	INLETNULL,													///< Unknown
	    	DIRECT,															///< Direct
	    	BATCH,															///< Batch (e.g. in MALDI)
	    	CHROMATOGRAPHY,											///< Chromatography (liquid)
	    	PARTICLEBEAM,												///< Particle beam
	    	MEMBRANESEPARATOR,									///< Membrane sparator
	    	OPENSPLIT,													///< Open split
	    	JETSEPARATOR,												///< Jet separator
	    	SEPTUM,															///< Septum
	    	RESERVOIR,													///< Reservoir
	    	MOVINGBELT,													///< Moving belt
	    	MOVINGWIRE,													///< Moving wire
	    	FLOWINJECTIONANALYSIS,							///< Flow injection analysis
	    	ELECTROSPRAYINLET,									///< Electro spray
	    	THERMOSPRAYINLET,										///< Thermo spray
	    	INFUSION,														///< Infusion
	    	CONTINUOUSFLOWFASTATOMBOMBARDMENT,	///< Continuous flow fast atom bombardment
	    	INDUCTIVELYCOUPLEDPLASMA,						///< Inductively coupled plasma
	    	MEMBRANE,														///< Membrane inlet
	    	NANOSPRAY,													///< Nanospray inlet
	    	SIZE_OF_INLETTYPE
	    };
			/// Names of inlet types
			static const std::string NamesOfInletType[SIZE_OF_INLETTYPE];

	    /// ionization method
	    enum IonizationMethod
	    {
	    	IONMETHODNULL,	///< Unknown
	    	ESI,						///< electrospray ionisation
	    	EI,							///< electron ionization
	    	CI,							///< chemical ionisation
	    	FAB,						///< fast atom bombardment
	    	TSP,						///< thermospray
	    	LD,							///< laser desorption
	    	FD,							///< field desorption
	    	FI,							///< flame ionization
	    	PD,							///< plasma desorption
	    	SI,							///< secondary ion MS
	    	TI,							///< thermal ionization
	    	API,						///< atmospheric pressure ionisation
	    	ISI,						///<
	    	CID,						///< collsion induced decomposition
	    	CAD,						///< collsion activated decomposition
	    	HN,							///<
	    	APCI,						///< atmospheric pressure chemical ionization
	    	APPI,						///< atmospheric pressure photo ionization
	    	ICP,						///< inductively coupled plasma
	    	NESI,						///< Nano electrospray ionization
	    	MESI,						///< Micro electrospray ionization
	    	SELDI,					///< Surface enhanced laser desorption ionization
	    	SEND,						///< Surface enhanced neat desorption
	    	FIB,						///< Fast ion bombardment
	    	MALDI,					///< Matrix-assisted laser desorption ionization
	    	MPI,						///< Multiphoton ionization
	    	DI,							///< desorption ionization
	    	FA,							///< flowing afterglow
	    	FII,						///< field ionization
	    	GD_MS,					///< glow discharge ionization
	    	NICI,						///< negative ion chemical ionization
	    	NRMS,						///< neutralization reionization mass spectrometry
	    	PI,							///< photoionization
	    	PYMS,						///< pyrolysis mass spectrometry
	    	REMPI,					///< resonance enhanced multiphoton ionization
	    	AI,							///< adiabatic ionization
	    	ASI,						///< associative ionization
	    	AD,							///< autodetachment
	    	AUI,						///< autoionization
	    	CEI,						///< charge exchange ionization
	    	CHEMI,					///< chemi-ionization
	    	DISSI,					///< dissociative ionization
	    	LSI,						///< liquid secondary ionization
	    	PEI,						///< penning ionization
	    	SOI,						///< soft ionization
	    	SPI,						///< spark ionization
	    	SUI,						///< surface ionization
	    	VI,							///< vertical ionization
	    	AP_MALDI,				///< atmospheric pressure matrix-assisted laser desorption ionization
	    	SILI,						///< desorption/ionization on silicon
	    	SALDI,					///< surface-assisted laser desorption ionization
	    	SIZE_OF_IONIZATIONMETHOD
	    };
			/// Names of inonization methods
			static const std::string NamesOfIonizationMethod[SIZE_OF_IONIZATIONMETHOD];
      
      /// Polarity of the ion source
      enum Polarity
      {
      	POLNULL,	///< Unknown
      	POSITIVE, ///< Positive polarity
      	NEGATIVE, ///< Negative polarity
      	SIZE_OF_POLARITY
      };
			/// Names of polarity of the ion source
			static const std::string NamesOfPolarity[SIZE_OF_POLARITY];
 			
 			/// Constructor
      IonSource();
      /// Copy constructor
      IonSource(const IonSource& source);
      /// Destructor
      ~IonSource();
 			
 			/// Assignemtn operator
      IonSource& operator= (const IonSource& source);

      /// Equality operator
      bool operator== (const IonSource& rhs) const;
      /// Equality operator
      bool operator!= (const IonSource& rhs) const;
			
			/// returns the inlet type 
      InletType getInletType() const;
      /// sets the  inlet type
      void setInletType(InletType inlet_type);
			
			/// returns the ionization method
      IonizationMethod getIonizationMethod() const;
      /// sets the ionization method
      void setIonizationMethod(IonizationMethod ionization_type);
			
			/// returns the ionization mode
      Polarity getPolarity() const;
      /// sets the ionization mode
      void setPolarity(Polarity polarity);

			/**
				@brief returns the position of this part in the whole Instrument. 
				
				Order can be ignored, as long the instrument has this default setup:
				- one ion source
				- one or many mass analyzers
				- one ion detector
				
				For more complex instuments, the order should be defined.
      */
      Int getOrder() const;
      /// sets the order
      void setOrder(Int order);			

    protected:
	    InletType inlet_type_;
	    IonizationMethod ionization_method_;
			Polarity polarity_;
    	Int order_;

	};

} // namespace OpenMS

#endif // OPENMS_METADATA_IONSOURCE_H
