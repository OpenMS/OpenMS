// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_METADATA_IONSOURCE_H
#define OPENMS_METADATA_IONSOURCE_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of a ion source ( Part of a MS Instrument )
		
		@ingroup Metadata
	*/
  class IonSource: public MetaInfoInterface
  {
    public:
    	/// inlet type
	    enum InletType{INLETNULL,DIRECT, BATCH, CHROMATOGRAPHY, PARTICLEBEAM, MEMBRANESEPARATOR, OPENSPLIT, JETSEPARATOR,SEPTUM, RESERVOIR, MOVINGBELT, MOVINGWIRE, FLOWINJECTIONANALYSIS, ELECTROSPRAYINLET,THERMOSPRAYINLET, INFUSION, CONTINUOUSFLOWFASTATOMBOMBARDMENT, INDUCTIVELYCOUPLEDPLASMA, SIZE_OF_INLETTYPE};
			/// Names of inlet types
			static const std::string NamesOfInletType[SIZE_OF_INLETTYPE];

	    /// ionization method
	    enum IonizationMethod{IONMETHODNULL,ESI,EI,CI,FAB,TSP,LD,FD,FI,PD,SI,TI,API,ISI,CID,CAD,HN,APCI,APPI,ICP,SIZE_OF_IONIZATIONMETHOD};
			/// Names of inonization methods
			static const std::string NamesOfIonizationMethod[SIZE_OF_IONIZATIONMETHOD];
      
      /// polarity of the ion source
      enum Polarity{POLNULL,POSITIVE, NEGATIVE,SIZE_OF_POLARITY};
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

    protected:
	    InletType inlet_type_;
	    IonizationMethod ionization_method_;
			Polarity polarity_;

	};

} // namespace OpenMS

#endif // OPENMS_METADATA_IONSOURCE_H
