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

#ifndef OPENMS_METADATA_INSTRUMENT_H
#define OPENMS_METADATA_INSTRUMENT_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/IonDetector.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Description of a MS instrument
		
		It contains information like vendor, model, ion source, mass analyzer(s), and ion detector.
		
		@ingroup Metadata
	*/
  class Instrument: public MetaInfoInterface
  {
    public:
    	/// Constructor
      Instrument();
      /// Copy constructor
      Instrument(const Instrument& source);
      /// Destructor
      ~Instrument();
      
      /// Assignement operator
      Instrument& operator= (const Instrument& source);

      /// Equality operator
      bool operator== (const Instrument& rhs) const;
      /// Equality operator
      bool operator!= (const Instrument& rhs) const;
			
			/// returns the name of the instrument
      const String& getName() const;
      /// sets the name of the instrument
      void setName(const String& name);
			
			/// returns the instrument vendor
      const String& getVendor() const;
      /// sets the instrument vendor
      void setVendor(const String& vendor);
			
			/// returns the instrument model
      const String& getModel() const;
      /// sets the instrument model
      void setModel(const String& model);
			
			/// returns a description of constumizations
      const String& getCustomizations() const;
      /// sets the a description of constumizations 
      void setCustomizations(const String& customizations);
			
			/// returns a const reference to the ion source
      const IonSource& getIonSource() const;
      /// returns a mutable reference to the ion source
      IonSource& getIonSource();
      /// sets the ion source
      void setIonSource(const IonSource& ion_source);
			
			/// returns a const reference to the mass analyer list
      const std::vector<MassAnalyzer>& getMassAnalyzers() const;
      /// returns a mutable reference to the mass analyzer list 
      std::vector<MassAnalyzer>& getMassAnalyzers();
      /// sets the mass analyzer list
      void setMassAnalyzers(const std::vector<MassAnalyzer>& mass_analyzers);
			
			/// returns a const reference to the ion detector
      const IonDetector& getIonDetector() const;
      /// returns a mutable reference to the ion detector
      IonDetector& getIonDetector();
      /// sets the ion detector
      void setIonDetector(const IonDetector& ion_detector);

    protected:
		  String name_;
		  String vendor_;
		  String model_;
		  String customizations_;
		  IonSource ion_source_;
		  std::vector<MassAnalyzer> mass_analyzers_;
		  IonDetector ion_detector_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_INSTRUMENT_H
