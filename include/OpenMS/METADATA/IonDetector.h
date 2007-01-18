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

#ifndef OPENMS_METADATA_IONDETECTOR_H
#define OPENMS_METADATA_IONDETECTOR_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of a ion detector ( part of a MS Instrument )
		
		
		
		@ingroup Metadata
	*/
  class IonDetector: public MetaInfoInterface
  {
    public:
    	/// detector type
			enum Type {TYPENULL,ELECTRONMULTIPLIER,PHOTOMULTIPLIER,FOCALPLANEARRAY,FARADAYCUP,CONVERSIONDYNODEELECTRONMULTIPLIER,CONVERSIONDYNODEPHOTOMULTIPLIER,MULTICOLLECTOR,CHANNELELECTRONMULTIPLIER,SIZE_OF_TYPE};
			/// Names of detector types
			static const std::string NamesOfType[SIZE_OF_TYPE];
			
			/// acquisition mode
			enum AcquisitionMode {ACQMODENULL,PULSECOUNTING,ADC,TDC,TRANSIENTRECORDER,SIZE_OF_ACQUISITIONMODE};
			/// Names of acquisition modes
			static const std::string NamesOfAcquisitionMode[SIZE_OF_ACQUISITIONMODE];

			/// Constructor
      IonDetector();
      /// Copy constructor
      IonDetector(const IonDetector& source);
      /// Destructor
      ~IonDetector();
			
			/// Assignment operator
      IonDetector& operator= (const IonDetector& source);

      /// Equality operator
      bool operator== (const IonDetector& rhs) const;
      /// Equality operator
      bool operator!= (const IonDetector& rhs) const;
			
			/// returns the detector type
      Type getType() const;
      /// sets the detector type
      void setType(Type type);
			
			/// returns the acquisition mode
      AcquisitionMode getAcquisitionMode() const;
      /// sets the acquisition mode
      void setAcquisitionMode(AcquisitionMode acquisition_mode);
			
			/// returns the resolution (in ns)
      float getResolution() const;
      /// sets the resolution (in ns)
      void setResolution(float resolution);
			
			/// retruns the analog-to-digital converter sampling frequency (in MHz)
      float getADCSamplingFrequency() const;
      /// sets the analog-to-digital converter sampling frequency (in MHz)
      void setADCSamplingFrequency(float ADC_sampling_frequency);

    protected:
	    Type type_;
	    AcquisitionMode acquisition_mode_; 
	    float resolution_;
	    float ADC_sampling_frequency_;
    
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_IONDETECTOR_H
