// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IONDETECTOR_H
#define OPENMS_METADATA_IONDETECTOR_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of a ion detector (part of a MS Instrument)
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI IonDetector
  	: public MetaInfoInterface
  {
    public:
    	/// Detector type
			enum Type
			{
				TYPENULL,														///< Unknown
				ELECTRONMULTIPLIER,									///< Electron multiplier
				PHOTOMULTIPLIER,										///< Photo multiplier
				FOCALPLANEARRAY,										///< Focal plane array
				FARADAYCUP,													///< Faraday cup
				CONVERSIONDYNODEELECTRONMULTIPLIER,	///< Conversion dynode electron multiplier
				CONVERSIONDYNODEPHOTOMULTIPLIER,		///< Conversion dynode photo multiplier
				MULTICOLLECTOR,											///< Multi-collector
				CHANNELELECTRONMULTIPLIER,					///< Channel electron multiplier
				CHANNELTRON,												///< channeltron
				DALYDETECTOR,												///< daly detector
				MICROCHANNELPLATEDETECTOR,					///< microchannel plate detector
				ARRAYDETECTOR,											///< array detector
				CONVERSIONDYNODE,										///< conversion dynode
				DYNODE,															///< dynode
				FOCALPLANECOLLECTOR,								///< focal plane collector
				IONTOPHOTONDETECTOR,								///< ion-to-photon detector
				POINTCOLLECTOR,											///< point collector
				POSTACCELERATIONDETECTOR,						///< postacceleration detector
				PHOTODIODEARRAYDETECTOR,						///< photodiode array detector
				INDUCTIVEDETECTOR,									///< inductive detector
				ELECTRONMULTIPLIERTUBE,							///< electron multiplier tube
				SIZE_OF_TYPE
			};
			/// Names of detector types
			static const std::string NamesOfType[SIZE_OF_TYPE];

			/// Acquisition mode
			enum AcquisitionMode
			{
				ACQMODENULL,				///< Unknown
				PULSECOUNTING,			///< Pulse counting
				ADC,								///< Analog-digital converter
				TDC,								///< Time-digital converter
				TRANSIENTRECORDER,	///< Transient recorder
				SIZE_OF_ACQUISITIONMODE
			};
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
      DoubleReal getResolution() const;
      /// sets the resolution (in ns)
      void setResolution(DoubleReal resolution);
			
			/// retruns the analog-to-digital converter sampling frequency (in Hz)
      DoubleReal getADCSamplingFrequency() const;
      /// sets the analog-to-digital converter sampling frequency (in Hz)
      void setADCSamplingFrequency(DoubleReal ADC_sampling_frequency);

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
	    Type type_;
	    AcquisitionMode acquisition_mode_; 
	    DoubleReal resolution_;
	    DoubleReal ADC_sampling_frequency_;
    	Int order_;
    
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_IONDETECTOR_H
