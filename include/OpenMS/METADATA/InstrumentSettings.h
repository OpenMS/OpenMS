// -*- mode: C++; tab-width: 2; -*-s
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

#ifndef OPENMS_METADATA_INSTRUMENTSETTINGS_H
#define OPENMS_METADATA_INSTRUMENTSETTINGS_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IonSource.h>

namespace OpenMS 
{
	/**
		@brief Description of the settings a MS Instrument was run with.
		
		
		
		@ingroup Metadata
	*/
  class InstrumentSettings: public MetaInfoInterface
  {
    public:
			/// scan mode
      enum ScanMode{
      	UNKNOWN,					///< Unknown scan method
      	ZOOM,							///< Zoom scan
      	FULL,							///< Full scan
      	SIM,							///< Selected ion monitoring
      	SRM,							///< Selected reaction monitoring
      	CRM,							///< Consecutive reaction monitoring
      	CNG,							///< Constant neutral gain scan
      	CNL,							///< Constant neutral loss scan
      	PRODUCT,					///< Product ion scan
      	PRECURSOR,				///< Precursor ion scan
      	SIZE_OF_SCANMODE};
			
			/// Names of scan modes
			static const std::string NamesOfScanMode[SIZE_OF_SCANMODE];
			
			/// Constructor
      InstrumentSettings();
      /// Copy constructor
      InstrumentSettings(const InstrumentSettings& source);
      /// Destructor
      ~InstrumentSettings();
      
      /// Assignment operator
      InstrumentSettings& operator= (const InstrumentSettings& source);

      /// Equality operator
      bool operator== (const InstrumentSettings& rhs) const;
      /// Equality operator
      bool operator!= (const InstrumentSettings& rhs) const;

			/// returns the scan mode
      ScanMode getScanMode() const;
      /// sets the scan mode
      void setScanMode(ScanMode scan_mode);
			
			/// returns the polarity
      IonSource::Polarity getPolarity() const;
      /// sets the polariy
      void setPolarity(IonSource::Polarity polarity);
			
			/// returns the scan begin in m/z dimension (default is 0.0)
      float getMzRangeStart() const;
      /// sets the scan begin in m/z dimension
      void setMzRangeStart(float mz_range_start);
			
			/// returns the scan end in m/z dimension (default is 0.0)
      float getMzRangeStop() const;
      /// sets the scan end in m/z dimension
      void setMzRangeStop(float mz_range_stop);

    protected:
      ScanMode scan_mode_;
      IonSource::Polarity polarity_;
      float mz_range_start_;
      float mz_range_stop_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_INSTRUMENTSETTINGS_H
