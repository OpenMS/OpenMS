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
  class InstrumentSettings
    : public MetaInfoInterface
  {
    public:
			/// scan mode
      enum ScanMode
      {
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
      	ER,								///< Enhanced resolution scan
      	SIZE_OF_SCANMODE
      };
			
			/// Names of scan modes
			static const std::string NamesOfScanMode[SIZE_OF_SCANMODE];
			
      /// Scan window type containing m/z begin and end
      struct ScanWindow
      {
        ScanWindow()
          : begin(0.0),
            end(0.0)
        {
        }
        
        bool operator==(const ScanWindow& rhs) const
        {
          return begin == rhs.begin && end == rhs.end;
        }
        
        DoubleReal begin;
        DoubleReal end;
      };        
        
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
			
			/// returns a const reference to the m/z scan windows
      const std::vector< ScanWindow >&  getScanWindows() const;
			/// returns a mutable reference to the m/z scan windows
      std::vector< ScanWindow >&  getScanWindows();
      /// sets the m/z scan windows
      void setScanWindows(std::vector< ScanWindow >  scan_windows);

    protected:
      ScanMode scan_mode_;
      IonSource::Polarity polarity_;
      std::vector< ScanWindow > scan_windows_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_INSTRUMENTSETTINGS_H
