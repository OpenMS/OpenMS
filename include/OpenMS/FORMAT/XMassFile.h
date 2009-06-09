// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMASSFILE_H
#define OPENMS_FORMAT_XMASSFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
// #include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>
#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>

namespace OpenMS
{
 	/**
 		@brief File adapter for 'XMass Analysis (fid)' files.
 		
 		MZ calculs are based on article :
      A database application for pre-processing, storage and comparison of mass spectra derived from patients and controls
      Mark K Titulaer, Ivar Siccama, Lennard J Dekker, Angelique LCT van Rijswijk, Ron MA Heeren, Peter A Sillevis Smitt, and Theo M Luider
      BMC Bioinformatics. 2006; 7: 403
      http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1594579&blobtype=pdf
  	
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI XMassFile
  {
// 	  private:
//		  PeakFileOptions options_;

    public:
      /// Default constructor
      XMassFile();
      
      template <typename SpectrumType>
      void load(const String& filename, SpectrumType& spectrum)
      {
        Internal::AcqusHandler acqus(filename.prefix(filename.length()-3) + String("acqus"));

        Internal::FidHandler fid(filename);
        if (!fid)
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}      
				
        //  Delete old spectrum
				spectrum.clear();
				
				//temporary variables
				typename SpectrumType::PeakType p;
				
// http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS-release/html/classOpenMS_1_1SpectrumSettings.html#5b198c3f5fab186ee03efe07e39a11dd
				
        while( spectrum.size() < acqus.getSize() )
        {
				  //fill peak
				  p.setPosition( (typename SpectrumType::PeakType::PositionType) acqus.getPosition(fid.getIndex()) );
				  p.setIntensity( (typename SpectrumType::PeakType::IntensityType) fid.getIntensity() );
				  spectrum.push_back(p);
        }
        fid.close();
      }

      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("'Xmass' file not writable."));
      }
      
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_XMASSFILE_H

