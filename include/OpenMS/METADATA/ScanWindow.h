// -*- mode: C++; tab-width: 2; -*-s
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SCANWINDOW_H
#define OPENMS_METADATA_SCANWINDOW_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Scan window description
		
		@ingroup Metadata
	*/
	struct OPENMS_DLLAPI ScanWindow
		: public MetaInfoInterface
	{
		///Default constructor
	  ScanWindow();
	  ///Copy constructor
	  ScanWindow(const ScanWindow& source);
	  ///Equality operator
	  bool operator==(const ScanWindow& source) const;
	  ///Equality operator
	  bool operator!=(const ScanWindow& source) const;
		///Assignment operator		
	  ScanWindow& operator=(const ScanWindow& source);
	  
	  ///Begin of the window
	  DoubleReal begin;
		///End of the window
	  DoubleReal end;
	};        
        
} // namespace OpenMS

#endif // OPENMS_METADATA_SCANWINDOW_H
