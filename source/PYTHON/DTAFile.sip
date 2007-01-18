// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Ulf Leser, Knut Reinert
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
// $Id: DTAFile.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//


class DTAFile
{
%TypeHeaderCode	
	#include <OpenMS/FORMAT/DTAFile.h>
	using namespace OpenMS;
%End

	public:

	typedef DDiscreteSpectrum1 SpectrumType;
	DTAFile();
	DTAFile(const std::string& /filename/);
	DTAFile(const DTAFile& /file/);
	~DTAFile();

	const std::string& filename() const;
	std::string& filename();

	DTAFile& operator >> (SpectrumType& /spectrum/);
	DTAFile& operator << (const SpectrumType& /spectrum/);
};
