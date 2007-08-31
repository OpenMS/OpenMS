// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bauer$
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticSeqan.h>

using namespace OpenMS;

SuffixArrayTrypticSeqan::SuffixArrayTrypticSeqan(const String & st,const String & sa_file_name) throw (Exception::InvalidValue,Exception::FileNotFound):SuffixArraySeqan(st, sa_file_name)
{
	//super(st,sa_file_name);
}

bool SuffixArrayTrypticSeqan::isDigestingEnd(const char aa1, const char aa2) const
{
	return (aa1 == 'K' || aa1 == 'R') && aa2 != 'P';
}

