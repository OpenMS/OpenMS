// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow$
// $Authors: Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_SEQANINCLUDEWRAPPER_H
#define OPENMS_DATASTRUCTURES_SEQANINCLUDEWRAPPER_H

#ifdef _MSC_VER // disable some seqan warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4244 4267 4390 4521 4522 4800)
#endif

#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/graph_align.h> 

#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif

#endif // OPENMS_DATASTRUCTURES_SEQANINCLUDEWRAPPER_H

