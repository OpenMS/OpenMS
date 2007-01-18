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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/Benchmark.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

using namespace OpenMS;

START_BENCHMARK(Core, 1.0, "$Id: Core_bench.C,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(DPeak<4> construction, 1.0)

	for (int count = 0; count < 200000; count++)
	{
		DPeak<4>* e_ptr;
		START_TIMER
			e_ptr = new DPeak<4>;
		STOP_TIMER
		delete e_ptr;
	}

END_SECTION

START_SECTION(DPeak<4> deconstruction, 1.0)

	for (int count = 0; count < 200000; count++)
	{
		DPeak<4>* e_ptr = new DPeak<4>;
		START_TIMER
			delete e_ptr;
		STOP_TIMER
	}

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_BENCHMARK
