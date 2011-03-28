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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PILISCrossValidation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PILISCrossValidation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PILISCrossValidation* ptr = 0;
PILISCrossValidation* nullPointer = 0;
START_SECTION(PILISCrossValidation())
{
	ptr = new PILISCrossValidation();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~PILISCrossValidation())
{
	delete ptr;
}
END_SECTION

START_SECTION((PILISCrossValidation(const PILISCrossValidation &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((PILISCrossValidation& operator=(const PILISCrossValidation &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setOptions(const Map<String, Option> &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setOption(const String &name, const Option &option)))
{
  // TODO
}
END_SECTION

START_SECTION((void apply(Param &PILIS_param, const PILISModel &base_model, const std::vector< Peptide > &peptides)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal scoreHits(const std::vector< std::vector< std::vector< RichPeakSpectrum > > > &sim_spectra, const std::vector< std::vector< RichPeakSpectrum > > &exp_spectra)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



