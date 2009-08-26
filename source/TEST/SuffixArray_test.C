// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SuffixArray, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((SuffixArray()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((SuffixArray(const String &st, const String &filename)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((SuffixArray(const SuffixArray &sa)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~SuffixArray()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual String toString()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void findSpec(std::vector< std::vector< std::pair< std::pair< SignedSize, SignedSize >, DoubleReal > > > &candidates, const std::vector< DoubleReal > &spec)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool save(const String &filename)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool open(const String &filename)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setTolerance(DoubleReal t)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual DoubleReal getTolerance() const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool isDigestingEnd(const char aa1, const char aa2) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setTags(const std::vector< String > &tags)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual const std::vector<String>& getTags()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setUseTags(bool use_tags)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool getUseTags()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setNumberOfModifications(Size number_of_mods)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual Size getNumberOfModifications()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void printStatistic()=0))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



