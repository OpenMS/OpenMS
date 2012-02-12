// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/DB/PersistentObject.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	// Test class
	class Dummy: public PersistentObject
	{
		public:
			
			Dummy()
				:PersistentObject(),
				subobjects_clear_(false)
				{
					
				}
		
			Dummy& operator=(const Dummy& rhs)
			{
				if (&rhs==this) return *this;
				
				PersistentObject::operator=(rhs);
				subobjects_clear_ = rhs.subobjects_clear_;
				
				return *this;
			}
			
			bool subobjectsClear() const
			{
				return subobjects_clear_;
			}
		
		protected:
			//emulation of subobjects
			bool subobjects_clear_;

			void clearChildIds_()
			{
				subobjects_clear_ = true;
			}

	};
}

START_TEST(PersistentObject, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Dummy* ptr = 0;
Dummy* nullPointer = 0;
START_SECTION((PersistentObject()))
	ptr = new Dummy();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~PersistentObject()))
	delete ptr;
END_SECTION

START_SECTION((const UID& getPersistenceId() const))
  Dummy tmp;
  TEST_EQUAL(tmp.getPersistenceId(),0)
END_SECTION

START_SECTION((void setPersistenceId(const UID& persistence_id)))
  Dummy tmp;
  tmp.setPersistenceId(4711);
  TEST_EQUAL(tmp.getPersistenceId(),4711)
END_SECTION

START_SECTION((void clearId(bool deep = true)))
  Dummy tmp;
  tmp.setPersistenceId(4711);
  tmp.clearId(false);
  TEST_EQUAL(tmp.getPersistenceId(),0)
  TEST_EQUAL(tmp.subobjectsClear(),false)
  
  tmp.setPersistenceId(4712);
  tmp.clearId(true);
  TEST_EQUAL(tmp.getPersistenceId(),0)
  TEST_EQUAL(tmp.subobjectsClear(),true)
END_SECTION

START_SECTION((PersistentObject& operator= (const PersistentObject& rhs)))
  Dummy tmp;
  tmp.clearId(true);
  tmp.setPersistenceId(4711);
  
  Dummy tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getPersistenceId(),4711)
  TEST_EQUAL(tmp2.subobjectsClear(),true)
  
  tmp2 = Dummy();
  TEST_EQUAL(tmp2.getPersistenceId(),0)
  TEST_EQUAL(tmp2.subobjectsClear(),false)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



