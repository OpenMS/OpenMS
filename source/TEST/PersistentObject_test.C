// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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



