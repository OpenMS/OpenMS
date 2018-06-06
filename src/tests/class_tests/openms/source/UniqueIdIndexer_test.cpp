// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

struct Dummy
  : UniqueIdInterface
{
    Dummy () : 
      dummy(0) 
    {}

    Size dummy;
};

class DummyVectorIndexed
: public vector<Dummy>,
  public UniqueIdIndexer<DummyVectorIndexed>
{
};

// this is used for testing purposes only
template < typename T >
struct CanAccessTheUniqueIdMap : T
{
  DummyVectorIndexed::UniqueIdMap & getUniqueIdMap() { return this->uniqueid_to_index_; }
};

// this saves us some typing
template < typename T >
CanAccessTheUniqueIdMap<T>& canAccessTheUniqueIdMap(T &rhs)
{
  return static_cast<CanAccessTheUniqueIdMap<T>&>(rhs);
}

}

START_TEST(UniqueIdIndexer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DummyVectorIndexed* ptr = nullptr;
DummyVectorIndexed* nullPointer = nullptr;
START_SECTION((UniqueIdIndexer()))
{
  ptr = new DummyVectorIndexed();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~UniqueIdIndexer()))
{
  delete ptr;
}
END_SECTION

START_SECTION((Size uniqueIdToIndex(UInt64 unique_id) const))
{
  DummyVectorIndexed dvi;
  Size num_uii = 10;
  dvi.resize(num_uii);
  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }

  for ( Size i = 0; i < num_uii; ++i )
  {
    TEST_EQUAL(dvi.uniqueIdToIndex(10*i+1000),i);
  }

  STATUS("shuffling ...");
  std::random_shuffle(dvi.begin(), dvi.end());

  for ( Size i = 0; i < num_uii; ++i )
  {
    Dummy const & current_dummy = dvi[i];
    TEST_EQUAL(dvi.uniqueIdToIndex(current_dummy.getUniqueId()),i);
    // TEST_EQUAL(dvidvi[current_dummy.dummy].getUniqueId()),current_dummy.dummy);
  }

  dvi.pop_back();
  dvi.pop_back();

  dvi.push_back(Dummy());
  dvi.back().setUniqueId(12345678);

  dvi.push_back(Dummy());
  dvi.push_back(Dummy());
  dvi.back().setUniqueId(12345678);
  dvi.push_back(Dummy());

  STATUS("shuffling ...");
  std::random_shuffle(dvi.begin(), dvi.end());

  TEST_EXCEPTION_WITH_MESSAGE(Exception::Postcondition,dvi.updateUniqueIdToIndex(),"Duplicate valid unique ids detected!   RandomAccessContainer has size()==12, num_valid_unique_id==10, uniqueid_to_index_.size()==9");

}
END_SECTION

START_SECTION((void updateUniqueIdToIndex() const))
{
  // see uniqueIdToIndex()
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((Size resolveUniqueIdConflicts()))
  DummyVectorIndexed dvi;
  Size num_uii = 10;
  dvi.resize(num_uii);
  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }
  
  TEST_EQUAL(dvi.resolveUniqueIdConflicts(), 0);
  
  // introduce three doubles
  Dummy a,b;
  a.setUniqueId(1000);
  b.setUniqueId(1000+30);
  dvi.push_back(a);
  dvi.push_back(b);
	dvi.push_back(b);
	TEST_EXCEPTION(Exception::Postcondition,dvi.updateUniqueIdToIndex());
	TEST_EQUAL(dvi.resolveUniqueIdConflicts(), 3);

END_SECTION

START_SECTION((void swap(UniqueIdIndexer &rhs)))
{

  DummyVectorIndexed dvi;
  Size num_uii = 10;

  dvi.resize(num_uii);

  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }

  dvi.updateUniqueIdToIndex();

  DummyVectorIndexed dvi2;

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi2).getUniqueIdMap().size(),0);

  std::swap(dvi, dvi2);

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),0);
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi2).getUniqueIdMap().size(),num_uii);

  dvi = dvi2;

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);
  canAccessTheUniqueIdMap(dvi).getUniqueIdMap().clear();
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),0);

  TEST_EQUAL(dvi.uniqueIdToIndex(4321234324124ull),Size(-1));

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
