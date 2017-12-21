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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/UnaryComposeFunctionAdapter.h>
///////////////////////////

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

using namespace OpenMS;
using namespace std;

// test class
class Element {
 public:
  Element(const string& a) : _a(a) {}
  const string& getA() const { return _a; }
 private:
  string _a;
};

START_TEST(UnaryComposeFunctionAdapter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef equal_to<string> T_EQUAL_TO;
typedef binder2nd<T_EQUAL_TO> T_BINDER_2ND;
typedef std::const_mem_fun_ref_t<const std::string&, Element> T_MEMBER_FUNCTIOM;

typedef UnaryComposeFunctionAdapter<T_BINDER_2ND, T_MEMBER_FUNCTIOM> UCFA;

UCFA* ptr = nullptr;
UCFA* null_ptr = nullptr;
START_SECTION((UnaryComposeFunctionAdapter(const OP1 &o1, const OP2 &o2)))
{
  ptr = new UCFA(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA));
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~UnaryComposeFunctionAdapter())
{
	delete ptr;
}
END_SECTION

START_SECTION((OP1::result_type operator()(const typename OP2::argument_type &x) const ))
{
  Element a("4"), b("3"), c("2"), d("1");

  vector<Element> elements;
  elements.push_back(a);
  elements.push_back(b);
  elements.push_back(c);
  elements.push_back(d);

  vector<Element>::const_iterator found = find_if(elements.begin(), elements.end(),
    unaryCompose(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA)));

  TEST_EQUAL(found != elements.end(), true)
  TEST_EQUAL(found - elements.begin(), 1)

  vector<Element>::const_iterator not_found = find_if(elements.begin(), elements.end(),
    unaryCompose(bind2nd(equal_to<string>(), "10"), mem_fun_ref(&Element::getA)));

  TEST_EQUAL(not_found == elements.end(), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



