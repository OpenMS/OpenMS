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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>

#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(<Factory>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Factory is singleton, therefore we don't test the constructor
START_SECTION(static FactoryProduct* create(const String& name))
	FilterFunctor* p = Factory<FilterFunctor>::create("TICFilter");
	TICFilter reducer;
	TEST_EQUAL(*p==reducer,true);
END_SECTION

START_SECTION( static void registerProduct(const String& name, const FunctionType creator) )
	Factory<FilterFunctor>::registerProduct(TICFilter::getProductName(), &TICFilter::create);
	FilterFunctor* ext = Factory<FilterFunctor>::create("TICFilter");
  FilterFunctor* nullPointer = nullptr;
  TEST_NOT_EQUAL(ext, nullPointer)
END_SECTION

START_SECTION(static bool isRegistered(const String& name))
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter"), true)
	TEST_EQUAL(Factory<FilterFunctor>::isRegistered("TICFilter_bla_bluff"), false)
END_SECTION

START_SECTION(static std::vector<String> registeredProducts())
	vector<String> list = Factory<FilterFunctor>::registeredProducts();
	TEST_EQUAL(list.size(),6)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


