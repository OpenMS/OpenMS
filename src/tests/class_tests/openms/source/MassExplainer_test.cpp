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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Compomer.h>

using namespace OpenMS;
using namespace std;

START_TEST(MassExplainer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassExplainer* ptr = nullptr;
MassExplainer* nullPointer = nullptr;
START_SECTION(MassExplainer())
{
  ptr = new MassExplainer();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MassExplainer())
{
	delete ptr;
}
END_SECTION

START_SECTION((MassExplainer(AdductsType adduct_base)))
{
	Adduct a(2, 1, 123.12, "Na", -0.5,0);
	MassExplainer::AdductsType va;
	va.push_back(a);  
	MassExplainer me(va);
	TEST_EQUAL(me.getAdductBase().size(), 1);
}
END_SECTION

START_SECTION((MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp)))
{
	MassExplainer me(5,10,2,-10.3);
	TEST_EQUAL(me.getAdductBase().size(), 4);
}
END_SECTION

START_SECTION((MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals)))
{
	MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me(va,5,10,2,-10.3,0);
	TEST_EQUAL(me.getAdductBase().size(), 2);
}
END_SECTION

START_SECTION((MassExplainer& operator=(const MassExplainer &rhs)))
{
  MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me(va,5,10,2,-10.3,0);
	MassExplainer me2;
	me2 = me;
	TEST_EQUAL(me.getAdductBase().size(), 2);

}
END_SECTION

START_SECTION((void setAdductBase(AdductsType adduct_base)))
{
	MassExplainer::AdductsType va;
	Adduct a1(2, 1, 123.12, "Na", -0.5,0);
	Adduct a2(3, 1, 123.12, "K", -0.7,0);
	va.push_back(a1);  
	va.push_back(a2);  
	MassExplainer me;
	me.setAdductBase(va);
	TEST_EQUAL(me.getAdductBase().size(),2);
}
END_SECTION

START_SECTION((AdductsType getAdductBase() const))
{
  NOT_TESTABLE; // tested above
}
END_SECTION

START_SECTION((const Compomer& getCompomerById(Size id) const))
{
  MassExplainer me;
	me.compute();
	TEST_EQUAL(me.getCompomerById(0).getID(), 0)
}
END_SECTION

START_SECTION((void compute()))
{
  NOT_TESTABLE; // tested above
}
END_SECTION

START_SECTION((SignedSize query(const Int net_charge, const float mass_to_explain, const float mass_delta, const float thresh_log_p, std::vector< Compomer >::const_iterator &firstExplanation, std::vector< Compomer >::const_iterator &lastExplanation) const))
{
  MassExplainer me;
	me.compute();
	
	
	MassExplainer::CompomerIterator s,e;
	SignedSize hits = me.query(2, 45.0, 13.0, -100000,s,e);
	
	std::cout << "hits: " << hits << std::endl;
	for (; s!=e; ++s)
	{
			std::cout << *s << std::endl;
	}
	
	TEST_EQUAL(hits, 5);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
