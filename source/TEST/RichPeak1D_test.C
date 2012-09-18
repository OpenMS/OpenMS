// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RichPeak1D.h>

///////////////////////////

START_TEST(RichPeak1D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RichPeak1D* d10_ptr = 0;
RichPeak1D* d10_nullPointer = 0;
START_SECTION((RichPeak1D()))
	d10_ptr = new RichPeak1D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~RichPeak1D()))
	delete d10_ptr;
END_SECTION

START_SECTION((RichPeak1D(const RichPeak1D &p)))
	RichPeak1D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak1D copy_of_p(p);

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((RichPeak1D& operator=(const RichPeak1D &rhs)))
	RichPeak1D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak1D copy_of_p;
	copy_of_p = p;

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((bool operator == (const RichPeak1D& rhs) const))
	RichPeak1D p1, p2;
	TEST_EQUAL(p1==p2, true)
	
	p1.setIntensity(5.0f);
	TEST_EQUAL(p1==p2, false)
	p2.setIntensity(5.0f);
	TEST_EQUAL(p1==p2, true)

	p1.setMetaValue("cluster_id",4711);
	TEST_EQUAL(p1==p2, false)
	p1.removeMetaValue("cluster_id");
	TEST_EQUAL(p1==p2, true)		
END_SECTION

START_SECTION((bool operator != (const RichPeak1D& rhs) const))
	RichPeak1D p1, p2;
	TEST_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5.0f);
	TEST_EQUAL(p1!=p2, true)
	p2.setIntensity(5.0f);
	TEST_EQUAL(p1!=p2, false)

	p1.setMetaValue("cluster_id",4711);
	TEST_EQUAL(p1!=p2, true)
	p1.removeMetaValue("cluster_id");
	TEST_EQUAL(p1!=p2, false)	
END_SECTION

START_SECTION(([EXTRA] meta info with copy constructor))
	RichPeak1D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak1D p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA] meta info with assignment))
	RichPeak1D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak1D p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
