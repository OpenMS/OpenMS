// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

///////////////////////////

START_TEST(Histogram, "$Id$")

/////////////////////////////////////////////////////////////

Histogram<float,float>* dis_ptr = nullptr;
Histogram<float,float>* dis_nullPointer = nullptr;

START_SECTION((Histogram()))
	dis_ptr = new Histogram<float,float>();
  TEST_NOT_EQUAL(dis_ptr, dis_nullPointer)
END_SECTION

START_SECTION((~Histogram()))
	delete dis_ptr;
END_SECTION

Histogram<float,float> d(0, 10, 1);

START_SECTION((Histogram(const Histogram& histogram)))
	Histogram<float, float> d2(d);
	TEST_EQUAL(d == d2, true)
END_SECTION

START_SECTION((BinSizeType minBound() const))
	TEST_REAL_SIMILAR(d.minBound(), 0.0)
END_SECTION

START_SECTION((BinSizeType maxBound() const))
	TEST_REAL_SIMILAR(d.maxBound(), 10.0)
END_SECTION

START_SECTION((BinSizeType binSize() const))
	TEST_REAL_SIMILAR(d.binSize(), 1)
END_SECTION

START_SECTION((Size size() const ))
	TEST_EQUAL(d.size(), 10)
END_SECTION

START_SECTION((Histogram(BinSizeType min, BinSizeType max, BinSizeType bin_size)))
	Histogram<float, float> d3(5.5f, 7.7f, 0.2f);
	TEST_REAL_SIMILAR(d3.minBound(), 5.5)
	TEST_REAL_SIMILAR(d3.maxBound(), 7.7)
	TEST_REAL_SIMILAR(d3.binSize(), 0.2)
END_SECTION

START_SECTION((ValueType minValue() const))
	TEST_REAL_SIMILAR(d.minValue(), 0.0)
END_SECTION

START_SECTION((ValueType maxValue() const))
	TEST_REAL_SIMILAR(d.maxValue(), 0.0)
END_SECTION

START_SECTION((ValueType operator [] (Size index) const))
	d.reset(4, 14, 2);
	TEST_EQUAL(d.size(),5);
	TEST_REAL_SIMILAR(d[0],0.0);
	TEST_REAL_SIMILAR(d[1],0.0);
	TEST_REAL_SIMILAR(d[2],0.0);
	TEST_REAL_SIMILAR(d[3],0.0);
	TEST_REAL_SIMILAR(d[4],0.0);
	TEST_EXCEPTION(Exception::IndexOverflow, d[5])
END_SECTION

START_SECTION((Size inc(BinSizeType val, ValueType increment=1)))
	Size bin_index = 123456;
	TEST_EXCEPTION(Exception::OutOfRange, d.inc(3.9f, 250.3f))
	TEST_EXCEPTION(Exception::OutOfRange, d.inc(14.1f, 250.3f))
		
	bin_index = d.inc(4.0f, 1.0);
	TEST_EQUAL(bin_index,0);
	bin_index = d.inc(5.9f, 1.0);
	TEST_EQUAL(bin_index,0);
	
	TEST_REAL_SIMILAR(d[0],2.0);
	TEST_REAL_SIMILAR(d[1],0.0);
	TEST_REAL_SIMILAR(d[2],0.0);
	TEST_REAL_SIMILAR(d[3],0.0);
	TEST_REAL_SIMILAR(d[4],0.0);
	
	bin_index = d.inc(8.0f, 45.0);
	TEST_EQUAL(bin_index,2);
	bin_index = d.inc(8.1f, 1.0);
	TEST_EQUAL(bin_index,2);
	bin_index = d.inc(9.9f, 4.0);
	TEST_EQUAL(bin_index,2);

	TEST_REAL_SIMILAR(d[0],2.0);
	TEST_REAL_SIMILAR(d[1],0.0);
	TEST_REAL_SIMILAR(d[2],50.0);
	TEST_REAL_SIMILAR(d[3],0.0);
	TEST_REAL_SIMILAR(d[4],0.0);

	bin_index = d.inc(12.0f, 1.0);
	TEST_EQUAL(bin_index,4);
	bin_index = d.inc(13.1f, 2.0);
	TEST_EQUAL(bin_index,4);
	bin_index = d.inc(14.0f, 3.0);	
	TEST_EQUAL(bin_index,4);

	TEST_REAL_SIMILAR(d[0],2.0);
	TEST_REAL_SIMILAR(d[1],0.0);
	TEST_REAL_SIMILAR(d[2],50.0);
	TEST_REAL_SIMILAR(d[3],0.0);
	TEST_REAL_SIMILAR(d[4],6.0);
END_SECTION

START_SECTION((ConstIterator begin() const))
	Histogram<float,float>::ConstIterator it = d.begin();
	TEST_REAL_SIMILAR(*it, 2.0)
END_SECTION

START_SECTION((ConstIterator end() const))
	Histogram<float,float>::ConstIterator it = d.begin();
	TEST_REAL_SIMILAR(*it,2.0);
	++it;
	TEST_REAL_SIMILAR(*it,0.0);
	++it;
	TEST_REAL_SIMILAR(*it,50.0);
	++it;
	TEST_REAL_SIMILAR(*it,0.0);
	++it;
	TEST_REAL_SIMILAR(*it,6.0);
	++it;
	TEST_EQUAL(it==d.end(),true);
END_SECTION

START_SECTION((ValueType binValue(BinSizeType val) const))
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(3.9f))
	TEST_REAL_SIMILAR(d.binValue(4.0f),2.0);
	TEST_REAL_SIMILAR(d.binValue(5.9f),2.0);
	TEST_REAL_SIMILAR(d.binValue(6.0f),0.0);
	TEST_REAL_SIMILAR(d.binValue(7.9f),0.0);
	TEST_REAL_SIMILAR(d.binValue(8.0f),50.0);
	TEST_REAL_SIMILAR(d.binValue(9.9f),50.0);
	TEST_REAL_SIMILAR(d.binValue(10.0f),0.0);
	TEST_REAL_SIMILAR(d.binValue(11.9f),0.0);
	TEST_REAL_SIMILAR(d.binValue(12.0f),6.0);
	TEST_REAL_SIMILAR(d.binValue(14.0f),6.0);
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(14.1f))
END_SECTION
	
START_SECTION((void reset(BinSizeType min, BinSizeType max, BinSizeType bin_size)))
	d.reset(1, 11, 2);
	TEST_REAL_SIMILAR(d.minBound(), 1)
	TEST_REAL_SIMILAR(d.maxBound(), 11)
	TEST_EQUAL(d.size(), 5)
	TEST_REAL_SIMILAR(d.binSize(), 2)
END_SECTION

START_SECTION((bool operator == (const Histogram& histogram) const))
	Histogram<float, float> dist(1, 11, 2);
	TEST_EQUAL(d == dist, true)
END_SECTION

START_SECTION((bool operator != (const Histogram& histogram) const))
	Histogram<float, float> dist(1, 12, 2);
	TEST_EQUAL(d != dist, true)
END_SECTION

START_SECTION((Histogram& operator = (const Histogram& histogram)))
	Histogram<float, float> dist;
	dist = d;
	TEST_EQUAL(d == dist, true)
END_SECTION

START_SECTION((void applyLogTransformation(BinSizeType multiplier)))
	TOLERANCE_ABSOLUTE(0.01)
	Histogram<float, float> dist(0,5,1);
	dist.inc(0.5,1);
	dist.inc(1.5,10);
	dist.inc(2.5,100);
	dist.inc(3.5,1000);
	dist.inc(4.5,10000);
	dist.applyLogTransformation(1.0);
	TEST_REAL_SIMILAR(dist.binValue(0.5),0.6931);
	TEST_REAL_SIMILAR(dist.binValue(1.5),2.3979);
	TEST_REAL_SIMILAR(dist.binValue(2.5),4.61512);
	TEST_REAL_SIMILAR(dist.binValue(3.5),6.90875);
	TEST_REAL_SIMILAR(dist.binValue(4.5),9.21044);
END_SECTION

START_SECTION((BinSizeType centerOfBin(Size bin_index) const))
	Histogram<float, float> dist(0,5,1);
	dist.inc(0.5,1);
	dist.inc(1.5,10);
	dist.inc(2.5,100);
	dist.inc(3.5,1000);
	dist.inc(4.5,10000);
	TEST_REAL_SIMILAR(dist.centerOfBin(0),0.5);
	TEST_REAL_SIMILAR(dist.centerOfBin(1),1.5);
	TEST_REAL_SIMILAR(dist.centerOfBin(2),2.5);
	TEST_REAL_SIMILAR(dist.centerOfBin(3),3.5);
	TEST_REAL_SIMILAR(dist.centerOfBin(4),4.5);
	TEST_EXCEPTION(Exception::IndexOverflow, dist.centerOfBin(5))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
