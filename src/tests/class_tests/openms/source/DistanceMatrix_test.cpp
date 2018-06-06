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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DistanceMatrix, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DistanceMatrix<double>* ptr = nullptr;
DistanceMatrix<double>* nullPointer = nullptr;
START_SECTION(DistanceMatrix())
{
	ptr = new DistanceMatrix<double>();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DistanceMatrix())
{
	delete ptr;
}
END_SECTION

DistanceMatrix<double> dm(8,1.0);

START_SECTION((DistanceMatrix(SizeType dimensionsize, Value value=Value())))
{
	TEST_EQUAL(dm.dimensionsize(), 8)
	TEST_EQUAL(dm(6,7),1)
}
END_SECTION

DistanceMatrix<double> dm2(dm);

START_SECTION((DistanceMatrix(const DistanceMatrix &source)))
{
	TEST_EQUAL(dm2.dimensionsize(), 8)
	TEST_EQUAL(dm2(2,3),1)
}
END_SECTION

START_SECTION(void resize(SizeType dimensionsize, Value value=Value()))
{
	dm2.resize(5);
	TEST_EQUAL(dm2.dimensionsize(),5)
}
END_SECTION

START_SECTION((SizeType dimensionsize() const))
{
	TEST_EQUAL(dm2.dimensionsize(),5)
}
END_SECTION

START_SECTION((void setValue(SizeType i, SizeType j, ValueType value)))
{
	dm.setValue(0,1,10);
	dm.setValue(0,2,9);
	dm.setValue(0,3,8);
	dm.setValue(0,4,7);
	dm.setValue(1,2,6);
	dm.setValue(1,3,5);
	dm.setValue(1,4,4);
	dm.setValue(2,3,3);
	dm.setValue(2,4,2);
	dm.setValue(3,4,0.5);
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),0.5)
	dm.setValue(3,4,1);
	TEST_EQUAL(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),1.0)
	//more tested below
}
END_SECTION

START_SECTION((const ValueType getValue(SizeType i, SizeType j) const))
{
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(0,2),9)
	TEST_EQUAL(dm.getValue(0,3),8)
	TEST_EQUAL(dm.getValue(0,4),7)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),5)
	TEST_EQUAL(dm.getValue(1,4),4)
	TEST_EQUAL(dm.getValue(2,3),3)
	TEST_EQUAL(dm.getValue(2,4),2)
	TEST_EQUAL(dm.getValue(3,4),1)
}
END_SECTION

START_SECTION((ValueType getValue(SizeType i, SizeType j)))
{
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(0,2),9)
	TEST_EQUAL(dm.getValue(0,3),8)
	TEST_EQUAL(dm.getValue(0,4),7)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),5)
	TEST_EQUAL(dm.getValue(1,4),4)
	TEST_EQUAL(dm.getValue(2,3),3)
	TEST_EQUAL(dm.getValue(2,4),2)
	TEST_EQUAL(dm.getValue(3,4),1)
}
END_SECTION

START_SECTION((void clear()))
{
	dm2.clear();
	TEST_EQUAL(dm2.dimensionsize(),0)
}
END_SECTION

START_SECTION((void setValueQuick(SizeType i, SizeType j, ValueType value)))
{
	dm.setValueQuick(0,1,1);
	dm.setValueQuick(0,2,2);
	dm.setValueQuick(0,3,3);
	dm.setValueQuick(0,4,4);
	dm.setValueQuick(1,2,5);
	dm.setValueQuick(1,3,6);
	dm.setValueQuick(1,4,7);
	dm.setValueQuick(2,3,8);
	dm.setValueQuick(2,4,9);
	dm.setValueQuick(3,4,10);
	TEST_EQUAL(dm.getValue(0,1),1)
	TEST_EQUAL(dm.getValue(0,2),2)
	TEST_EQUAL(dm.getValue(0,3),3)
	TEST_EQUAL(dm.getValue(0,4),4)
	TEST_EQUAL(dm.getValue(1,2),5)
	TEST_EQUAL(dm.getValue(1,3),6)
	TEST_EQUAL(dm.getValue(1,4),7)
	TEST_EQUAL(dm.getValue(2,3),8)
	TEST_EQUAL(dm.getValue(2,4),9)
	TEST_EQUAL(dm.getValue(3,4),10)
}
END_SECTION

START_SECTION((const ValueType operator()(SizeType i, SizeType j) const))
{
	TEST_EQUAL(dm.getValue(0,1),dm(0,1))
	TEST_EQUAL(dm.getValue(0,2),dm(0,2))
	TEST_EQUAL(dm.getValue(0,3),dm(0,3))
	TEST_EQUAL(dm.getValue(0,4),dm(0,4))
	TEST_EQUAL(dm.getValue(1,2),dm(1,2))
	TEST_EQUAL(dm.getValue(1,3),dm(1,3))
	TEST_EQUAL(dm.getValue(1,4),dm(1,4))
	TEST_EQUAL(dm.getValue(2,3),dm(2,3))
	TEST_EQUAL(dm.getValue(2,4),dm(2,4))
	TEST_EQUAL(dm.getValue(3,4),dm(3,4))
}
END_SECTION

START_SECTION((ValueType operator()(SizeType i, SizeType j)))
{
	TEST_EQUAL(dm.getValue(0,1),dm(0,1))
	TEST_EQUAL(dm.getValue(0,2),dm(0,2))
	TEST_EQUAL(dm.getValue(0,3),dm(0,3))
	TEST_EQUAL(dm.getValue(0,4),dm(0,4))
	TEST_EQUAL(dm.getValue(1,2),dm(1,2))
	TEST_EQUAL(dm.getValue(1,3),dm(1,3))
	TEST_EQUAL(dm.getValue(1,4),dm(1,4))
	TEST_EQUAL(dm.getValue(2,3),dm(2,3))
	TEST_EQUAL(dm.getValue(2,4),dm(2,4))
	TEST_EQUAL(dm.getValue(3,4),dm(3,4))
}
END_SECTION

START_SECTION((void reduce(SizeType j)))
{
	dm.reduce(2);
	TEST_EQUAL(dm.getValue(0,1),1)
	TEST_EQUAL(dm.getValue(0,2),3)
	TEST_EQUAL(dm.getValue(0,3),4)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),7)
	TEST_EQUAL(dm.getValue(2,3),10)
	TEST_EQUAL(dm.dimensionsize(),7)
}
END_SECTION

START_SECTION((std::pair<SizeType,SizeType> getMinElementCoordinates() const))
{
	dm.updateMinElement();
	pair<Size,Size> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,1)
	TEST_EQUAL(min.second,0)
}
END_SECTION

START_SECTION((void updateMinElement()))
{
	dm.setValueQuick(2,3,0.5);
	dm.updateMinElement();
	std::pair<Size,Size> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,3)
	TEST_EQUAL(min.second,2)
}
END_SECTION

DistanceMatrix<double> dm3(dm);

START_SECTION(bool operator==(DistanceMatrix< ValueType > const &rhs) const)
{
	TEST_EQUAL((dm==dm3),true)
}
END_SECTION


START_SECTION((template <typename Value> std::ostream & operator<<(std::ostream &os, const DistanceMatrix< Value > &matrix)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



