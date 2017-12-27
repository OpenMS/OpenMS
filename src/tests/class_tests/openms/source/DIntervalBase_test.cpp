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
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

START_TEST(DIntervalBase, "$Id$")

/////////////////////////////////////////////////////////////
	
//1D check
DIntervalBase<1>* ptr1 = nullptr;
DIntervalBase<1>* nullPointer1 = nullptr;

START_SECTION((DIntervalBase()))
	ptr1 = new DIntervalBase<1>;
  TEST_NOT_EQUAL(ptr1, nullPointer1)
END_SECTION

START_SECTION((~DIntervalBase()))
	delete ptr1;
END_SECTION

//2D check
DIntervalBase<2>* ptr2 = nullptr;
DIntervalBase<2>* nullPointer2 = nullptr;

START_SECTION([EXTRA] DIntervalBase())
	ptr2 = new DIntervalBase<2>;
  TEST_NOT_EQUAL(ptr2, nullPointer2)
END_SECTION

START_SECTION([EXTRA] ~DIntervalBase())
	delete ptr2;
END_SECTION

//misc stuff for testing
typedef DIntervalBase<2> I2;
typedef DIntervalBase<2>::PositionType I2Pos;
I2Pos p1;
p1[0]=5.0;
p1[1]=17.5;
I2Pos p2;
p2[0]=65.0;
p2[1]=-57.5;

START_SECTION((PositionType const& maxPosition() const))
  TEST_EQUAL(I2::empty.maxPosition()==I2Pos::minNegative(), true);
  TEST_EQUAL(I2::zero.maxPosition()==I2Pos::zero(), true);
END_SECTION

START_SECTION((PositionType const& minPosition() const))
  TEST_EQUAL(I2::empty.minPosition()==I2Pos::maxPositive(), true);
  TEST_EQUAL(I2::zero.minPosition()==I2Pos::zero(), true);
END_SECTION

START_SECTION((void setMinMax(PositionType const & min, PositionType const & max)))
  I2 tmp(I2::empty);
  tmp.setMinMax(p1,p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);
END_SECTION

START_SECTION((void setMin(PositionType const & position)))
  I2 tmp(I2::empty);
  tmp.setMin(p1);
  TEST_EQUAL(tmp.minPosition(),p1);
  TEST_EQUAL(tmp.maxPosition(),p1);
  tmp.setMin(p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);
END_SECTION

START_SECTION((void setMax(PositionType const & position)))
  I2 tmp(I2::empty);
  tmp.setMax(p1);
  TEST_EQUAL(tmp.minPosition(),p1);
  TEST_EQUAL(tmp.maxPosition(),p1);
  tmp.setMax(p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],-57.5);
END_SECTION

START_SECTION((bool operator==(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp==tmp,true);
	TEST_EQUAL(tmp==I2::empty,true);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp==I2::empty,false);
END_SECTION

START_SECTION((bool operator!=(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp!=tmp,false);
	TEST_EQUAL(tmp!=I2::empty,false);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp!=I2::empty,true);
END_SECTION

START_SECTION((DIntervalBase(const DIntervalBase& rhs)))
	I2 tmp(p1,p2);
	I2 tmp2(tmp);
	TEST_EQUAL(tmp==tmp2,true);
END_SECTION

START_SECTION((DIntervalBase( PositionType const & minimum, PositionType const & maximum )))
	I2 tmp(p1,p2);
	I2 tmp2(tmp.minPosition(), tmp.maxPosition());
	TEST_EQUAL(tmp==tmp2,true);
END_SECTION

START_SECTION((DIntervalBase& operator=(const DIntervalBase & rhs)))
	I2 tmp(p1,p2);
	I2 tmp2;
	TEST_EQUAL(tmp==tmp2,false);
	tmp2 = tmp;
	TEST_EQUAL(tmp==tmp2,true);
	tmp2 = tmp = I2::empty;
	TEST_EQUAL(tmp==tmp2,true);
	TEST_EQUAL(tmp==I2::empty,true);
END_SECTION

START_SECTION((void clear()))
	I2 tmp;
	TEST_EQUAL(tmp==I2::empty,true);
	tmp.setMax(p1);
	TEST_EQUAL(tmp==I2::empty,false);
	tmp.clear();
	TEST_EQUAL(tmp==I2::empty,true);
  TEST_EQUAL(tmp.maxPosition()==I2Pos::minNegative(), true);
  TEST_EQUAL(tmp.minPosition()==I2Pos::maxPositive(), true);
END_SECTION

START_SECTION((PositionType center() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.center());
  TEST_REAL_SIMILAR(pos[0],35.0);
  TEST_REAL_SIMILAR(pos[1],-20.0);
END_SECTION

START_SECTION((PositionType diagonal() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.diagonal());
  TEST_REAL_SIMILAR(pos[0],60.0);
  TEST_REAL_SIMILAR(pos[1],75.0);
END_SECTION

START_SECTION((CoordinateType width() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.width(),60.0)
END_SECTION

START_SECTION((CoordinateType height() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.height(),75.0)
END_SECTION

START_SECTION((CoordinateType maxX() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.maxX(),65.0)
END_SECTION

START_SECTION((CoordinateType maxY() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.maxY(),17.5)
END_SECTION

START_SECTION((CoordinateType minX() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.minX(),5.0)
END_SECTION

START_SECTION((CoordinateType minY() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.minY(),-57.5)
END_SECTION

START_SECTION((void setMinX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinX(57.67);
	TEST_REAL_SIMILAR(tmp.minX(),57.67)
END_SECTION

START_SECTION((void setMaxX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxX(57.67);
	TEST_REAL_SIMILAR(tmp.maxX(),57.67)
END_SECTION

START_SECTION((void setMinY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinY(57.67);
	TEST_REAL_SIMILAR(tmp.minY(),57.67)
END_SECTION

START_SECTION((void setMaxY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxY(57.67);
	TEST_REAL_SIMILAR(tmp.maxY(),57.67)
END_SECTION



START_SECTION((template <UInt D2> void assign(const DIntervalBase< D2 > rhs)))
DIntervalBase<2>::PositionType p1;
p1[0]=5.0;
p1[1]=17.5;
DIntervalBase<2>::PositionType p2;
p2[0]=65.0;
p2[1]=-57.5;
DIntervalBase<2> i2(p1,p2);

DIntervalBase<3> tmp;
tmp.assign(i2);
TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);

DIntervalBase<1> tmp2;
tmp2.assign(i2);
TEST_REAL_SIMILAR(tmp2.minPosition()[0],5.0);
TEST_REAL_SIMILAR(tmp2.maxPosition()[0],65.0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


