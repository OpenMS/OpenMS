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
// $Authors: $
// --------------------------------------------------------------------------
//


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

// Includes go here....
#include <OpenMS/DATASTRUCTURES/Matrix.h>

///////////////////////////

START_TEST(Matrix, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Matrix<int>* ptr = nullptr;
Matrix<int>* nullPointer = nullptr;
START_SECTION((Matrix()))
{
	ptr = new Matrix<int>;
  TEST_NOT_EQUAL(ptr, nullPointer);

  Matrix<int> mi1;
	TEST_EQUAL(mi1.size(),0);
	TEST_EQUAL(mi1.cols(),0);
	TEST_EQUAL(mi1.rows(),0);
	TEST_EQUAL(mi1.empty(),true);
  STATUS("mi1:\n"<< mi1);
}
END_SECTION;

START_SECTION((~Matrix()))
{
	delete ptr;
}
END_SECTION;

Matrix<int> mi;

START_SECTION((void resize(size_type i, size_type j, value_type value = value_type())))
{
	mi.resize(2,2,3);
  STATUS("mi1:\n"<< mi);
	mi.resize(2,3,7);
  STATUS("mi1:\n"<< mi);
	TEST_EQUAL(mi(0,0),3);
	TEST_EQUAL(mi(0,1),3);
	TEST_EQUAL(mi(0,2),3);
	TEST_EQUAL(mi(1,0),3);
	TEST_EQUAL(mi(1,1),7);
	TEST_EQUAL(mi(1,2),7);
}
END_SECTION

START_SECTION((void resize(std::pair<Size, Size> const & size_pair, value_type value = value_type())))
{
	std::pair<Size, Size> const sizepair(2,2);
	mi.resize(sizepair,3);
  STATUS("mi1:\n"<< mi);
	mi.resize(2,3,7);
  STATUS("mi1:\n"<< mi);
	TEST_EQUAL(mi(0,0),3);
	TEST_EQUAL(mi(0,1),3);
	TEST_EQUAL(mi(0,2),3);
	TEST_EQUAL(mi(1,0),3);
	TEST_EQUAL(mi(1,1),7);
	TEST_EQUAL(mi(1,2),7);
}
END_SECTION

START_SECTION((Matrix(const Matrix & source)))
{
  Matrix<int> mi2(mi);
  STATUS("mi2:\n"<< mi2);
	TEST_EQUAL(mi2.cols(),3);
	TEST_EQUAL(mi2.rows(),2);
	TEST_EQUAL(mi2(0,0),3);
	TEST_EQUAL(mi2(0,1),3);
	TEST_EQUAL(mi2(0,2),3);
	TEST_EQUAL(mi2(1,0),3);
	TEST_EQUAL(mi2(1,1),7);
	TEST_EQUAL(mi2(1,2),7);
}
END_SECTION

START_SECTION((Matrix& operator = (const Matrix & rhs)))
{
	Matrix<int> mi3;
  STATUS("mi3:\n"<<mi3);
	mi3=mi;
  STATUS("mi3:\n"<<mi3);
	TEST_EQUAL(mi3.cols(),3);
	TEST_EQUAL(mi3.rows(),2);
	TEST_EQUAL(mi3(0,0),3);
	TEST_EQUAL(mi3(0,1),3);
	TEST_EQUAL(mi3(0,2),3);
	TEST_EQUAL(mi3(1,0),3);
	TEST_EQUAL(mi3(1,1),7);
	TEST_EQUAL(mi3(1,2),7);
}
END_SECTION

mi(1,1)=17;

START_SECTION((const_reference getValue(size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr.getValue(1,1),17);
}
END_SECTION

START_SECTION((const_reference operator() (size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr(1,1),17);
}
END_SECTION

START_SECTION((reference getValue(size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,2));
	mi.getValue(1,2)=33;
	STATUS(mi.getValue(1,2));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,2),33);
}
END_SECTION

START_SECTION((reference operator() (size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,0));
	mi.getValue(1,0) = 44;
	STATUS(mi.getValue(1,0));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,0), 44);
}
END_SECTION

START_SECTION(container_type row(size_type const i) const)
{
	Matrix<int>::container_type row = mi.row(0);
	TEST_EQUAL(row.size(), 3)
	TEST_EQUAL(row[0], 3)
	TEST_EQUAL(row[1], 3)
	TEST_EQUAL(row[2], 3)
	row = mi.row(1);
	TEST_EQUAL(row[0], 44)
	TEST_EQUAL(row[1], 17)
	TEST_EQUAL(row[2], 33)
}
END_SECTION

START_SECTION(container_type col(size_type const i) const)
{
	Matrix<int>::container_type col = mi.col(0);
	TEST_EQUAL(col.size(), 2)
	TEST_EQUAL(col[0], 3)
	TEST_EQUAL(col[1], 44)
	col = mi.col(1);
	TEST_EQUAL(col.size(), 2)
	TEST_EQUAL(col[0], 3)
	TEST_EQUAL(col[1], 17)
	col = mi.col(2);
	TEST_EQUAL(col.size(), 2)
	TEST_EQUAL(col[0], 3)
	TEST_EQUAL(col[1], 33)
}
END_SECTION

START_SECTION((void clear()))
{
	Matrix<int> mi4(mi);
  STATUS("mi4:\n"<<mi4);
	mi4.clear();
  STATUS("mi4:\n"<<mi4);
  TEST_EQUAL(mi4.empty(),true);
}
END_SECTION

START_SECTION((void setValue(size_type const i, size_type const j, value_type value)))
{
	mi.setValue(1,1,18);
	STATUS("mi:\n"<<mi);
	TEST_EQUAL(mi(1,1),18);
}
END_SECTION;

Matrix<int> mi5(4,5,6);

START_SECTION((Matrix(const SizeType rows, const SizeType cols, ValueType value = ValueType())))
{
  STATUS("mi5:\n"<<mi5);
	TEST_EQUAL(mi5.size(),20);
}
END_SECTION;

START_SECTION((SizeType cols() const))
{
	TEST_EQUAL(mi5.rows(),4);
}
END_SECTION;

START_SECTION((SizeType rows() const))
{
	TEST_EQUAL(mi5.cols(),5);
}
END_SECTION;

Matrix<float> mf(6,7,8);

START_SECTION((SizeType colIndex(SizeType index) const))
{
	TEST_EQUAL(mf.colIndex(30),2);
}
END_SECTION;

START_SECTION((SizeType const index(SizeType row, SizeType col) const))
{
	TEST_EQUAL(mf.index(5,5),40);
}
END_SECTION;

START_SECTION((SizeType rowIndex(SizeType index) const))
{
  TEST_EQUAL(mf.rowIndex(30),4);
}
END_SECTION;

START_SECTION((std::pair<Size,Size> const indexPair(Size index) const))
{
	std::pair<Size,Size> result = mf.indexPair(30);
  TEST_EQUAL(result.first,4);
	TEST_EQUAL(result.second,2);
}
END_SECTION

START_SECTION((std::pair<Size,Size> sizePair() const))
{
	Matrix<float> const mf(6,7,8);
	TEST_EQUAL(mf.sizePair().first,6);
	TEST_EQUAL(mf.sizePair().second,7);
}
END_SECTION

START_SECTION((bool operator==(Matrix const &rhs) const))
{
	Matrix<int> mi1(4,5,6);
	mi1(2,3)=17;
	Matrix<int> mi2(4,5,6);
	TEST_NOT_EQUAL(mi1,mi2);
	mi1(2,3)=6;
	TEST_EQUAL(mi1,mi2);

	Matrix<int> mi3(5,4,6);
	Matrix<int> mi4(4,4,6);
	Matrix<int> mi5(5,5,6);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi3);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi4);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi5);(void) comparison);
}
END_SECTION

START_SECTION((bool operator<(Matrix const &rhs) const))
{
	Matrix<int> mi1(4,5,6);
	TEST_EQUAL(mi1<mi1,false);
	mi1(2,3)=17;
	TEST_EQUAL(mi1<mi1,false);
	Matrix<int> mi2(4,5,6);
	TEST_EQUAL(mi1<mi2,false);
	TEST_EQUAL(mi2<mi1,true);
	mi2(2,3)=18;
	TEST_EQUAL(mi1<mi2,true);

	Matrix<int> mi3(5,4,6);
	Matrix<int> mi4(4,4,6);
	Matrix<int> mi5(5,5,6);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi3);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi4);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi5);(void) comparison);
}
END_SECTION

START_SECTION((template <int ROWS, int COLS> void setMatrix(const ValueType matrix[ROWS][COLS])))
{
	double test_matrix[4][4] = {
		{0, 2.5,   3, 0.1},
		{0,   1, 5.9, 0.2},
		{0,   2, 5.6, 0.1},
		{0,   2,   3, 0.1}
	};

	Matrix<double> myMatrix;
	myMatrix.setMatrix<4,4>(test_matrix);
	for (size_t i=0;i<4;++i)
	{
		for (size_t j=0;j<4;++j)
		{
			TEST_EQUAL( myMatrix(i,j), test_matrix[i][j] )
		}
	}

}
END_SECTION

START_SECTION((template <typename  Value > std::ostream & operator<<(std::ostream &os, const Matrix< Value > &matrix)))
{
	Matrix<int> mi(2,3,6);
	mi(1,2)=112;
	mi(0,0)=100;
	mi(1,1)=111;
	mi(0,2)=103;
	std::ostringstream os;
	os << mi;
	// Uh, finally I got the whitespace right
	char matrix_dump[] =
	"   100      6    103 \n"
	"     6    111    112 \n";
	TEST_EQUAL(os.str(),matrix_dump);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


