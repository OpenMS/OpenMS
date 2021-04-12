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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

template<class T>
ostream& operator<<(ostream& os, const vector<T>& vec)
{
	if (vec.empty())
	{
		os << "()";
		return os;
	}

	os << "(";
	typename vector<T>::const_iterator i = vec.begin();
	while (true)
	{
		os << *i++;
		if (i == vec.end())
			break;
		os << ",";
	}
	os << ")";
	return os;
}

DRange<1> makeRange(double a, double b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

START_TEST(PeakFileOptions, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakFileOptions* ptr = nullptr;
PeakFileOptions* nullPointer = nullptr;
START_SECTION((PeakFileOptions()))
	ptr = new PeakFileOptions();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~PeakFileOptions()))
	delete ptr;
END_SECTION

START_SECTION((void setCompression(bool compress)))
	PeakFileOptions tmp;
	tmp.setCompression(true);
	TEST_EQUAL(tmp.getCompression(), true);
END_SECTION

START_SECTION((bool getCompression() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getCompression(), false);
END_SECTION

START_SECTION((void setMetadataOnly(bool only)))
	PeakFileOptions tmp;
	tmp.setMetadataOnly(true);
	TEST_EQUAL(tmp.getMetadataOnly(), true);
END_SECTION

START_SECTION((bool getMetadataOnly() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMetadataOnly(), false);
END_SECTION

START_SECTION((void setWriteSupplementalData(bool write)))
	PeakFileOptions tmp;
	tmp.setWriteSupplementalData(false);
	TEST_EQUAL(tmp.getWriteSupplementalData(), false);
END_SECTION

START_SECTION((bool getWriteSupplementalData() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getWriteSupplementalData(), true);
END_SECTION

START_SECTION((void setRTRange(const DRange<1>& range)))
	PeakFileOptions tmp;
	tmp.setRTRange(makeRange(2, 4));
	TEST_EQUAL(tmp.hasRTRange(), true);
	TEST_EQUAL(tmp.getRTRange(), makeRange(2, 4));
END_SECTION

START_SECTION((bool hasRTRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasRTRange(), false);
END_SECTION

START_SECTION((const DRange<1>& getRTRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getRTRange(), DRange<1>());
END_SECTION

START_SECTION((void setMZRange(const DRange<1>& range)))
	PeakFileOptions tmp;
	tmp.setMZRange(makeRange(3, 5));
	TEST_EQUAL(tmp.hasMZRange(), true);
	TEST_EQUAL(tmp.getMZRange(), makeRange(3, 5));
END_SECTION

START_SECTION((bool hasMZRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasMZRange(), false);
END_SECTION

START_SECTION((const DRange<1>& getMZRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMZRange(), DRange<1>());
END_SECTION

START_SECTION((void setIntensityRange(const DRange<1>& range)))
	PeakFileOptions tmp;
	tmp.setIntensityRange(makeRange(3, 5));
	TEST_EQUAL(tmp.hasIntensityRange(), true);
	TEST_EQUAL(tmp.getIntensityRange(), makeRange(3, 5));
END_SECTION

START_SECTION((bool hasIntensityRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasIntensityRange(), false);
END_SECTION

START_SECTION((const DRange<1>& getIntensityRange() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getIntensityRange(), DRange<1>());
END_SECTION

START_SECTION((void setMSLevels(const vector<Int>& levels)))
	PeakFileOptions tmp;
	vector<Int> levels;
	levels.push_back(1);
	levels.push_back(3);
	levels.push_back(5);
	tmp.setMSLevels(levels);
	TEST_EQUAL(tmp.hasMSLevels(), true);
	TEST_EQUAL(tmp.getMSLevels()==levels,true);
END_SECTION

START_SECTION((void addMSLevel(int level)))
	PeakFileOptions tmp;
	tmp.addMSLevel(1);
	tmp.addMSLevel(3);
	tmp.addMSLevel(5);
	TEST_EQUAL(tmp.hasMSLevels(), true);
	TEST_EQUAL(tmp.getMSLevels().size(), 3);
	vector<Int> levels;
	levels.push_back(1);
	levels.push_back(3);
	levels.push_back(5);
	TEST_EQUAL(tmp.getMSLevels()==levels,true);
END_SECTION

START_SECTION((void clearMSLevels()))
	PeakFileOptions tmp;
	vector<Int> levels;
	levels.push_back(1);
	levels.push_back(3);
	levels.push_back(5);
	tmp.setMSLevels(levels);
	TEST_EQUAL(tmp.getMSLevels()==levels,true);

	// now clear the ms levels
	tmp.clearMSLevels();
	TEST_EQUAL(tmp.hasMSLevels(), false);
	TEST_EQUAL(tmp.getMSLevels()==vector<Int>(),true);
END_SECTION

START_SECTION((bool hasMSLevels() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.hasMSLevels(), false);
END_SECTION

START_SECTION((bool containsMSLevel(int level) const))
	PeakFileOptions tmp;
	vector<Int> levels;
	levels.push_back(1);
	levels.push_back(3);
	levels.push_back(5);
	tmp.setMSLevels(levels);
	TEST_EQUAL(tmp.containsMSLevel(3), true);
	TEST_EQUAL(tmp.containsMSLevel(2), false);
END_SECTION

START_SECTION((const vector<Int>& getMSLevels() const))
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMSLevels()==vector<Int>(),true);
END_SECTION

START_SECTION(Size getMaxDataPoolSize() const)
{
	PeakFileOptions tmp;
	TEST_EQUAL(tmp.getMaxDataPoolSize()!=0,true);
}
END_SECTION

START_SECTION(void setMaxDataPoolSize(Size size))
{
	PeakFileOptions tmp;
	tmp.setMaxDataPoolSize(250);
	TEST_EQUAL(tmp.getMaxDataPoolSize()==250,true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
