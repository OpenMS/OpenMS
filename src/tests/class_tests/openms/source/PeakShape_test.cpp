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
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

///////////////////////////

START_TEST(PeakShape, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

PeakShape* peakshape_ptr = nullptr;
PeakShape* peakshape_nullPointer = nullptr;

START_SECTION((PeakShape()))
  peakshape_ptr = new PeakShape;
  TEST_NOT_EQUAL(peakshape_ptr, peakshape_nullPointer)
END_SECTION

START_SECTION((virtual ~PeakShape()))
		delete peakshape_ptr;
END_SECTION
	
START_SECTION((PeakShape& operator = (const PeakShape& rhs)))
		PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy;
    peakshape_copy = peakshape;

    TEST_REAL_SIMILAR(peakshape_copy.height, 10003.232) 
    TEST_REAL_SIMILAR(peakshape_copy.mz_position, 0.323)
    TEST_REAL_SIMILAR(peakshape_copy.left_width, 2.998)
    TEST_REAL_SIMILAR(peakshape_copy.right_width, 2.776)
    TEST_REAL_SIMILAR(peakshape_copy.area, 8329832.141)
    TEST_EQUAL(peakshape_copy.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((PeakShape(const PeakShape& rhs)))
    PeakShape peakshape;
    peakshape.height = 10003.232;
    peakshape.mz_position = 0.323;
    peakshape.left_width = 2.998;
    peakshape.right_width = 2.776;
    peakshape.area = 8329832.141;
    peakshape.type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape_copy(peakshape);
   
    TEST_REAL_SIMILAR(peakshape.height,10003.232) 
    TEST_REAL_SIMILAR(peakshape.mz_position, 0.323)
    TEST_REAL_SIMILAR(peakshape.left_width, 2.998)
    TEST_REAL_SIMILAR(peakshape.right_width,2.776)
    TEST_REAL_SIMILAR(peakshape.area, 8329832.141)
    TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION
MSSpectrum spec;
spec.resize(100);
for(Int i = 0; i<100;++i)
{
  spec[i].setMZ(i*0.1);
  spec[i].setIntensity(100.); 
}
START_SECTION((PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakIterator left_, PeakIterator right_, Type type_)))
    double height = 100.0;
    double mz_position = 0.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    
    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_EQUAL(peakshape.iteratorsSet(), true)
    TEST_REAL_SIMILAR(peakshape.height,height) 
	  TEST_REAL_SIMILAR(peakshape.mz_position, mz_position)
		TEST_REAL_SIMILAR(peakshape.left_width, left_width)
    TEST_REAL_SIMILAR(peakshape.right_width, right_width)
    TEST_REAL_SIMILAR(peakshape.area, area)
    TEST_REAL_SIMILAR(peakshape.r_value, 0.0)
		TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, Type type_)))
    double height = 100.0;
    double mz_position = 0.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);

    TEST_EQUAL(peakshape.iteratorsSet(), false)
    TEST_REAL_SIMILAR(peakshape.height,height) 
	  TEST_REAL_SIMILAR(peakshape.mz_position, mz_position)
		TEST_REAL_SIMILAR(peakshape.left_width, left_width)
    TEST_REAL_SIMILAR(peakshape.right_width, right_width)
    TEST_REAL_SIMILAR(peakshape.area, area)
    TEST_REAL_SIMILAR(peakshape.r_value, 0.0)
		TEST_EQUAL(peakshape.type, PeakShape::LORENTZ_PEAK)
END_SECTION

START_SECTION((bool iteratorsSet() const))
    double height = 100.0;
    double mz_position = 0.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);

    
    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape2(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_EQUAL(peakshape2.iteratorsSet(), true)
    TEST_EQUAL(peakshape.iteratorsSet(), false)
END_SECTION

START_SECTION((PeakIterator getRightEndpoint() const))
    double height = 100.0;
    double mz_position = 4.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getMZ(), (spec.begin()+30)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getIntensity(), (spec.begin()+30)->getIntensity())
END_SECTION


START_SECTION((void setRightEndpoint(PeakIterator right_endpoint)))
    double height = 100.0;
    double mz_position = 4.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);

    peakshape.setLeftEndpoint(it1);
    peakshape.setRightEndpoint(it2); 
    TEST_EQUAL(peakshape.iteratorsSet(), true)
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getMZ(), (spec.begin()+30)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getRightEndpoint()->getIntensity(), (spec.begin()+30)->getIntensity())
END_SECTION

START_SECTION((PeakIterator getLeftEndpoint() const))
    double height = 100.0;
    double mz_position = 4.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape::PeakIterator it2 = spec.begin()+30;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,it1,it2,
		type);

    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getMZ(), (spec.begin()+2)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getIntensity(), (spec.begin()+2)->getIntensity())
END_SECTION


START_SECTION((void setLeftEndpoint(PeakIterator left_endpoint)))
    double height = 100.0;
    double mz_position = 4.0;
    double left_width = 3.0;
    double right_width = 3.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;

    PeakShape::PeakIterator it1 = spec.begin()+2;
    PeakShape peakshape(height,
    mz_position,
    left_width,
    right_width,
		area,
		type);
    peakshape.setLeftEndpoint(it1);

    TEST_EQUAL(peakshape.iteratorsSet(), false)
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getMZ(), (spec.begin()+2)->getMZ())
    TEST_REAL_SIMILAR(peakshape.getLeftEndpoint()->getIntensity(), (spec.begin()+2)->getIntensity())
END_SECTION


START_SECTION((double getSymmetricMeasure() const))
		double height = 100.0;
    double mz_position = 0.0;
    double left_width = 3.0;
    double right_width = 9.0;
    double area = 309.23292;
    PeakShape::Type type = PeakShape::SECH_PEAK;

    PeakShape peakshape(height,
												mz_position,
												left_width,
												right_width,
												area,
												type);

    double sym_value = peakshape.getSymmetricMeasure();
    TEST_REAL_SIMILAR(sym_value,3.0/9.0)
END_SECTION

START_SECTION((double operator() (double x) const))
    double height = 100.0;
    double mz_position = 0.0;
    double left_width = 4.0;
    double right_width = 4.0;
    double area = 100;
    PeakShape::Type type = PeakShape::LORENTZ_PEAK;
    
    PeakShape peakshape(height,
												mz_position,
												left_width,
												right_width,
												area,
												type);
   
    TEST_REAL_SIMILAR(peakshape.getFWHM(),.5)
END_SECTION

START_SECTION((double getFWHM() const))
  double height = 100.0;
  double mz_position = 0.0;
  double left_width = 4.0;
  double right_width = 4.0;
  double area = 100;
  PeakShape::Type type = PeakShape::LORENTZ_PEAK;
    
  PeakShape p(height,
							mz_position,
							left_width,
							right_width,
							area,
							type);


  TEST_REAL_SIMILAR(p.getFWHM(),1/right_width + 1/left_width)
END_SECTION

START_SECTION(bool operator==(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1==p2,true)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1==p2,false)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1==p2,true)
END_SECTION

START_SECTION(bool operator!=(const PeakShape &rhs) const)
	PeakShape p1,p2;
	TEST_EQUAL(p1!=p2,false)
	
	p1.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,true)
	
	p2.mz_position = 14.4;
	TEST_EQUAL(p1!=p2,false)
END_SECTION

START_SECTION(([PeakShape::PositionLess] bool operator()(const PeakShape &a, const PeakShape &b)))
  PeakShape p1(0.,123.,0.,0.,0.,PeakShape::LORENTZ_PEAK);
  PeakShape p2(0.,124.,0.,0.,0.,PeakShape::LORENTZ_PEAK);
  PeakShape::PositionLess comp;
  TEST_EQUAL(comp(p1,p2),true);
  TEST_EQUAL(comp(p2,p1),false);
END_SECTION

  
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

