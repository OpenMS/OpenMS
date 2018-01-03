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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/MultiGradient.h>
#include <QtGui/QColor>
#include <OpenMS/CONCEPT/Types.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MultiGradient, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MultiGradient* d10_ptr = nullptr;
MultiGradient* d10_nullPointer = nullptr;
START_SECTION((MultiGradient()))
	d10_ptr = new MultiGradient();
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~MultiGradient()))
	delete d10_ptr;
END_SECTION

START_SECTION((InterpolationMode getInterpolationMode() const))
	TEST_EQUAL(MultiGradient().getInterpolationMode(),MultiGradient::IM_LINEAR)
END_SECTION

START_SECTION((void setInterpolationMode(InterpolationMode mode)))
	MultiGradient mg;
	mg.setInterpolationMode(MultiGradient::IM_STAIRS);
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_STAIRS)
END_SECTION

START_SECTION((Size size() const))
	MultiGradient mg;
	TEST_EQUAL(mg.size(),2);
END_SECTION

START_SECTION((UInt position(UInt index)))
	MultiGradient mg;
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),100);
END_SECTION

START_SECTION((QColor color(UInt index)))
	MultiGradient mg;
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::black,true);
END_SECTION

START_SECTION((void insert(double position, QColor color)))
	MultiGradient mg;
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),50);
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::red,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),50);
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::red,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.insert(25,Qt::green);
	mg.insert(75,Qt::blue);
	TEST_EQUAL(mg.size(),5);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);
	TEST_EQUAL(mg.position(2),50);
	TEST_EQUAL(mg.position(3),75);
	TEST_EQUAL(mg.position(4),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::red,true);
	TEST_EQUAL(mg.color(3)==Qt::blue,true);
	TEST_EQUAL(mg.color(4)==Qt::black,true);
	mg.insert(76,Qt::magenta);
	TEST_EQUAL(mg.size(),6);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);
	TEST_EQUAL(mg.position(2),50);
	TEST_EQUAL(mg.position(3),75);
	TEST_EQUAL(mg.position(4),76);
	TEST_EQUAL(mg.position(5),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::red,true);
	TEST_EQUAL(mg.color(3)==Qt::blue,true);
	TEST_EQUAL(mg.color(4)==Qt::magenta,true);
	TEST_EQUAL(mg.color(5)==Qt::black,true);
END_SECTION

START_SECTION((bool remove(double position)))
	MultiGradient mg;
	mg.insert(25,Qt::green);
	mg.insert(50,Qt::red);
	mg.insert(75,Qt::blue);
	mg.remove(50);
	TEST_EQUAL(mg.size(),4);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),25);
	TEST_EQUAL(mg.position(2),75);
	TEST_EQUAL(mg.position(3),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::green,true);
	TEST_EQUAL(mg.color(2)==Qt::blue,true);
	TEST_EQUAL(mg.color(3)==Qt::black,true);
	mg.remove(25);
	TEST_EQUAL(mg.size(),3);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),75);
	TEST_EQUAL(mg.position(2),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::blue,true);
	TEST_EQUAL(mg.color(2)==Qt::black,true);
	mg.remove(75);
	TEST_EQUAL(mg.size(),2);
	TEST_EQUAL(mg.position(0),0);
	TEST_EQUAL(mg.position(1),100);
	TEST_EQUAL(mg.color(0)==Qt::white,true);
	TEST_EQUAL(mg.color(1)==Qt::black,true);
END_SECTION

START_SECTION((bool exists(double position)))
	MultiGradient mg;
	mg.insert(25,Qt::green);
	mg.insert(50,Qt::red);
	mg.insert(75,Qt::blue);
	TEST_EQUAL(mg.exists(0),true);
	TEST_EQUAL(mg.exists(1),false);
	TEST_EQUAL(mg.exists(25),true);
	TEST_EQUAL(mg.exists(49),false);
	TEST_EQUAL(mg.exists(50),true);
	TEST_EQUAL(mg.exists(51),false);
	TEST_EQUAL(mg.exists(75),true);
	TEST_EQUAL(mg.exists(99),false);
	TEST_EQUAL(mg.exists(100),true);
END_SECTION

START_SECTION((QColor interpolatedColorAt(double position) const))
	MultiGradient mg;
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(191,191,191),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==QColor(127,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(63,63,63),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(255,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(127,0,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);
	mg.insert(50,Qt::green);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(127,255,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::green,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(0,127,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);
	mg.insert(50,Qt::blue);
	TEST_EQUAL(mg.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25)==QColor(127,127,255),true);
	TEST_EQUAL(mg.interpolatedColorAt(50)==Qt::blue,true);
	TEST_EQUAL(mg.interpolatedColorAt(75)==QColor(0,0,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(100)==Qt::black,true);

	MultiGradient mg2;
	mg2.setInterpolationMode(MultiGradient::IM_STAIRS);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(25)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);
	mg2.insert(50,Qt::red);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(49)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(51)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(99)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);
END_SECTION

START_SECTION((QColor interpolatedColorAt(double position, double min, double max) const))
	MultiGradient mg;
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.interpolatedColorAt(0,0,100)==Qt::white,true);
	TEST_EQUAL(mg.interpolatedColorAt(25,0,100)==QColor(255,127,127),true);
	TEST_EQUAL(mg.interpolatedColorAt(50,0,100)==Qt::red,true);
	TEST_EQUAL(mg.interpolatedColorAt(75,0,100)==QColor(127,0,0),true);
	TEST_EQUAL(mg.interpolatedColorAt(100,0,100)==Qt::black,true);

	MultiGradient mg2;
	mg2.setInterpolationMode(MultiGradient::IM_STAIRS);
	mg2.insert(50,Qt::red);
	TEST_EQUAL(mg2.interpolatedColorAt(0)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(49)==Qt::white,true);
	TEST_EQUAL(mg2.interpolatedColorAt(50)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(51)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(99)==Qt::red,true);
	TEST_EQUAL(mg2.interpolatedColorAt(100)==Qt::black,true);
END_SECTION

START_SECTION((void activatePrecalculationMode(double min, double max, UInt steps)))
NOT_TESTABLE
END_SECTION

START_SECTION((QColor precalculatedColorAt(double position) const ))
	MultiGradient mg;
	mg.insert(0,Qt::white);
	mg.insert(100,Qt::blue);
	mg.activatePrecalculationMode(-50.0,50.0,100);

	//Test precalclulated Values
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).red(),255);
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).green(),255);
	TEST_EQUAL(mg.precalculatedColorAt(-50.0).blue(),255);

	TEST_EQUAL(mg.precalculatedColorAt(-25.0).red(),193);
	TEST_EQUAL(mg.precalculatedColorAt(-25.0).green(),193);
	TEST_EQUAL(mg.precalculatedColorAt(-25.0).blue(),255);

	TEST_EQUAL(mg.precalculatedColorAt(0.0).red(),128);
	TEST_EQUAL(mg.precalculatedColorAt(0.0).green(),128);
	TEST_EQUAL(mg.precalculatedColorAt(0.0).blue(),255);

	TEST_EQUAL(mg.precalculatedColorAt(25.0).red(),64);
	TEST_EQUAL(mg.precalculatedColorAt(25.0).green(),64);
	TEST_EQUAL(mg.precalculatedColorAt(25.0).blue(),255);

	TEST_EQUAL(mg.precalculatedColorAt(50.0).red(),2);
	TEST_EQUAL(mg.precalculatedColorAt(50.0).green(),2);
	TEST_EQUAL(mg.precalculatedColorAt(50.0).blue(),255);
END_SECTION

START_SECTION((void deactivatePrecalculationMode()))
	MultiGradient mg;
	mg.activatePrecalculationMode(-50,50,100);
	mg.deactivatePrecalculationMode();
	NOT_TESTABLE
END_SECTION

START_SECTION((std::string toString() const))
	MultiGradient mg;
	TEST_EQUAL(mg.toString(),"Linear|0,#ffffff;100,#000000")
	mg.setInterpolationMode(MultiGradient::IM_STAIRS);
	mg.insert(50,Qt::red);
	TEST_EQUAL(mg.toString(),"Stairs|0,#ffffff;50,#ff0000;100,#000000")
END_SECTION

START_SECTION((void fromString(const std::string& gradient)))
	MultiGradient mg;
	mg.fromString("Linear|0,#ff0000;100,#000000");
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_LINEAR)
	TEST_EQUAL(mg.size(),2)
	TEST_EQUAL(mg.color(0)==Qt::red, true);
	TEST_EQUAL(mg.color(1)==Qt::black, true);
	TEST_EQUAL(mg.position(0), 0);
	TEST_EQUAL(mg.position(1), 100);
	mg.fromString("Stairs|0,#ffffff;50,#ff0000;100,#000000");
	TEST_EQUAL(mg.getInterpolationMode(),MultiGradient::IM_STAIRS)
	TEST_EQUAL(mg.size(),3)
	TEST_EQUAL(mg.color(0)==Qt::white, true);
	TEST_EQUAL(mg.color(1)==Qt::red, true);
	TEST_EQUAL(mg.color(2)==Qt::black, true);
	TEST_EQUAL(mg.position(0), 0);
	TEST_EQUAL(mg.position(1), 50);
	TEST_EQUAL(mg.position(2), 100);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



