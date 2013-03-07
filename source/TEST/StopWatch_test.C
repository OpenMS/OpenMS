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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/SYSTEM/StopWatch.h>

#include <ctime>
/////////////////////////////////////////////////////////////

using namespace OpenMS;

void wait (int seconds)
{
   clock_t endwait = clock () + seconds * CLOCKS_PER_SEC;
   while (clock () < endwait) {}
}

START_TEST(StopWatch, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((StopWatch& operator = (const StopWatch& stop_watch)))
  StopWatch s1;
  s1.start();
  wait(1);
  s1.stop();

  StopWatch s2;
  TEST_EQUAL(s1!=s2, true)
  s2 = s1;
  TEST_EQUAL(s1==s2, true)

END_SECTION

START_SECTION((StopWatch()))
  NOT_TESTABLE; // tested above
END_SECTION

START_SECTION((StopWatch(const StopWatch& stop_watch)))
  StopWatch s1;
  s1.start();
  wait(1);
  s1.stop();

  StopWatch s2(s1);
  TEST_EQUAL(s1==s2, true)
  
END_SECTION

START_SECTION((bool isRunning() const))
  StopWatch w;

  w.start();
  TEST_EQUAL(w.isRunning(), true);

  w.stop();

END_SECTION

START_SECTION((bool operator != (const StopWatch& stop_watch) const))
  StopWatch s, s2;
  StopWatch s3(s2);
  TEST_EQUAL(s2==s3, true);


  s.start();
  s2.start();
  wait(3);
  s.stop();

  wait(3);
  s2.stop();

  TEST_EQUAL(s!=s2, true)

  TEST_EQUAL(s<=s2, true)

  TEST_EQUAL(s2>=s, true)

  TEST_EQUAL(s2!=s3, true)
  s3 = s2;
  TEST_EQUAL(s2==s3, true)

  s2.start();
  wait(1);
  TEST_EQUAL(s2==s3, false)
  s2.stop();

END_SECTION

START_SECTION((bool operator < (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // since we do not have control over system time...
END_SECTION

START_SECTION((bool operator <= (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // tested above
END_SECTION

START_SECTION((bool operator == (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // tested above
END_SECTION

START_SECTION((bool operator > (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // since we do not have control over system time...
END_SECTION

START_SECTION((bool operator >= (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // tested above
END_SECTION

START_SECTION((bool start()))
  NOT_TESTABLE; // see below
END_SECTION

START_SECTION((bool stop()))
  StopWatch s;
  s.start();
  wait(3);
  s.stop();

  TEST_EQUAL(s.getClockTime() > 2, true)
  TEST_EQUAL(s.getClockTime() < 4, true)
  
  DoubleReal t1 = s.getCPUTime();
  DoubleReal t2 = s.getClockTime();
  DoubleReal t3 = s.getSystemTime();
  DoubleReal t4 = s.getUserTime();
  // wait some more
  wait(3);
  // ... and see if time is still the old one
  TEST_EQUAL(s.getCPUTime(), t1)
  TEST_EQUAL(s.getClockTime(), t2)
  TEST_EQUAL(s.getSystemTime(), t3)
  TEST_EQUAL(s.getUserTime(), t4)

END_SECTION

START_SECTION((DoubleReal getCPUTime() const ))
  StopWatch s;
  s.start();
  wait(3);
  s.stop();

  //std::cerr << "CPU TIME: " << s.getCPUTime() << " "<< s.getUserTime() << "  "<< s.getSystemTime() << "\n\n";

  TEST_EQUAL(s.getCPUTime() > 0.1, true) // waiting costs CPU time... just not sure how much...
  TEST_EQUAL(s.getClockTime() > 2, true) // and must consume wall time
  TEST_EQUAL(s.getClockTime() < 4, true)
  TEST_EQUAL(s.getUserTime() > 0.1, true) //  and some user time
  TEST_EQUAL(s.getUserTime() < 4, true)
  TEST_EQUAL(s.getSystemTime() < 0.5, true) // but no system time
END_SECTION

START_SECTION((DoubleReal getClockTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((DoubleReal getSystemTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((DoubleReal getUserTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((void clear()))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((void reset()))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((virtual ~StopWatch()))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((static String toString(DoubleReal time)))
  
  TEST_EQUAL(StopWatch::toString(0), "0.00 s")

  TEST_EQUAL(StopWatch::toString(1), "1.00 s")
  TEST_EQUAL(StopWatch::toString(1.5), "1.50 s")
  TEST_EQUAL(StopWatch::toString(100.5), "01:40 m")
  TEST_EQUAL(StopWatch::toString(3600*24*5 + 3600*9 + 5), "5d 09:00:05 h")
  TEST_EQUAL(StopWatch::toString(160.5), "02:40 m")
  TEST_EQUAL(StopWatch::toString(3600*23 + 160.5), "23:02:40 h")


END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
