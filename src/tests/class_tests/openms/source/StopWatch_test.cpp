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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/SYSTEM/StopWatch.h>

#include <chrono>
/////////////////////////////////////////////////////////////

using namespace OpenMS;

void wait(double seconds)
{
  auto start = std::chrono::system_clock::now();
  while (true)
  {
   double s = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
   if (s > seconds) break;
  };
}

START_TEST(StopWatch, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((StopWatch& operator = (const StopWatch& stop_watch)))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION((StopWatch()))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION((StopWatch(const StopWatch& stop_watch)))
  StopWatch s1, s2;
  s1.start();
  wait(0.01);
  TEST_EQUAL(s1 != s2, true) // before stop
  s1.stop();
  TEST_EQUAL(s1 != s2, true)
  s2 = s1;
  TEST_EQUAL(s1 == s2, true)

  StopWatch s3(s1);
  TEST_EQUAL(s1==s3, true)
  
  StopWatch s4;
  s1.reset();
  TEST_EQUAL(s1 == s4, true)

  s1.start();
  s2.start();
  wait(0.01);
  s1.stop();

  wait(0.01);
  s2.stop();

  TEST_EQUAL(s1 != s2, true)
  TEST_EQUAL(s1 <= s2, true)
  TEST_EQUAL(s2 >= s1, true)
END_SECTION

START_SECTION((bool isRunning() const))
  StopWatch w;
  TEST_EQUAL(w.isRunning(), false);
  w.start();
  TEST_EQUAL(w.isRunning(), true);
  w.stop();
  TEST_EQUAL(w.isRunning(), false); 
END_SECTION

START_SECTION((bool operator != (const StopWatch& stop_watch) const))
  NOT_TESTABLE; // tested above
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
  StopWatch s1;
  s1.start();
  TEST_EXCEPTION(Exception::Precondition, s1.start()); // cannot start twice
END_SECTION

START_SECTION((bool stop()))
  const double t_wait = 0.2;
  const double t_wait_more = 0.1;
  StopWatch s, s_nostop, s_reset, s_resume;
  s.start();
  s_nostop.start();
  s_reset.start();
  s_resume.resume();
  wait(t_wait);
  s.stop();
  s_resume.stop();
  TEST_EXCEPTION(Exception::Precondition, s.stop()); // cannot stop twice

  TEST_EQUAL(s.getClockTime() > 0.1, true)
  TEST_EQUAL(s.getClockTime() < 0.3, true)
  
  double t1 = s.getCPUTime();
  double t2 = s.getClockTime();
  double t3 = s.getSystemTime();
  double t4 = s.getUserTime();
  s_reset.reset();
  TEST_EQUAL(s_reset.isRunning(), true); // keeps on running
  s_resume.resume();
  // wait some more
  wait(t_wait_more);
  // ... and see if time is still the old one
  TEST_EQUAL(s.getCPUTime(), t1)
  TEST_EQUAL(s.getClockTime(), t2)
  TEST_EQUAL(s.getSystemTime(), t3)
  TEST_EQUAL(s.getUserTime(), t4)
  TEST_EQUAL(s.getCPUTime(), t1)

  TEST_EQUAL(s.getCPUTime() > t_wait / 2, true) // waiting costs CPU time in our implementation... just not sure how much...
  TEST_EQUAL(s.getClockTime() > t_wait * 0.95, true) // and must consume wall time
  TEST_EQUAL(s.getClockTime() < t_wait * 3, true) // be a bit more loose if e.g. a VM is busy
  TEST_EQUAL(s.getUserTime() > t_wait / 2, true) //  and some user time
  TEST_EQUAL(s.getUserTime() < t_wait * 2, true)
  TEST_EQUAL(s.getSystemTime() < t_wait, true) // and usually quite few system time
                                               //(not guaranteed on VMs, therefore do a trivial check)

  // the watch that never stopped should be ahead...
  TEST_EQUAL(s.getCPUTime() < s_nostop.getCPUTime(), true) 
  TEST_EQUAL(s.getClockTime() < s_nostop.getClockTime(), true)
  TEST_EQUAL(s.getUserTime() < s_nostop.getUserTime(), true)
  TEST_EQUAL(s.getSystemTime() <= s_nostop.getSystemTime(), true)

  s.reset(); // was stopped, so remains stopped
  TEST_EQUAL(s.isRunning(), false);
  TEST_EQUAL(s == StopWatch(), true);

  // kept on running the whole time after reset above .. should accumulate time
  TEST_EQUAL(s_reset.getCPUTime() > 0, true);

  // don't stop the timer.. just keep running and query on the fly
  TEST_EQUAL(s_resume.getCPUTime() > (t_wait_more + t_wait) / 2, true) // waiting costs CPU time in our implementation... just not sure how much...
  TEST_EQUAL(s_resume.getClockTime() > (t_wait_more + t_wait) * 0.95, true) //  must consume wall time
END_SECTION

START_SECTION((void clear()))
  StopWatch s;
  s.start();
  s.clear();
  TEST_EQUAL(s.isRunning(), false);
  TEST_EQUAL(s == StopWatch(), true);
  END_SECTION

START_SECTION((void reset()))
  NOT_TESTABLE; // done above to save Test time
END_SECTION

START_SECTION((void resume()))
  StopWatch s1;
  s1.start();
  TEST_EXCEPTION(Exception::Precondition, s1.resume()); // cannot start twice
END_SECTION

START_SECTION((double getCPUTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((double getClockTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((double getSystemTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((double getUserTime() const ))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((virtual ~StopWatch()))
  NOT_TESTABLE; // done above
END_SECTION

START_SECTION((static String toString(double time)))
  
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
