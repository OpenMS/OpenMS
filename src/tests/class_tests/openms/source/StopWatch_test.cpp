// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
  TEST_FALSE(s1 == s2) // before stop
  s1.stop();
  TEST_FALSE(s1 == s2)
  s2 = s1;
  TEST_TRUE(s1 == s2)

  StopWatch s3(s1);
  TEST_TRUE(s1 == s3)
  
  StopWatch s4;
  s1.reset();
  TEST_TRUE(s1 == s4)

  s1.start();
  s2.start();
  wait(0.01);
  s1.stop();

  wait(0.01);
  s2.stop();

  TEST_FALSE(s1 == s2)
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
  std::cout << "Usertime: " << s.getUserTime() << "\n";
#ifdef OPENMS_WINDOWSPLATFORM
  // workaround for Windows-CI on VMs which report usertime = 0 ...
  TEST_EQUAL(s.getUserTime() >= 0, true)//  and some user time
#else
  TEST_EQUAL(s.getUserTime() > t_wait / 2, true)//  and some user time
#endif
  TEST_EQUAL(s.getUserTime() < t_wait * 2, true)
  std::cout << "Systemtime: " << s.getSystemTime() << "\n";
  TEST_EQUAL(s.getSystemTime() < t_wait * 2, true)// and usually quite few system time
                                                  // (not guaranteed on VMs, therefore do a trivial check)

  // the watch that never stopped should be ahead...
  TEST_EQUAL(s.getCPUTime() < s_nostop.getCPUTime(), true) 
  TEST_EQUAL(s.getClockTime() < s_nostop.getClockTime(), true)
  std::cout << "compare: " << s.getUserTime() << " <> " << s_nostop.getUserTime() << "\n";
#ifdef OPENMS_WINDOWSPLATFORM
  // workaround for Windows-CI on VMs which report usertime = 0 ...
  TEST_EQUAL(s.getUserTime() <= s_nostop.getUserTime(), true)
#else
  TEST_EQUAL(s.getUserTime() < s_nostop.getUserTime(), true)
#endif
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
