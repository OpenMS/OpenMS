// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/CONCEPT/Exception.h>

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h>
#endif
#ifdef OPENMS_HAS_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h>
#include <sys/timeb.h>
#endif


namespace OpenMS
{
#ifdef OPENMS_WINDOWSPLATFORM
  const long long StopWatch::SecondsTo100Nano_ = 10000000LL;
#else
  const PointerSizeInt StopWatch::cpu_speed_ = sysconf(_SC_CLK_TCK);
#endif

  void StopWatch::clear()
  { // stopped when running
    *this = StopWatch(); // default init
  }

  void StopWatch::start()
  {
    if (is_running_)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StopWatch is already started!");
    }

    clear();
    last_start_ = snapShot_();
    is_running_ = true;
  }

  void StopWatch::stop()
  {
    if (!is_running_)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StopWatch cannot be stopped if not running!");
    }
    
    TimeDiff_ now = snapShot_();
    auto diff = now - last_start_;
    accumulated_times_ += diff;

    is_running_ = false;
  }

  void StopWatch::resume()
  {
    if (is_running_)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StopWatch cannot be resumed if already running!");
    }

    last_start_ = snapShot_();
    is_running_ = true;
  }

  void StopWatch::reset()
  {
    if (is_running_ == false)
    {
      clear();
    }
    else
    {
      clear();
      start();
    }
  }

  StopWatch::TimeDiff_ StopWatch::snapShot_() const
  {
    TimeDiff_ t;

#ifdef OPENMS_WINDOWSPLATFORM
    LARGE_INTEGER lpFrequency; ///< counts of QueryPerformanceCounter per second; fixed at boot time;
    QueryPerformanceFrequency(&lpFrequency);

    LARGE_INTEGER tms;
    //QueryPerformanceCounter returns values that represent time in units of 1 / (the frequency of the performance counter obtained from QueryPerformanceFrequency)
    QueryPerformanceCounter(&tms);
    t.start_time = tms.QuadPart / lpFrequency.QuadPart;
    const double secToUsec = 1e6;
    t.start_time_usec = (PointerSizeInt)((double)(tms.QuadPart - (t.start_time * lpFrequency.QuadPart)) / (double)(lpFrequency.QuadPart) * secToUsec);

    FILETIME ct, et, kt, ut; 
    // ct is creation time of process, but et is end-time (which is undefined for running processes like ours);
    // Thus so cannot be used to measure wall time and we need QueryPerformanceCounter from above
    GetProcessTimes(GetCurrentProcess(), &ct, &et, &kt, &ut);
    ULARGE_INTEGER kernel_time;
    kernel_time.HighPart = kt.dwHighDateTime;
    kernel_time.LowPart = kt.dwLowDateTime;
    ULARGE_INTEGER user_time;
    user_time.HighPart = ut.dwHighDateTime;
    user_time.LowPart = ut.dwLowDateTime;
    
    t.user_ticks = (TimeType)user_time.QuadPart;
    t.kernel_ticks = (TimeType)kernel_time.QuadPart;
#else

    struct timeval timeval_buffer;
    struct timezone timezone_buffer;
    gettimeofday(&timeval_buffer, &timezone_buffer);
    t.start_time = timeval_buffer.tv_sec; // seconds since 1970-01-01 00:00
    t.start_time_usec = timeval_buffer.tv_usec; // additional(!) usec 

    struct tms tms_buffer;
    times(&tms_buffer);
    t.user_ticks = tms_buffer.tms_utime; // reports value in CPU clock ticks
    t.kernel_ticks = tms_buffer.tms_stime; // reports value in CPU clock ticks
#endif

    return t;
  }

  //getClockTime returns the current amount of real (clock) time
  //accumulated by this stop_watch.  If the stop_watch is stopped, this is just
  //the total accumulated time.  If the stop_watch is running, this is the
  //accumulated time + the time since the stop_watch was last started.
  double StopWatch::getClockTime() const
  {
    if (is_running_ == false)
    {
      /* stop_watch is currently off, so just return accumulated time */
      return accumulated_times_.clockTime();
    }

    /* stop_watch is currently running, so add the elapsed time since */
    /* the stop_watch was last started to the accumulated time        */
    auto now = snapShot_();
    auto diff = now - last_start_;
       
    /* convert into floating point number of seconds */
    return accumulated_times_.clockTime() + diff.clockTime();
  }

  //getUserTime reports the current amount of user cpu time
  //accumulated by this StopWatch.  If the stop_watch is currently off,
  //this is just the accumulated time.  If the StopWatch is running, this
  //is the accumulated time plus the time since the stop_watch was last started.
  double StopWatch::getUserTime() const
  {
    if (is_running_ == false)
    {
      /* stop_watch is currently off, so just return accumulated time */
      return accumulated_times_.userTime();
    }

    /* stop_watch is currently running, so add the elapsed time since */
    /* the stop_watch was last started to the accumulated time        */
    auto now = snapShot_();
    auto diff = now - last_start_;

    /* convert into floating point number of seconds */
    return accumulated_times_.userTime() + diff.userTime();
  }

  // system_time reports the current amount of system cpu time
  // accumulated by this StopWatch.  If the stop_watch is currently off,
  // this is just the accumulated time.  If the StopWatch is running, this
  // is the accumulated time plus the time since the stop_watch was last started
  double StopWatch::getSystemTime() const
  {
    if (is_running_ == false)
    {
      /* stop_watch is currently off, so just return accumulated time */
      return accumulated_times_.kernelTime();
    }

    /* stop_watch is currently running, so add the elapsed time since */
    /* the stop_watch was last started to the accumulated time        */
    auto now = snapShot_();
    auto diff = now - last_start_;

    /* convert into floating point number of seconds */
    return accumulated_times_.kernelTime() + diff.kernelTime();
  }

  bool StopWatch::operator==(const StopWatch& rhs) const
  {
    return accumulated_times_ == rhs.accumulated_times_
           && last_start_ == rhs.last_start_
           && is_running_ == rhs.is_running_;
  }


  String StopWatch::toString(const double time_in_seconds)
  {
    int d(0), h(0), m(0), s(0);

    TimeType time_i = (TimeType)time_in_seconds; // trunc to integer

    // compute days
    d = int(time_i / (3600*24));
    time_i -= d*(3600*24);

    // hours
    h = int(time_i / 3600);
    time_i -= h*3600;

    // minutes
    m = int(time_i / 60);
    time_i -= m*60;

    s = int(time_i);


    String s_d = String(d);
    String s_h = String(h).fillLeft('0', 2) + ":";
    String s_m = String(m).fillLeft('0', 2) + ":";
    String s_s = String(s).fillLeft('0', 2); // if we show seconds in combination with minutes, we round to nominal 

    return ( (d>0 ? s_d + "d " + s_h + s_m + s_s + " h" :
             (h>0 ?              s_h + s_m + s_s + " h" :
             (m>0 ?                    s_m + s_s + " m" :
             (      String::number(time_in_seconds, 2) + " s"))))); // second (shown by itself with no minutes) has two digits after decimal

  }

  String StopWatch::toString() const
  {
    return(
      StopWatch::toString(this->getClockTime()) + " (wall), " +
      StopWatch::toString(this->getCPUTime()) + " (CPU), " + 
      StopWatch::toString(this->getSystemTime()) + " (system), " +
      StopWatch::toString(this->getUserTime()) + " (user)"
    );
  }

  double StopWatch::getCPUTime() const
  {
    return getUserTime() + getSystemTime();
  }

  bool StopWatch::isRunning() const
  {
    return is_running_;
  }

  bool StopWatch::operator!=(const StopWatch & stop_watch) const
  {
    return !(*this == stop_watch);
  }

  bool StopWatch::operator<(const StopWatch & stop_watch) const
  {
    return getCPUTime() < stop_watch.getCPUTime();
  }

 
  bool StopWatch::operator<=(const StopWatch & stop_watch) const
  {
    return !(stop_watch < *this);
  }

  bool StopWatch::operator>=(const StopWatch & stop_watch) const
  {
    return !(*this < stop_watch);
  }

  bool StopWatch::operator>(const StopWatch & stop_watch) const
  {
    return stop_watch < *this;
  }

  inline double StopWatch::TimeDiff_::userTime() const
  {
    return ticksToSeconds_(user_ticks);
  }

  inline double StopWatch::TimeDiff_::kernelTime() const
  {
    return ticksToSeconds_(kernel_ticks);
  }

  inline double StopWatch::TimeDiff_::getCPUTime() const
  {
    return userTime() + kernelTime();
  }

  inline double StopWatch::TimeDiff_::clockTime() const
  {
    return (double)start_time + (double)start_time_usec / 1e6;
  }

  inline double StopWatch::TimeDiff_::ticksToSeconds_(TimeType in) const
  {
#ifdef OPENMS_WINDOWSPLATFORM
    return in / double(StopWatch::SecondsTo100Nano_);
#else
    return in / double(StopWatch::cpu_speed_); // technically, this is inaccurate since CPU speed may not be constant (turbo-boost)... but finding a better solution is hard...
#endif
  }

  StopWatch::TimeDiff_ StopWatch::TimeDiff_::operator-(const StopWatch::TimeDiff_& earlier) const
  {
    TimeDiff_ diff(*this);
    diff.kernel_ticks -= earlier.kernel_ticks;
    diff.user_ticks -= earlier.user_ticks;
    diff.start_time -= earlier.start_time;
    diff.start_time_usec -= earlier.start_time_usec;

    /* Adjust for the fact that the usec may be negative.     */
    /* If they are, take away 1 second and add 1 million      */
    /* microseconds until they are positive.                  */
    while (diff.start_time_usec < 0L)
    {
      --diff.start_time;
      diff.start_time_usec += 1000000L;
    }
    return diff;
  }

  StopWatch::TimeDiff_& StopWatch::TimeDiff_::operator+=(const StopWatch::TimeDiff_& other)
  {
    user_ticks += other.user_ticks;
    kernel_ticks += other.kernel_ticks;
    start_time += other.start_time;
    start_time_usec += other.start_time_usec;

    while (start_time_usec > 1000000L)
    {
      ++start_time;
      start_time_usec -= 1000000L;
    }

    return *this;
  }

  bool StopWatch::TimeDiff_::operator==(const TimeDiff_& rhs) const
  {
    return user_ticks == rhs.user_ticks &&
      kernel_ticks == rhs.kernel_ticks &&
      start_time == rhs.start_time &&
      start_time_usec == rhs.start_time_usec;
  }

} // namespace OpenMS
