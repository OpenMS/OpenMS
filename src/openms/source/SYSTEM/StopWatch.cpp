// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/SYSTEM/StopWatch.h>

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h>
#endif
#ifdef OPENMS_HAS_TIME_H
#include <ctime>
#endif
#ifdef OPENMS_HAS_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef OPENMS_HAS_SYS_TIMES_H
#include <sys/times.h>
#endif
#ifdef OPENMS_HAS_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h>
#include <sys/timeb.h>
#endif


namespace OpenMS
{

  PointerSizeInt StopWatch::cpu_speed_ = 0L;

#ifdef OPENMS_WINDOWSPLATFORM
  PointerSizeInt StopWatch::clock_speed_ = 0L;
#endif

  StopWatch::StopWatch() :
    is_running_(false),
    start_secs_(0),
    start_usecs_(0),
    start_user_time_(0),
    start_system_time_(0),
    current_secs_(0),
    current_usecs_(0),
    current_user_time_(0),
    current_system_time_(0)
  {
#ifdef OPENMS_HAS_SYSCONF
    if (cpu_speed_ == 0L)
    {
      cpu_speed_ = sysconf(_SC_CLK_TCK);
    }
#endif

#ifdef OPENMS_WINDOWSPLATFORM
    if (cpu_speed_ == 0L)
    {
      LARGE_INTEGER ticks;
      if (QueryPerformanceFrequency(&ticks))
      {
        cpu_speed_ = (PointerSizeInt) ticks.QuadPart;
      }
      else
      {
        cpu_speed_ = 1L;
      }
      clock_speed_ = CLOCKS_PER_SEC;
    }
#endif
  }

  StopWatch::StopWatch(const StopWatch & stop_watch) :
    is_running_(stop_watch.is_running_),
    start_secs_(stop_watch.start_secs_),
    start_usecs_(stop_watch.start_usecs_),
    start_user_time_(stop_watch.start_user_time_),
    start_system_time_(stop_watch.start_system_time_),
    current_secs_(stop_watch.current_secs_),
    current_usecs_(stop_watch.current_usecs_),
    current_user_time_(stop_watch.current_user_time_),
    current_system_time_(stop_watch.current_system_time_)
  {
  }

  StopWatch::~StopWatch()
  {
  }

  void StopWatch::clear()
  {
    is_running_ = false;
    current_secs_ = 0L;
    current_usecs_ = 0L;
    current_user_time_ = 0L;
    current_system_time_ = (TimeType)0;
  }

  bool StopWatch::start()
  {
    if (is_running_ == true)
    {
      /* tried to start a running stop_watch */
      return false;
    }

#ifdef OPENMS_WINDOWSPLATFORM
    LARGE_INTEGER tms;
    FILETIME kt, ut, ct, et;

    QueryPerformanceCounter(&tms);
    HANDLE my_id = GetCurrentProcess();
    GetProcessTimes(my_id, &ct, &et, &kt, &ut);
    ULARGE_INTEGER kernel_time;
    kernel_time.HighPart = kt.dwHighDateTime;
    kernel_time.LowPart = kt.dwLowDateTime;
    ULARGE_INTEGER user_time;
    user_time.HighPart = ut.dwHighDateTime;
    user_time.LowPart = ut.dwLowDateTime;

    start_secs_  = tms.QuadPart / cpu_speed_;
    start_usecs_ = (PointerSizeInt)((double)(tms.QuadPart - (start_secs_ * cpu_speed_)) / (double)(cpu_speed_) * 1000000.0);

    start_user_time_ = (TimeType) (user_time.QuadPart / 10);
    start_system_time_ = (TimeType) (kernel_time.QuadPart / 10);

#else

    struct tms tms_buffer;
    struct timeval timeval_buffer;
    struct timezone timezone_buffer;

    gettimeofday(&timeval_buffer, &timezone_buffer);
    times(&tms_buffer);

    start_secs_ = timeval_buffer.tv_sec;
    start_usecs_ = timeval_buffer.tv_usec;
    start_user_time_ = tms_buffer.tms_utime;
    start_system_time_ = tms_buffer.tms_stime;
#endif

    is_running_ = true;

    return true;
  }

  bool StopWatch::stop()
  {
    if (is_running_ == false) /* tried to stop a stopped stop_watch */
    {
      return false;
    }
#ifdef OPENMS_WINDOWSPLATFORM
    LARGE_INTEGER tms;

    QueryPerformanceCounter(&tms);
    FILETIME kt, ut, ct, et;

    HANDLE my_id = GetCurrentProcess();
    GetProcessTimes(my_id, &ct, &et, &kt, &ut);

    ULARGE_INTEGER kernel_time;
    kernel_time.HighPart = kt.dwHighDateTime;
    kernel_time.LowPart = kt.dwLowDateTime;
    ULARGE_INTEGER user_time;
    user_time.HighPart = ut.dwHighDateTime;
    user_time.LowPart = ut.dwLowDateTime;

    PointerSizeInt secs_to_add = tms.QuadPart / cpu_speed_;
    current_secs_ += secs_to_add - start_secs_;
    PointerSizeInt usecs_to_add = (PointerSizeInt)((double)(tms.QuadPart - secs_to_add * cpu_speed_) / (double)(cpu_speed_) * 1000000.0);
    current_usecs_ += usecs_to_add - start_usecs_;

    current_user_time_ += (TimeType) (user_time.QuadPart / 10 - start_user_time_);
    current_system_time_ += (TimeType) (kernel_time.QuadPart / 10 - start_system_time_);
#else
    struct tms tms_buffer;
    struct timeval timeval_buffer;
    struct timezone timezone_buffer;

    gettimeofday(&timeval_buffer, &timezone_buffer);
    times(&tms_buffer);

    current_secs_ += timeval_buffer.tv_sec - start_secs_;
    current_usecs_ += timeval_buffer.tv_usec - start_usecs_;

    current_user_time_ += tms_buffer.tms_utime - start_user_time_;
    current_system_time_ += tms_buffer.tms_stime - start_system_time_;
#endif

    is_running_ = false;

    return true;
  }

  void StopWatch::reset()
  {
    if (is_running_ == false)
    {
      clear();
    }
    else
    {
      stop();
      clear();
      start();
    }
  }

  //getClockTime returns the current amount of real (clock) time
  //accumulated by this stop_watch.  If the stop_watch is stopped, this is just
  //the total accumulated time.  If the stop_watch is running, this is the
  //accumulated time + the time since the stop_watch was last started.
  double StopWatch::getClockTime() const
  {
    PointerSizeInt elapsed_seconds = 0;
    PointerSizeInt elapsed_useconds = 0;

    if (is_running_ == false)
    {
      /* stop_watch is currently off, so just return accumulated time */
      elapsed_seconds = current_secs_;
      elapsed_useconds = current_usecs_;
    }
    else
    {
      /* stop_watch is currently running, so add the elapsed time since */
      /* the stop_watch was last started to the accumulated time        */
#ifdef OPENMS_WINDOWSPLATFORM
      LARGE_INTEGER tms;
      if (QueryPerformanceCounter(&tms))
      {
        PointerSizeInt secs_to_add = tms.QuadPart / cpu_speed_;
        elapsed_seconds = current_secs_ + secs_to_add - start_secs_;
        PointerSizeInt usecs_to_add = (PointerSizeInt)((double)(tms.QuadPart - secs_to_add * cpu_speed_) / (double)(cpu_speed_) * 1000000.0);
        elapsed_useconds  = current_usecs_ + usecs_to_add - start_usecs_;
      }
#else
      struct timeval timeval_buffer;
      struct timezone timezone_buffer;

      gettimeofday(&timeval_buffer, &timezone_buffer);

      elapsed_seconds = current_secs_ + timeval_buffer.tv_sec - start_secs_;
      elapsed_useconds = current_usecs_ + timeval_buffer.tv_usec - start_usecs_;
#endif
    }


    /* Adjust for the fact that the useconds may be negative. */
    /* If they are, take away 1 second and add 1 million      */
    /* microseconds until they are positive.                  */
    while (elapsed_useconds < 0L)
    {
      elapsed_useconds += 1000000L;
      --elapsed_seconds;
    }

    /* convert into floating point number of seconds */
    return (double)((double)elapsed_seconds + (double)elapsed_useconds / 1000000.0);
  }

  //getUserTime reports the current amount of user cpu time
  //accumulated by this StopWatch.  If the stop_watch is currently off,
  //this is just the accumulated time.  If the StopWatch is running, this
  //is the accumulated time plus the time since the stop_watch was last started.
  double StopWatch::getUserTime() const
  {
    double temp_value;

#ifdef OPENMS_WINDOWSPLATFORM
    FILETIME kt, ut, ct, et;
#else
    struct tms tms_buffer;
#endif
    if (is_running_ == false)
    {
      /* stop_watch is off, just return accumulated time */
      temp_value = (double)current_user_time_;
    }
    else
    {
      /* stop_watch is on, add current running time to accumulated time */
#ifdef OPENMS_WINDOWSPLATFORM
      HANDLE my_id = GetCurrentProcess();
      GetProcessTimes(my_id, &ct, &et, &kt, &ut);

      ULARGE_INTEGER kernel_time;
      kernel_time.HighPart = kt.dwHighDateTime;
      kernel_time.LowPart = kt.dwLowDateTime;
      ULARGE_INTEGER user_time;
      user_time.HighPart = ut.dwHighDateTime;
      user_time.LowPart = ut.dwLowDateTime;

      temp_value = (double)(current_user_time_ - start_user_time_ + user_time.QuadPart / 10.0);
#else
      times(&tms_buffer);
      temp_value = (double)(current_user_time_ - start_user_time_ + tms_buffer.tms_utime);
#endif
    }

#ifdef OPENMS_WINDOWSPLATFORM
    return (double)(temp_value / 1000000.0);

#else
    /* convert from clock ticks to seconds using the */
    /* cpu-speed value obtained in the constructor   */
    return (double)(temp_value / (double)cpu_speed_);

#endif
  }

  // system_time reports the current amount of system cpu time
  // accumulated by this StopWatch.  If the stop_watch is currently off,
  // this is just the accumulated time.  If the StopWatch is running, this
  // is the accumulated time plus the time since the stop_watch was last started
  double StopWatch::getSystemTime() const
  {
    double temp_value(0.0);

    if (is_running_ == false)
    {
      /* stop_watch is off, just return accumulated time */
      temp_value = (double)current_system_time_;
    }
    else
    {
      /* stop_watch is on, return accumulated plus current */
#ifdef OPENMS_WINDOWSPLATFORM
      //struct tms tms_buffer;
      FILETIME kt, ut, ct, et;
      //times(&tms_buffer);
      HANDLE my_id = GetCurrentProcess();
      GetProcessTimes(my_id, &ct, &et, &kt, &ut);

      ULARGE_INTEGER kernel_time;
      kernel_time.HighPart = kt.dwHighDateTime;
      kernel_time.LowPart = kt.dwLowDateTime;
      ULARGE_INTEGER user_time;
      user_time.HighPart = ut.dwHighDateTime;
      user_time.LowPart = ut.dwLowDateTime;
      temp_value = (double)((double)(current_system_time_ - start_system_time_) + kernel_time.QuadPart / 10.0);
#endif
    }

    /* convert from clock ticks to seconds using the */
    /* cpu-speed value obtained by the constructor   */
    return (double)(temp_value / 1000000.0);
  }

  StopWatch & StopWatch::operator=(const StopWatch & stop_watch)
  {
    if (this == &stop_watch)
    {
      return *this;
    }

    is_running_ = stop_watch.is_running_;
    start_secs_ = stop_watch.start_secs_;
    start_usecs_ = stop_watch.start_usecs_;
    start_user_time_ = stop_watch.start_user_time_;
    start_system_time_ = stop_watch.start_system_time_;
    current_secs_ = stop_watch.current_secs_;
    current_usecs_ = stop_watch.current_usecs_;
    current_user_time_ = stop_watch.current_user_time_;
    current_system_time_ = stop_watch.current_system_time_;

    return *this;
  }

  bool StopWatch::operator==(const StopWatch & stop_watch) const
  {
    return start_secs_ == stop_watch.start_secs_
           && start_usecs_ == stop_watch.start_usecs_
           && start_user_time_ == stop_watch.start_user_time_
           && start_system_time_ == stop_watch.start_system_time_
           && current_secs_ == stop_watch.current_secs_
           && current_usecs_ == stop_watch.current_usecs_
           && current_user_time_ == stop_watch.current_user_time_
           && current_system_time_ == stop_watch.current_system_time_;
  }


  String StopWatch::toString(double time)
  {
    int d(0), h(0), m(0);
    double s(0), s_single(0);

    TimeType time_i = (TimeType) time; // trunc to integer

    // compute days
    d = int(time_i / (3600*24));
    time_i -= d*(3600*24);
    time -= d*(3600*24);

    // hours
    h = int(time_i / 3600);
    time_i -= h*3600;
    time -= h*3600;

    // minutes
    m = int(time_i / 60);
    time_i -= m*60;
    time -= m*60;

    s_single = time;
    s = (double) time_i;


    String s_d = String(d);
    String s_h = String(h).fillLeft('0', 2) + ":";
    String s_m = String(m).fillLeft('0', 2) + ":";
    String s_s = String(s).fillLeft('0', 2); // if we show seconds in combination with minutes, we round to nominal 

    String s_s_single = String::number(s_single, 2); // second (shown by itself with no minutes) has two digits after decimal

    return ( (d>0 ? s_d + "d " + s_h + s_m + s_s + " h" :
             (h>0 ?              s_h + s_m + s_s + " h" :
             (m>0 ?                    s_m + s_s + " m" :
             (                        s_s_single + " s")))));

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

} // namespace OpenMS
