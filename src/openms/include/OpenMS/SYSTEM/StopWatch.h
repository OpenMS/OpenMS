// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#ifdef OPENMS_HAS_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef OPENMS_HAS_TIME_H
#include <ctime>
#endif

namespace OpenMS
{
  /**
      @brief This class is used to determine the current process' CPU (user and/or kernel) and wall time.

      CPU time is measured as sum of all threads of the current process.

      To read a time, the stopwatch must be started beforehand, but not necessarily stopped.
      
      You can stop() the timer and resume() it later to omit intermediate steps which should not count towards 
      the measured times.


      @ingroup System
  */
  class OPENMS_DLLAPI StopWatch
  {
public:
    /**	Starting, Stopping and Resetting the stop watch
    */
    //@{

    /**
        @brief Start the stop watch.
        
        If the watch holds data from previous measurements, these will be reset before starting up,
        i.e. it is not possible to resume by start(), stop(), start().
        Use resume(), stop(), resume() instead.
    
        @throw Exception::Precondition if the StopWatch is already running
    */
    void start();

    /** 
        @brief Stop the stop watch (can be resumed later).
        If the stop watch was not running an exception is thrown.

        @throw Exception::Precondition if the StopWatch is not running
    */
    void stop();

    /**
      @brief Resume a stopped StopWatch

      @throw Exception::Precondition if the StopWatch is already running
    */
    void resume();

    /** 
        @brief Clear the stop watch but keep running.
        
        The stop watch is reset to 0, but not stopped (if running).
        @see clear
    */
    void reset();

    /**	Clear and stop the stop watch.
        This sets the stop watch to zero and stops it when running.
        @see	reset
*/
    void clear();

    //@}

    /**	@name Readout of the StopWatch
    */
    //@{

    /**	Get clock time.
            Return the accumulated wall clock (real) time in seconds.
    */
    double getClockTime() const;

    /**	Get user time.
            Return the accumulated user time in seconds across all threads.
    */
    double getUserTime() const;

    /**	Get user time.
            Return the accumulated system time in seconds across all threads.
    */
    double getSystemTime() const;

    /**	Get CPU time.
            Return the accumulated CPU time in seconds.
            CPU time is the sum of user time and system time across all threads.
    */
    double getCPUTime() const;

    /**	@name	Predicates
    */
    //@{

    /**	Return true if the stop watch is running.
            @return bool <b>true</b> if the stop watch is running, <b>false</b> otherwise
    */
    bool isRunning() const;

    /**	Equality operator.
            Return <b>true</b> if two stop watches are equal, i.e. they contain exactly
            the same time intervals for clock, user and system time and have the
            same running status.
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> on equality, <b>false</b> otherwise
    */
    bool operator==(const StopWatch & stop_watch) const;

    /**	Inequality operator.
            Return <b>false</b> if two stop watches differ in any way, i.e. they differ
            in either the clock, user, or system time or have a different
            running status.
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> on inequality, <b>false</b> otherwise
    */
    bool operator!=(const StopWatch & stop_watch) const;

    /**	Lesser than operator.
            Return true, if the stop watch is in all timings lesser than the
            stop watch to be compared with (clock, user and system time).
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> if all times are lesser
    */
    bool operator<(const StopWatch & stop_watch) const;

    /**	Lesser or equal operator.
            Return true, if the stop watch is in all timings lesser or equal than the
            stop watch to be compared with (clock, user and system time).
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> if all times are lesser or equal
    */
    bool operator<=(const StopWatch & stop_watch) const;

    /**	Greater or equal operator.
            Return true, if the stop watch is in all timings greater or equal than the
            stop watch to be compared with (clock, user and system time).
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> if all times are greater or equal
    */
    bool operator>=(const StopWatch & stop_watch) const;

    /**	Greater operator.
            Return true, if the stop watch is in all timings greater than the
            stop watch to be compared with (clock, user and system time).
            @param stop_watch the stop watch to compare with
            @return bool <b>true</b> if all times are greater
    */
    bool operator>(const StopWatch & stop_watch) const;

    //@}

    /**
      @brief get a compact representation of the current time status.
      
      The output will be something like: 
      2.10 s (wall), 1.67 s (CPU), 0.12 s (system), 1.54 s (user)
      
    */
    String toString() const;

    /**
      custom string formatting of time, using only the minimal number of units required (e.g., does not report hours when seconds suffice).
    */
    static String toString(const double time_in_seconds);

private:
  #ifdef OPENMS_WINDOWSPLATFORM
    typedef UInt64 TimeType; ///< do not use clock_t on Windows, since its not big enough for larger time intervals
    static const long long SecondsTo100Nano_;  ///< 10 million; convert from 100 nanosecond ticks to seconds (factor of 1 billion/100 = 10 million)
  #else
    typedef clock_t TimeType;
    static const PointerSizeInt cpu_speed_; ///< POSIX API returns CPU ticks, so we need to divide by CPU speed
  #endif

    struct TimeDiff_
    {
      TimeType user_ticks{ 0 }; ///< platform dependent value (might be CPU ticks or time intervals)
      TimeType kernel_ticks{ 0 }; ///< platform dependent value (might be CPU ticks or time intervals)
      PointerSizeInt start_time{ 0 }; ///< time in seconds (relative or absolute, depending on usage)
      PointerSizeInt start_time_usec{ 0 }; ///< time in microseconds (relative or absolute, depending on usage)

      double userTime() const;
      double kernelTime() const;
      double getCPUTime() const;
      double clockTime() const;

      TimeDiff_ operator-(const TimeDiff_& earlier) const;
      TimeDiff_& operator+=(const TimeDiff_& other);
      bool operator==(const TimeDiff_& rhs) const;

      private:
        double ticksToSeconds_(TimeType in) const;
    };


    /// get the absolute times for current system, user and kernel times
    TimeDiff_ snapShot_() const;

    /// currently accumulated times between start to stop intervals (initially 0),
    /// not counting the currently running interval which started at last_start_
    TimeDiff_ accumulated_times_;

    /// point in time of last start()
    TimeDiff_ last_start_;

    /// state of stop watch, either true(on) or false(off)
    bool is_running_ = false;

  };

}



