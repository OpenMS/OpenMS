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

#ifndef OPENMS_SYSTEM_STOPWATCH_H
#define OPENMS_SYSTEM_STOPWATCH_H

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
      @brief StopWatch Class.

      This class is used to determine the current process time.

      @ingroup System
  */
  class OPENMS_DLLAPI StopWatch
  {
public:

    /**	@name	Constructors and Destructors
    */
    //@{

    /**	Default constructor.
            Create a new stop watch. The stop watch is stopped.
    */
    StopWatch();

    /**	Copy constructor.
            Create a new stop watch from an existing stop watch.
    */
    StopWatch(const StopWatch & stop_watch);

    /**	Destructor.
            Destructs a stop watch object.
    */
    virtual ~StopWatch();

    //@}
    /**	Starting, Stopping and Resetting the stop watch
    */
    //@{

    /**	Clear and stop the stop watch.
            This sets the stop watch to zero and stops it when running.
            @see	reset
    */
    void clear();

    /** Start the stop watch.
            The stop watch is started. If the stop watch is already running, <b>false</b>
            is returned.
            @return bool <b>false</b> if the stop watch was already running, <b>true</b> otherwise
    */
    bool start();

    /** Stop the stop watch.
            The stop watch is stopped. If the stop watch was not running, <b>false</b>
            is returned.
            @return bool <b>false</b> if the was not running, <b>true</b> otherwise
    */
    bool stop();

    /** Clear the stop watch without stopping.
            The stop watch is cleared, but not stopped (if running).
            @see clear
    */
    void reset();

    //@}

    /**	@name Readout of the StopWatch
    */
    //@{

    /**	Get clock time.
            Return the accumulated clock (real) time in seconds.
    */
    double getClockTime() const;

    /**	Get user time.
            Return the accumulated user time in seconds.
    */
    double getUserTime() const;

    /**	Get user time.
            Return the accumulated system time in seconds.
    */
    double getSystemTime() const;

    /**	Get CPU time.
            Return the accumulated CPU time in seconds.
            CPU time is the sum of user time and system time.
    */
    double getCPUTime() const;

    //@}

    /**	@name	Assignment
    */
    //@{

    /**	Assignment operator.
            Assigns a stop watch from another. The two stop watch will then run
            synchronously.
            @return StopWatch <tt>*this</tt>
    */
    StopWatch & operator=(const StopWatch & stop_watch);

    //@}

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
    static String toString(double time);

private:

    static PointerSizeInt cpu_speed_;

#ifdef OPENMS_WINDOWSPLATFORM
    static PointerSizeInt clock_speed_;
    typedef OPENMS_UINT64_TYPE TimeType; // do not use clock_t on Windows, since its not big enough for larger time intervals
#else
    typedef clock_t TimeType;
#endif

    // state of stop watch, either true(on) or false(off)
    bool is_running_;

    // clock seconds value when the stop watch was last started
    PointerSizeInt start_secs_;

    // clock useconds value when the stop watch was last started
    PointerSizeInt start_usecs_;

    // user time when the stop watch was last started
    TimeType start_user_time_;

    // system time when the stop watch was last started
    TimeType start_system_time_;

    // current accumulated clock seconds
    PointerSizeInt current_secs_;

    // current accumulated clock useconds
    PointerSizeInt current_usecs_;

    // current accumulated user time
    TimeType current_user_time_;

    // current accumulated user time
    TimeType current_system_time_;
  };

}



#endif // OPENMS_SYSTEM_STOPWATCH_H
