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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

	/**
	@brief Some functions to get system information

	Supports current memory and peak memory consumption.

	*/
	class OPENMS_DLLAPI SysInfo
	{
		public:
			/// Get memory consumption in KiloBytes (KB)
      /// On Windows, this is equivalent to 'Peak Working Set (Memory)' in Task Manager.
      /// On other OS this might be very unreliable, depending on operating system and kernel version.
      ///
			/// @param mem_virtual Total virtual memory currently allocated by this process
			/// @return True on success, false otherwise. If false is returned, then @p mem_virtual is set to 0.
			static bool getProcessMemoryConsumption(size_t& mem_virtual);
  
      /// Get peak memory consumption in KiloBytes (KB)
      /// On Windows, this is equivalent to 'Working Set (Memory)' in Task Manager.
      /// On other OS this might be very unreliable, depending on operating system and kernel version.
      ///
      /// @param mem_virtual Total virtual memory allocated by this process
      /// @return True on success, false otherwise. If false is returned, then @p mem_virtual is set to 0.
      static bool getProcessPeakMemoryConsumption(size_t& mem_virtual);

      /**
        @brief A convenience class to report either absolute or delta (between two timepoints) RAM usage

        Working RAM and Peak RAM usage are recorded at two time points ('before' and 'after').
        @note Peak RAM is only supported on WindowsOS; other OS will only report Working RAM usage
        
        When constructed, MemUsage automatically queries the present RAM usage (first timepoint), i.e. calls @ref before().
        Data for the second timepoint can be recorded using @ref after().

        @ref delta() reports the difference between the timepoints (before -> after);
        @ref usage() reports only the second timepoint's absolute value (after).

        When @ref delta() or @ref usage() are called, and the second timepoint is not recorded yet, this will be done internally.

      */
      struct OPENMS_DLLAPI MemUsage
      {
        size_t mem_before, mem_before_peak, mem_after, mem_after_peak;

        /// C'tor, calls @ref before() automatically
        MemUsage();

        /// forget all data (you need to call @ref before() again)
        void reset();
        /// record data for the first timepoint
        void before();
        /// record data for the second timepoint
        void after();
        /// get difference in memory usage between the two timepoints
        /// @ref after() will be called unless it was called earlier
        String delta(const String& event = "delta");

        /// get current memory usage (i.e. 'after')
        /// @ref after() will be called unless it was called earlier
        String usage();

      private:
        // convert difference to string
        String diff_str_(size_t mem_before, size_t mem_after);

      };
  };
}

