// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <vector>

namespace OpenMS
{

  /// convert bytes to a human readable unit (TiB, GiB, MiB, KiB), e.g. "45.34 MiB"
  OPENMS_DLLAPI std::string bytesToHumanReadable(UInt64 bytes);

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
			/// @param[out] mem_virtual Total virtual memory currently allocated by this process
			/// @return True on success, false otherwise. If false is returned, then @p mem_virtual is set to 0.
			static bool getProcessMemoryConsumption(size_t& mem_virtual);
  
      /// Get peak memory consumption in KiloBytes (KB)
      /// On Windows, this is equivalent to 'Working Set (Memory)' in Task Manager.
      /// On other OS this might be very unreliable, depending on operating system and kernel version.
      ///
      /// @param[out] mem_virtual Total virtual memory allocated by this process
      /// @return True on success, false otherwise. If false is returned, then @p mem_virtual is set to 0.
      static bool getProcessPeakMemoryConsumption(size_t& mem_virtual);

      /**
        @brief Get currently available physical system memory in KiloBytes (KB).

        On Linux this prefers MemAvailable from /proc/meminfo so file-system
        cache that can be reclaimed by the kernel is counted as available.

        @param[out] mem_available Available physical system memory
        @return True on success, false otherwise. If false is returned, then @p mem_available is set to 0.
      */
      static bool getFreeSystemMemory(size_t& mem_available);

      /// `@brief` Get the current process ID
      ///
      /// `@return` The process ID of the current process
      static Int64 getProcessId();

      /**
        @brief The number of physical cores the process may run on (at least 1).

        Counts cores, not their hardware threads (SMT), which std::thread::hardware_concurrency() and
        omp_get_num_procs() count. On Linux: the cores of the CPUs in the affinity mask of the calling thread (taskset,
        cpusets), or of @p cpus if given (e.g. the CPUs of the OpenMP places, if the OpenMP runtime has bound the
        calling thread already), at most the CPU quota of the cgroup of the process and of those above it as far as
        visible (cgroup v2 cpu.max, rounded up). On macOS: hw.physicalcpu. On Windows: the cores of the processors in
        the affinity mask of the process, in the processor group of the calling thread. Falls back to the hardware
        threads where this cannot be found out.

        @param cpus Linux only: the logical CPUs (numbers as in /sys/devices/system/cpu) to count the cores of instead
                    of those of the affinity mask
      */
      static UInt getAvailablePhysicalCores(const std::vector<int>& cpus = {});
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
        std::string delta(const std::string& event = "delta");

        /// get current memory usage (i.e. 'after')
        /// @ref after() will be called unless it was called earlier
        std::string usage();

      private:
        // convert difference to string
        std::string diff_str_(size_t mem_before, size_t mem_after);

      };
  };
}
