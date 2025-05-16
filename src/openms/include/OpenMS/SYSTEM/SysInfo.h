// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <array>

#if defined(_M_X64) || defined(__x86_64__) || defined(__amd64__)
    #define OPENMS_ARCH_X64 1
#endif

namespace OpenMS
{

  /// convert bytes to a human readable unit (TiB, GiB, MiB, KiB), e.g. "45.34 MiB"
  OPENMS_DLLAPI std::string bytesToHumanReadable(UInt64 bytes);

	/**
	@brief Some functions to get system information

	Supports current memory and peak memory consumption and CPU capabilities (for x64)

	*/
	class OPENMS_DLLAPI SysInfo
	{
		public:
    
      /// A hierachical list of instruction sets supported by x86_64 CPUs.
      /// If an instruction is supported, all lower instructions are supported too.
      enum class X64_SIMDLevel : int
      {
          SSE2           = 2,  // 2:  SSE2
          SSE3           = 3,  // 3:  SSE3
          SSSE3          = 4,  // 4:  SSSE3
          SSE41          = 5,  // 5:  SSE4.1
          SSE42          = 6,  // 6:  SSE4.2
          AVX            = 7,  // 7:  AVX
          AVX2           = 8,  // 8:  AVX2
          AVX512F        = 9,  // 9:  AVX512F
          AVX512_BW_DQ_VL= 10  // 10: AVX512BW/DQ/VL
      };


      static inline const std::array<std::string, 11> X64_SIMDLevelNames = {
        "NA",
        "NA",
        "SSE2",          
        "SSE3",          
        "SSSE3",         
        "SSE4.1",        
        "SSE4.2",        
        "AVX",           
        "AVX2",          
        "AVX512F",       
        "AVX512BW/DQ/VL" 
      };
    
      /// @brief Get the maximum CPU capabilities of the current CPU this code runs on (only valid on X64 currently; reports '?' for ARM and AARCH until implemented)
      static String maxSIMDCapabilityAsString();

      // only valid on X64 ; to not use on ARM
      #ifdef OPENMS_ARCH_X64
        /// report the maximum CPU capabilities of the current CPU this code runs on
        static X64_SIMDLevel getHighestRuntimeCPUFeature();
      #endif

      /// @brief Check if the CPU supports the given instruction set and abort the program if not (to avoid runtime errors)
      static void fatalCPUCapabilityCheck();
    
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

