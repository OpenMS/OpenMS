// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// OpenMSConfig.h — Shared visibility/export macros for all OpenMS sub-libraries
//
// This header replaces the auto-generated export header from generate_export_header().
// It defines OPENMS_DLLAPI for use across all three sub-libraries (Core, IO, Algo)
// while maintaining backward compatibility with existing #include directives.

#ifndef OPENMS_EXPORT_H
#define OPENMS_EXPORT_H

#ifdef OPENMS_STATIC_DEFINE
  #define OPENMS_DLLAPI
  #define OPENMS_NO_EXPORT
#else
  #ifdef _MSC_VER
    // Windows: export when building any OpenMS sub-library, import when consuming
    #if defined(OpenMS_Core_EXPORTS) || defined(OpenMS_IO_EXPORTS) || defined(OpenMS_Algo_EXPORTS) || defined(OpenMS_EXPORTS)
      #define OPENMS_DLLAPI __declspec(dllexport)
    #else
      #define OPENMS_DLLAPI __declspec(dllimport)
    #endif
    #define OPENMS_NO_EXPORT
    // Suppress DLL interface warnings for STL types
    #pragma warning(disable: 4251)
  #else
    // Unix: use visibility attribute
    #define OPENMS_DLLAPI __attribute__((visibility("default")))
    #define OPENMS_NO_EXPORT __attribute__((visibility("hidden")))
  #endif
#endif

#ifndef OPENMS_DEPRECATED
  #define OPENMS_DEPRECATED __attribute__((__deprecated__))
#endif

#ifndef OPENMS_DEPRECATED_EXPORT
  #define OPENMS_DEPRECATED_EXPORT OPENMS_DLLAPI OPENMS_DEPRECATED
#endif

#ifndef OPENMS_DEPRECATED_NO_EXPORT
  #define OPENMS_DEPRECATED_NO_EXPORT OPENMS_NO_EXPORT OPENMS_DEPRECATED
#endif

// NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if)
#if 0 /* DEFINE_NO_DEPRECATED */
  #ifndef OPENMS_NO_DEPRECATED
    #define OPENMS_NO_DEPRECATED
  #endif
#endif

#endif // OPENMS_EXPORT_H
