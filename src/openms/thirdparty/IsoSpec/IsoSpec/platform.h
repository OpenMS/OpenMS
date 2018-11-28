
#pragma once

#if !defined(ISOSPEC_BUILDING_R)
#define ISOSPEC_BUILDING_R false
#endif

#if !defined(ISOSPEC_BUILDING_CPP)
#define ISOSPEC_BUILDING_CPP true
#endif

#if !defined(ISOSPEC_BUILDING_PYTHON)
#define ISOSPEC_BUILDING_PYTHON false
#endif

#if !defined(ISOSPEC_BUILDING_OPENMS)
#define ISOSPEC_BUILDING_OPENMS false
#endif

#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY true
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS false /* CYGWIN doesn't really count as Windows for our purposes, we'll be using UNIX API anyway */
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN true
#define ISOSPEC_TEST_GOT_MMAN true
#elif defined(__MINGW32__) || defined(_WIN32)
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY false
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS true
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN false
#define ISOSPEC_TEST_GOT_MMAN true
#else
#define ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY false /* Well, probably... */
#define ISOSPEC_TEST_WE_ARE_ON_WINDOWS false
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN false
#define ISOSPEC_TEST_GOT_MMAN false
#endif

#if !defined(ISOSPEC_USE_PTHREADS)
#define ISOSPEC_USE_PTHREADS false /* TODO: possibly put a macro here to detect whether we */
#endif                             /* can/should use pthreads - or rip them out altogether.
                                    * Investigate whether the performance advantage of pthreads on
                                    * some platforms (*cough* CYGWIN *cough*) is still large
                                    * enough to justify keeping both implementations around */

#if !defined(ISOSPEC_WE_ARE_ON_UNIX_YAY)
#define ISOSPEC_WE_ARE_ON_UNIX_YAY ISOSPEC_TEST_WE_ARE_ON_UNIX_YAY
#endif

#if !defined(ISOSPEC_WE_ARE_ON_WINDOWS)
#define ISOSPEC_WE_ARE_ON_WINDOWS ISOSPEC_TEST_WE_ARE_ON_WINDOWS
#endif

#if !defined(ISOSPEC_GOT_SYSTEM_MMAN)
#define ISOSPEC_GOT_SYSTEM_MMAN ISOSPEC_TEST_GOT_SYSTEM_MMAN
#endif

#if !defined(ISOSPEC_GOT_MMAN)
#define ISOSPEC_GOT_MMAN ISOSPEC_TEST_GOT_MMAN
#endif


// Note: __GNUC__ is defined by clang and gcc
#ifdef __GNUC__
#define ISOSPEC_LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define ISOSPEC_UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)
// For aggressive inlining 
#define ISOSPEC_FORCE_INLINE __attribute__ ((always_inline)) inline
#elif defined _MSC_VER
#define ISOSPEC_LIKELY(condition) condition
#define ISOSPEC_UNLIKELY(condition) condition
#define ISOSPEC_FORCE_INLINE __forceinline inline
#else
#define ISOSPEC_LIKELY(condition) condition
#define ISOSPEC_UNLIKELY(condition) condition
#define ISOSPEC_FORCE_INLINE inline
#endif


#if ISOSPEC_GOT_MMAN
    #if ISOSPEC_GOT_SYSTEM_MMAN
        #include <sys/mman.h>
    #else
        #include "mman.h"
    #endif
#else
    #include <stdlib.h>     /* malloc, free, rand */
#endif


#if defined(OPENMS_DLLAPI) /* IsoSpec is being built as a part of OpenMS: use their visibility macros */
#define ISOSPEC_EXPORT_SYMBOL OPENMS_DLLAPI
#else /* it's a can of worms we don't yet want to open ourselves though... */
#define ISOSPEC_EXPORT_SYMBOL
#endif
