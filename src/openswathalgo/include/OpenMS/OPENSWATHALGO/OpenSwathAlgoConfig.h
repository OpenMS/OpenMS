
#ifndef OPENSWATHALGO_DLLAPI_H
#define OPENSWATHALGO_DLLAPI_H

#ifdef OPENSWATHALGO_STATIC_DEFINE
#  define OPENSWATHALGO_DLLAPI
#  define OPENSWATHALGO_NO_EXPORT
#else
#  ifndef OPENSWATHALGO_DLLAPI
#    ifdef OpenSwathAlgo_EXPORTS
        /* We are building this library */
#      define OPENSWATHALGO_DLLAPI __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define OPENSWATHALGO_DLLAPI __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef OPENSWATHALGO_NO_EXPORT
#    define OPENSWATHALGO_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef OPENSWATHALGO_DEPRECATED
#  define OPENSWATHALGO_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef OPENSWATHALGO_DEPRECATED_EXPORT
#  define OPENSWATHALGO_DEPRECATED_EXPORT OPENSWATHALGO_DLLAPI OPENSWATHALGO_DEPRECATED
#endif

#ifndef OPENSWATHALGO_DEPRECATED_NO_EXPORT
#  define OPENSWATHALGO_DEPRECATED_NO_EXPORT OPENSWATHALGO_NO_EXPORT OPENSWATHALGO_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef OPENSWATHALGO_NO_DEPRECATED
#    define OPENSWATHALGO_NO_DEPRECATED
#  endif
#endif

#endif
