
#ifndef SUPERHIRN_DLLAPI_H
#define SUPERHIRN_DLLAPI_H

#ifdef SUPERHIRN_STATIC_DEFINE
#  define SUPERHIRN_DLLAPI
#  define SUPERHIRN_NO_EXPORT
#else
#  ifndef SUPERHIRN_DLLAPI
#    ifdef SuperHirn_EXPORTS
        /* We are building this library */
#      define SUPERHIRN_DLLAPI __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define SUPERHIRN_DLLAPI __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef SUPERHIRN_NO_EXPORT
#    define SUPERHIRN_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef SUPERHIRN_DEPRECATED
#  define SUPERHIRN_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SUPERHIRN_DEPRECATED_EXPORT
#  define SUPERHIRN_DEPRECATED_EXPORT SUPERHIRN_DLLAPI SUPERHIRN_DEPRECATED
#endif

#ifndef SUPERHIRN_DEPRECATED_NO_EXPORT
#  define SUPERHIRN_DEPRECATED_NO_EXPORT SUPERHIRN_NO_EXPORT SUPERHIRN_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SUPERHIRN_NO_DEPRECATED
#    define SUPERHIRN_NO_DEPRECATED
#  endif
#endif

#endif /* SUPERHIRN_DLLAPI_H */
