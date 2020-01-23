
#ifndef OPENMS_DLLAPI_H
#define OPENMS_DLLAPI_H

#ifdef OPENMS_STATIC_DEFINE
#  define OPENMS_DLLAPI
#  define OPENMS_NO_EXPORT
#else
#  ifndef OPENMS_DLLAPI
#    ifdef OpenMS_EXPORTS
        /* We are building this library */
#      define OPENMS_DLLAPI __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define OPENMS_DLLAPI __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef OPENMS_NO_EXPORT
#    define OPENMS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef OPENMS_DEPRECATED
#  define OPENMS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef OPENMS_DEPRECATED_EXPORT
#  define OPENMS_DEPRECATED_EXPORT OPENMS_DLLAPI OPENMS_DEPRECATED
#endif

#ifndef OPENMS_DEPRECATED_NO_EXPORT
#  define OPENMS_DEPRECATED_NO_EXPORT OPENMS_NO_EXPORT OPENMS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef OPENMS_NO_DEPRECATED
#    define OPENMS_NO_DEPRECATED
#  endif
#endif

#endif /* OPENMS_DLLAPI_H */
