
#ifndef OPENMS_EXPORT_H
#define OPENMS_EXPORT_H

#ifdef OPENMS_STATIC_DEFINE
#  define OPENMS_EXPORT
#  define OPENMS_NO_EXPORT
#else
#  ifndef OPENMS_EXPORT
#    ifdef IsoSpec_EXPORTS
        /* We are building this library */
#      define OPENMS_EXPORT 
#    else
        /* We are using this library */
#      define OPENMS_EXPORT 
#    endif
#  endif

#  ifndef OPENMS_NO_EXPORT
#    define OPENMS_NO_EXPORT 
#  endif
#endif

#ifndef OPENMS_DEPRECATED
#  define OPENMS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef OPENMS_DEPRECATED_EXPORT
#  define OPENMS_DEPRECATED_EXPORT OPENMS_EXPORT OPENMS_DEPRECATED
#endif

#ifndef OPENMS_DEPRECATED_NO_EXPORT
#  define OPENMS_DEPRECATED_NO_EXPORT OPENMS_NO_EXPORT OPENMS_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef OPENMS_NO_DEPRECATED
#    define OPENMS_NO_DEPRECATED
#  endif
#endif

#endif /* OPENMS_EXPORT_H */
