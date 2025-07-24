
#ifndef OPENMS_GUI_DLLAPI_H
#define OPENMS_GUI_DLLAPI_H

#ifdef OPENMS_GUI_STATIC_DEFINE
#  define OPENMS_GUI_DLLAPI
#  define OPENMS_GUI_NO_EXPORT
#else
#  ifndef OPENMS_GUI_DLLAPI
#    ifdef OpenMS_GUI_EXPORTS
        /* We are building this library */
#      define OPENMS_GUI_DLLAPI __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define OPENMS_GUI_DLLAPI __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef OPENMS_GUI_NO_EXPORT
#    define OPENMS_GUI_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef OPENMS_GUI_DEPRECATED
#  define OPENMS_GUI_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef OPENMS_GUI_DEPRECATED_EXPORT
#  define OPENMS_GUI_DEPRECATED_EXPORT OPENMS_GUI_DLLAPI OPENMS_GUI_DEPRECATED
#endif

#ifndef OPENMS_GUI_DEPRECATED_NO_EXPORT
#  define OPENMS_GUI_DEPRECATED_NO_EXPORT OPENMS_GUI_NO_EXPORT OPENMS_GUI_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef OPENMS_GUI_NO_DEPRECATED
#    define OPENMS_GUI_NO_DEPRECATED
#  endif
#endif

#endif /* OPENMS_GUI_DLLAPI_H */
