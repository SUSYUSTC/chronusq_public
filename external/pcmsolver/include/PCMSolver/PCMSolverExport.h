
#ifndef PCMSolver_EXPORT_H
#define PCMSolver_EXPORT_H

#ifdef PCMSolver_STATIC_DEFINE
#  define PCMSolver_EXPORT
#  define PCMSolver_NO_EXPORT
#else
#  ifndef PCMSolver_EXPORT
#    ifdef PCMSolver_EXPORTS
        /* We are building this library */
#      define PCMSolver_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define PCMSolver_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef PCMSolver_NO_EXPORT
#    define PCMSolver_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef PCMSolver_DEPRECATED
#  define PCMSolver_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef PCMSolver_DEPRECATED_EXPORT
#  define PCMSolver_DEPRECATED_EXPORT PCMSolver_EXPORT PCMSolver_DEPRECATED
#endif

#ifndef PCMSolver_DEPRECATED_NO_EXPORT
#  define PCMSolver_DEPRECATED_NO_EXPORT PCMSolver_NO_EXPORT PCMSolver_DEPRECATED
#endif

#if 1 /* DEFINE_NO_DEPRECATED */
#  ifndef PCMSolver_NO_DEPRECATED
#    define PCMSolver_NO_DEPRECATED
#  endif
#endif

#endif
