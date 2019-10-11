
FIND_PATH(MKL_INCLUDE_DIR mkl_cblas.h
    HINTS $ENV{MKL_HOME}/include $ENV{MKLROOT}/include
)

IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
    SET(INTEL_LIBDIR "intel64")
ELSE ( CMAKE_SYSTEM_PROCESSOR  STREQUAL "x86_64" )
    SET(INTEL_LIBDIR "ia32")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

FIND_LIBRARY(MKL_LAPACK
    NAMES mkl_lapack
    HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
)


IF ( MKL_LAPACK_FOUND )
    # old MKL versions
    IF ( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
        FIND_LIBRARY(MKL
            NAMES mkl_ia64
            HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
                  $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
        )
    ELSE ( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
        FIND_LIBRARY(MKL
            NAMES mkl_ia32
            HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
                  $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
        )
    ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

    FIND_LIBRARY(MKL_GUIDE
        NAMES guide
        HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
              $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
    )

    SET (MKL ${MKL} ${MKL_LAPACK} ${MKL_GUIDE} )
ELSE (MKL_LAPACK_FOUND )
    # newer MKL version
    SET (MKL_LAPACK "")
    FIND_LIBRARY(MKL_INTEL
        NAMES mkl_intel_lp64
        HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
              $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
    )
    FIND_LIBRARY(MKL_SEQUENTIAL
        NAMES mkl_sequential
        HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
              $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
    )
    FIND_LIBRARY(MKL_CORE
        NAMES mkl_core
        HINTS $ENV{MKL_HOME}/lib/${INTEL_LIBDIR} 
              $ENV{MKLROOT}/lib/${INTEL_LIBDIR}
    )
    SET (MKL ${MKL_INTEL} ${MKL_SEQUENTIAL} ${MKL_CORE} )
ENDIF ( MKL_LAPACK_FOUND )


SET( MKL_BLAS_INCLUDE_FILE ${MKL_INCLUDE_DIR}/mkl_blas.h )
SET( MKL_LAPACK_INCLUDE_FILE ${MKL_INCLUDE_DIR}/mkl_lapack.h )
GET_FILENAME_COMPONENT(MKL_LIB_DIR ${MKL_CORE} PATH)

IF (MKL_INCLUDE_DIR)
    SET(MKL_FOUND ON)
ENDIF (MKL_INCLUDE_DIR)

IF (MKL_FOUND OR MKL_INTEL_FOUND )
    IF (NOT MKL_FIND_QUIETLY)
        MESSAGE(STATUS "Found MKL: ${MKL_INCLUDE_DIR}")
    ENDIF (NOT MKL_FIND_QUIETLY)
ELSE(MKL_FOUND OR MKL_INTEL_FOUND )
    IF (MKL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find MKL")
    ENDIF (MKL_FIND_REQUIRED)
ENDIF (MKL_FOUND OR MKL_INTEL_FOUND )

MARK_AS_ADVANCED(MKL_INCLUDE_DIR)
MARK_AS_ADVANCED(MKL_LIB_DIR)
MARK_AS_ADVANCED(MKL_LAPACK)
MARK_AS_ADVANCED(MKL)
MARK_AS_ADVANCED(MKL_GUIDE)
MARK_AS_ADVANCED(MKL_INTEL)
MARK_AS_ADVANCED(MKL_SEQUENTIAL)
MARK_AS_ADVANCED(MKL_CORE)