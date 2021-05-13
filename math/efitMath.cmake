######################################################################
#
# CMake to find solver dependencies
#
######################################################################

set(ENABLE_Lapack FALSE CACHE BOOL "Whether to link with Lapack")

set(math_libs "")

#---------------------------------------------------------------------
#
# Find BLAS / LAPACK
# https://cmake.org/cmake/help/latest/module/FindBLAS.html
#
#---------------------------------------------------------------------

set(BLAS_STATIC ON)
if (EXISTS ${BLAS_LIBRARIES})
    set(HAVE_BLAS TRUE)
else()
    set(HAVE_BLAS FALSE)
    find_package(BLAS REQUIRED)
    if (ENABLE_BLAS)
        set(HAVE_BLAS TRUE)
    endif()
endif()
set(math_libs ${BLAS_LIBRARIES} ${math_libs})

if (EXISTS ${LAPACK_LIBRARIES})
    set(HAVE_LAPACK TRUE)
else()
    set(HAVE_LAPACK FALSE)
    find_package(LAPACK REQUIRED)
endif()
set(math_libs ${LAPACK_LIBRARIES} ${math_libs})
