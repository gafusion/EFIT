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
set(HAVE_BLAS FALSE)
find_package(BLAS REQUIRED)
if (ENABLE_BLAS)
    set(HAVE_BLAS TRUE)
endif()
set(math_libs ${BLAS_LIBRARIES} ${math_libs})
set(HAVE_LAPACK FALSE)
find_package(LAPACK REQUIRED)
set(math_libs ${LAPACK_LIBRARIES} ${math_libs})
