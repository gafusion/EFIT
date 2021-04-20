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

set(BLA_STATIC ON)
find_package(BLAS REQUIRED)
set(math_libs ${BLAS_LIBRARIES} ${math_libs})
if (ENABLE_Lapack)
  find_package(LAPACK REQUIRED)
  set(math_libs ${LAPACK_LIBRARIES} ${math_libs})
endif()
