#ifndef CONFIG_H
#define CONFIG_H


#cmakedefine USEMPI
#cmakedefine USE_NETCDF
#cmakedefine USE_HDF5
#cmakedefine USE_SNAP
#cmakedefine USE_MSE
#cmakedefine USE_MDS
#cmakedefine MPI_THREAD_FUNNELED
#cmakedefine OBJ_MEM_PROF
#cmakedefine HAVE_OPENMP
#cmakedefine HAVE_BLAS
#cmakedefine HAVE_LAPACK
#cmakedefine HAVE_NETCDF
#cmakedefine HAVE_SNAP
#cmakedefine HAVE_MSE
#cmakedefine HAVE_MDSPLUS
#cmakedefine TIME_LEVEL1
#cmakedefine TIME_LEVEL2
#cmakedefine DEBUG_LEVEL1
#cmakedefine DEBUG_LEVEL2
#cmakedefine DEBUG_LEVEL3
#cmakedefine DEBUG_MSELS
#cmakedefine DEBUG_PLTS
#cmakedefine PROJECT_REV ${PROJECT_REV}
#cmakedefine PROJECT_URL ${PROJECT_URL}
#cmakedefine __cray
#cmakedefine __ifort
#cmakedefine __gfortran
#cmakedefine __pgi
#cmakedefine __xlf

#endif 
