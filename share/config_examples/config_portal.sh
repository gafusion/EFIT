# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module load gcc/9.3.0 cmake/3.17.1 netcdf-fortran-4.5.2 netcdf-c-4.7.3

    cmake \
    -DMPICMD:STRING='srun -n ' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/usr/pppl/gcc/9.3-pkgs/netcdf-fortran-4.5.2' \
    -DNetCDF_C_LIBRARY='/usr/pppl/gcc/9.3-pkgs/netcdf-c-4.7.3'\
    ..
