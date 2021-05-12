    #!/bin/tcsh
    module load cmake
    module load env/gcc9.2

    cmake \
    -DHAVE_NETCDF:BOOL=TRUE \
    -DNETCDF_DIR:PATH=/fusion/usc/opt/env/gcc9.2 \
    ..
