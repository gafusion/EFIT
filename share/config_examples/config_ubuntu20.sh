#  The following packages are required to install and run the code:
#    git (optional)
#    build-essential
#    cmake
#    gfortran

    cmake \
    -DBLAS_LIBRARIES:PATH='/usr/lib/x86_64-linux-gnu/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib/x86_64-linux-gnu/liblapack.a' \
    ..
