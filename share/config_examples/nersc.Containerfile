FROM docker.io/ubuntu:jammy

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /opt

RUN \
    apt-get update        && \
    apt-get upgrade --yes && \
    apt-get install --yes    \
        build-essential      \
        cmake                \
        gfortran             \
        libblas-dev          \
        liblapack-dev        \
        libcurl4-gnutls-dev  \
        libxml2-dev          \
        m4                   \
        pkg-config           \
        python3-dev          \
        python3-pip          \
        wget                 \
        vim              &&  \
    apt-get clean all    &&  \
    rm -rf /var/lib/apt/lists/*

#install mpich
ARG mpich=4.0.2
ARG mpich_prefix=mpich-$mpich
ENV MPI_DIR=/opt/mpi
RUN \
    wget https://www.mpich.org/static/downloads/$mpich/$mpich_prefix.tar.gz && \
    tar xvzf $mpich_prefix.tar.gz                                           && \
    cd $mpich_prefix                                                        && \
    FFLAGS="-fPIC -fallow-argument-mismatch" FCFLAGS="-fallow-argument-mismatch" CPPFLAGS="-fPIC" ./configure --prefix=$MPI_DIR && \
    make -j 16                                                              && \
    make install                                                            && \
    make clean                                                              && \
    cd ..                                                                   && \
    rm -rf $mpich_prefix

ENV PATH=$MPI_DIR/bin:$PATH
ENV LD_LIBRARY_PATH=$MPI_DIR/lib:$LD_LIBRARY_PATH

RUN /sbin/ldconfig

#install hdf5
ENV HDF5_DIR=/opt/hdf5
ARG h5_version=hdf5-1_8_23
RUN wget https://github.com/HDFGroup/hdf5/archive/refs/tags/$h5_version.tar.gz && \
    tar xvzf $h5_version.tar.gz                                                && \
    cd hdf5-$h5_version                                                        && \
    FFLAGS="-fPIC" CPPFLAGS="-fPIC" ./configure \
    --prefix=$HDF5_DIR --enable-fortran --enable-static --enable-shared --with-pic && \
    make && \
    make install

ENV LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH

#install netcdf c
ENV NETCDF_DIR=/opt/netcdf
ARG netcdf_c_version=4.8.0
RUN wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v$netcdf_c_version.tar.gz && \
    tar xvzf v$netcdf_c_version.tar.gz                                           && \
    cd netcdf-c-$netcdf_c_version                                                && \
    ./configure --prefix=$NETCDF_DIR --disable-hdf5 --disable-netcdf-4 && \
    make && \
    make install

ENV LD_LIBRARY_PATH=$NCDIR/lib:$LD_LIBRARY_PATH

#install netcdf fortran
ARG netcdf_fortran_version=4.5.3
RUN wget https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v$netcdf_fortran_version.tar.gz && \
    tar xvzf v$netcdf_fortran_version.tar.gz                                           && \
    cd netcdf-fortran-$netcdf_fortran_version                                                && \
    CPPFLAGS=-I${NETCDF_DIR}/include LDFLAGS=-L${NETCDF_DIR}/lib ./configure --prefix=$NETCDF_DIR --disable-fortran-type-check --disable-netcdf-4 && \
    make && \
    make install

ENV LD_LIBRARY_PATH=$NFDIR/lib:$LD_LIBRARY_PATH

RUN /sbin/ldconfig

ENV EFIT_ROOT=/opt/efit-main

ADD efit /opt/efit

ARG mpicmd_string="srun -n"
RUN mkdir -p $EFIT_ROOT/build \
        && cd $EFIT_ROOT/build \
        && cmake \
            -DCMAKE_C_COMPILER=mpicc \
            -DCMAKE_Fortran_COMPILER=mpifort \
            -DENABLE_PARALLEL:BOOL=ON \
            -DMPICMD:STRING="$mpicmd_string" \
            -DBLAS_LIBRARIES:PATH="/usr/lib/x86_64-linux-gnu/blas/libblas.a" \
            -DLAPACK_LIBRARIES:PATH="/usr/lib/x86_64-linux-gnu/lapack/liblapack.a" \
            -DENABLE_NETCDF:BOOL=ON \
            -DNetCDF_DIR=$NETCDF_DIR \
            -DENABLE_HDF5:BOOL=ON \
            -DHDF5_ROOT:PATH=$HDF5_DIR \
            -DCMAKE_BUILD_TYPE:STRING=Release \
            -DTEST_EFUND:BOOL=False \
            -DCMAKE_BUILD_TYPE:STRING=Debug \
            /opt/efit

RUN cd $EFIT_ROOT/build && make