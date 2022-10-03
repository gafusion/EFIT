#!/bin/sh

dockerdir=`dirname $0`
sourcedir=`dirname $dockerdir`

rm -rf CMake*

    cmake \
        -DCMAKE_INSTALL_PREFIX:PATH=/software \
        -DCMAKE_BUILD_TYPE:STRING=Debug \
        -DCMAKE_COLOR_MAKEFILE:BOOL=TRUE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -DENABLE_DOCS:BOOL=TRUE \
        $sourcedir
