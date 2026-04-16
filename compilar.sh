#!/bin/bash

if [ $# -ne 1 ]; then
       echo "Usage: $0 <install_dir>"
fi

RUN_PATH=$1
SOURCE_DIR="$(dirname $0)"

echo "MPICH source code in $SOURCE_DIR"

# get into source dir
pushd .
cd $SOURCE_DIR

echo "Installing mpich runtime in dorectory ${RUN_PATH}..."

# clean previous compilations
make clean 

# configure and compile source code
./configure --enable-threads=multiple --enable-romio --prefix=${RUN_PATH} --with-file-system=ufs+nfs CPPFLAGS=-DUSE_MPI_VERSIONS --enable-g=none --enable-fast=O0 --enable-shared --disable-fortran CC=gcc CXX=g++

make -j8

# fix pkgconfig configuration file
#sed 's/-rpath /-Wl,--no-undefined,-rpath-link=/' src/packaging/pkgconfig/mpich.pc > src/packaging/pkgconfig/mpich.pc.bak
#cp src/packaging/pkgconfig/mpich.pc.bak src/packaging/pkgconfig/mpich.pc
#rm src/packaging/pkgconfig/mpich.pc.bak

make install

popd
