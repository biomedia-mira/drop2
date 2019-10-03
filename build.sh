#!/bin/bash

ROOT=$(pwd)
# prepare 3-party libs
mkdir 3rdParty
cd 3rdParty
export THIRD_PARTY_DIR=$(pwd)

# get Eigen
cd $THIRD_PARTY_DIR
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz --progress=bar:force:noscroll
mkdir eigen
tar xf 3.3.7.tar.gz -C eigen --strip-components=1

# download and install ITK
cd $THIRD_PARTY_DIR
wget https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz --progress=bar:force:noscroll
tar xf InsightToolkit-4.13.2.tar.gz
cd InsightToolkit-4.13.2
mkdir build
cd build
cmake \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_INSTALL_PREFIX=$THIRD_PARTY_DIR/itk \
   ..
make -j$(nproc)
make install

# install last Boost libs
apt-get install libboost-all-dev libtbb-dev

# building the
cd $ROOT
mkdir build
cd build
cmake ..
make -j$(nproc)
