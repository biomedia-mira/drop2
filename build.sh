#!/bin/bash

ITK_VERSION="5.0.0"
EIGEN_VERSION="3.3.7"

ROOT=$(pwd)
# prepare 3-party libs
export THIRD_PARTY_DIR=$ROOT/3rdParty
mkdir $THIRD_PARTY_DIR

# get Eigen
cd $THIRD_PARTY_DIR
wget https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz --progress=bar:force:noscroll
mkdir eigen
tar xf eigen-${EIGEN_VERSION}.tar.gz -C eigen --strip-components=1

# download and install ITK
cd $THIRD_PARTY_DIR
wget https://sourceforge.net/projects/itk/files/itk/${ITK_VERSION%.*}/InsightToolkit-${ITK_VERSION}.tar.gz --progress=bar:force:noscroll
tar xf InsightToolkit-${ITK_VERSION}.tar.gz
cd InsightToolkit-${ITK_VERSION}
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
apt-get install -y libboost-all-dev libtbb-dev

# building the
cd $ROOT
mkdir build
cd build
cmake ..
make -j$(nproc)

unset THIRD_PARTY_DIR
