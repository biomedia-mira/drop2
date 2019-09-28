# drop library #

drop is a cross-platform C++ library for linear and nonlinear image registration.

## Dependencies ##

drop depends on several libraries which need to be installed.

* [mia library](https://bitbucket.org/bglocker/mia)
* [mrfopt library](https://bitbucket.org/bglocker/mrfopt)
* [Intel TBB](https://www.threadingbuildingblocks.org/)

## Build instructions ##

drop comes with a CMake configuration file. In order to be able to find the mia and mrfopt library, you should have mia, mrfopt and drop under the same parent folder. In this parent folder, you need to create a CMakeLists.txt file containing the following lines:

```
#!text

cmake_minimum_required(VERSION 3.0)

project(common)

add_subdirectory(mia)
add_subdirectory(mrfopt)
add_subdirectory(drop)

```

From the folder where the CMakeLists.txt is located, you can then do the following to build drop and mia:

```
#!bash

$ mkdir build
$ cd build
$ cmake ..
$ make

```