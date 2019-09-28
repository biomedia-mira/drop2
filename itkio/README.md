# itkio library #

itkio is a cross-platform C++ library for simplified IO of image data using the ITK library.

## Dependencies ##

itkio depends on several libraries which need to be installed.

* [Mia library](https://bitbucket.org/bglocker/mia)
* [ITK](http://itk.org)
* [Boost](http://www.boost.org/) (components [system](http://www.boost.org/doc/libs/1_57_0/libs/system/doc/index.html), [filesystem](http://www.boost.org/doc/libs/1_57_0/libs/filesystem/doc/index.htm))

In order to be able to find Eigen and ITK during compilation, set the environment variable THIRD_PARTY_DIR to the parent folder containing both Eigen and ITK.

*Example*: If Eigen and ITK are located in /home/<user>/Coding/3rdparty, one can set the environment variable temporarily via

```
#!bash

$ export THIRD_PARTY_DIR=/home/<user>/Coding/3rdparty/
```

## Build instructions ##

iktio comes with a CMake configuration file. In order to be able to find the mia library, you should have both mia and itkio under the same parent folder. In this parent folder, you need to create a CMakeLists.txt file containing the following lines:

```
#!text

cmake_minimum_required(VERSION 3.0)

project(common)

add_subdirectory(mia)
add_subdirectory(itkio)

```

From the folder where the CMakeLists.txt is located, you can then do the following to build itkio and mia:

```
#!bash

$ mkdir build
$ cd build
$ cmake ..
$ make

```

## Building ITK ##

Here is an example for how to build ITK. Download the sources [here](https://itk.org/ITK/resources/software.html) and unzip in your 3rdparty folder (e.g., /home/<user>/Coding/3rdparty).

From the ITK folder, do the following:

```
#!bash

$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/home/<user>/Coding/3rdparty/itk ..
$ make -j4
$ make install

```

## License ##

Copyright 2016 Ben Glocker ([b.glocker@imperial.ac.uk](mailto:b.glocker@imperial.ac.uk))

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.