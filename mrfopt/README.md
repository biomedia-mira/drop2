# mrfopt library #

mrfopt is a cross-platform C++ library with an easy interface to third-party higher-order MRF optimization methods.

## Dependencies ##

mrfopt depends on third-party code which is included in the ext folder:

* [Kolmogorov's QPBO](http://pub.ist.ac.at/~vnk/software.html)
* [Ishikawa's ELC](http://www.f.waseda.jp/hfs/software.html)

## Build instructions ##

mrfopt comes with a CMake configuration file. To build the library, you can do the following from the base folder:

```
#!bash

$ mkdir build
$ cd build
$ cmake ..
$ make
```