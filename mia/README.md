# mia library #

mia is a lightweight, cross-platform C++ image processing library.

## Dependencies ##

mia depends on the [Eigen library](eigen.tuxfamily.org). In order to be able to find Eigen during compilation, set the environment variable THIRD_PARTY_DIR to the parent folder of the Eigen base folder.

*Example*: If Eigen is located in /home/user/Coding/3rdparty/eigen, one can set the environment variable temporarily via

```
#!bash

$ export THIRD_PARTY_DIR=/home/user/Coding/3rdparty/
```


## Build instructions ##

mia comes with a CMake configuration file. To build the library, you can do the following from the base folder:

```
#!bash

$ mkdir build
$ cd build
$ cmake ..
$ make
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
