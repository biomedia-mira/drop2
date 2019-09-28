# drop2

drop2 is a cross-platform C++ intensity-based image registration tool. It works for 2D and 3D images and has both linear and non-linear registration methods implemented.

Linear (i.e., rigid and affine) registration is based on Downhill simplex optimization.

Non-linear registration is based on free form deformations with efficient discrete MRF optimization as introduced in:

```
@article{glocker2008dense,
  title={Dense image registration through MRFs and efficient linear programming},
  author={Glocker, Ben and Komodakis, Nikos and Tziritas, Georgios and Navab, Nassir and Paragios, Nikos},
  journal={Medical image analysis},
  volume={12},
  number={6},
  pages={731--741},
  year={2008},
  publisher={Elsevier}
}
```

If you make use of the code, please cite this paper in any resulting publications.

## Dependencies ##

drop2 depends on several third-party libraries:

* [Eigen](eigen.tuxfamily.org)
* [Intel TBB](https://www.threadingbuildingblocks.org/)
* [Boost](http://www.boost.org/)
* [ITK](http://itk.org)

## Build instructions ##

Make sure all dependencies are installed, and download and install [ITK](http://itk.org) and [Eigen](http://eigen.tuxfamily.org):

```
sudo apt-get install libboost-all-dev libtbb-dev
```

drop2 comes with a CMake configuration file. From the folder where the CMakeLists.txt is located (same as this README), do the following to build drop2:

```
#!bash

$ mkdir build
$ cd build
$ export THIRD_PARTY_DIR=<folder_containing_eigen>
$ cmake ..
$ make

```

The actual registration tool is called dropreg and will be located in build/drop/apps/dropreg.
