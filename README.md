# drop2

*drop2* is a cross-platform C++ intensity-based image registration tool. It works for 2D and 3D images and has both linear and non-linear registration methods implemented.

*Linear* (i.e., rigid and affine) registration is based on Downhill simplex optimization. *Non-linear* (i.e., deformable) registration is based on free form deformations with an efficient discrete MRF optimization scheme as introduced in:

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

@article{glocker2011deformable,
  title={Deformable medical image registration: setting the state of the art with discrete methods},
  author={Glocker, Ben and Sotiras, Aristeidis and Komodakis, Nikos and Paragios, Nikos},
  journal={Annual review of biomedical engineering},
  volume={13},
  pages={219--244},
  year={2011},
  publisher={Annual Reviews}
}
```

If you make use of the code, please cite one of these papers in any resulting publications.

## Dependencies

*drop2* depends on several third-party libraries:

* [Eigen](eigen.tuxfamily.org)
* [Intel TBB](https://www.threadingbuildingblocks.org/) (tested up to v.4.4)
* [Boost](http://www.boost.org/) (tested up to v1.58)
* [ITK](http://itk.org) (tested up to [v4.13.2](https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz)

## Build instructions

Make sure all dependencies are installed, and download and install [ITK](https://itk.org/Wiki/ITK/Getting_Started/Build/Linux) and [Eigen](http://eigen.tuxfamily.org).

You can download and install ITK via:

```
#!bash

$ mkdir 3rdparty
$ cd 3rdparty
$ wget https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz
$ tar xvf InsightToolkit-4.13.2.tar.gz
$ cd InsightToolkit-4.13.2
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=../../itk ..
$ make -j4
```

You can install Boost and TBB via `apt-get`:

```
#!bash

$ sudo apt-get install libboost-all-dev libtbb-dev
```

Note, you might have to specify a specific version via `apt-get install <package>=<version>`.

*drop2* comes with a CMake configuration file. From the top folder where `CMakeLists.txt` is located (same as this README), do the following to build all internal libraries and executables:

```
#!bash

$ mkdir build
$ cd build
$ export THIRD_PARTY_DIR=<folder_containing_eigen_and_itk>
$ cmake ..
$ make -j4

```

The actual registration tool is called `dropreg` and will be located in `build/drop/apps/dropreg`.

## Usage

Note, *drop2* is implemented with the intention of making image registration easy-to-use. Some internal optimization parameters are hard-coded and cannot be changed via command line arguments. Run `./dropreg -h` to see a list of arguments.

Detailed instructions and examples are coming soon...
