# drop2

[![Build Status](https://travis-ci.org/biomedia-mira/drop2.svg?branch=master)](https://travis-ci.org/biomedia-mira/drop2)

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

If you make use of *drop2*, it would be great if you cite one of these papers in any resulting publications.

## Dependencies

*drop2* depends on several third-party libraries:

* [Eigen](eigen.tuxfamily.org)
* [Intel TBB](https://www.threadingbuildingblocks.org/) (tested up to v.4.4)
* [Boost](http://www.boost.org/) (tested up to v1.58)
* [ITK](http://itk.org) (tested up to [v4.13.2](https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz))

## Build instructions

Eigen is a header-only library and can be simply installed via:

```
#!bash

$ mkdir 3rdparty
$ cd 3rdparty
$ wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
$ mkdir -p eigen && tar xvf 3.3.7.tar.gz -C eigen --strip-components=1
```

You can download and install ITK in the same `3rdparty` folder via:

```
#!bash

$ wget https://sourceforge.net/projects/itk/files/itk/4.13/InsightToolkit-4.13.2.tar.gz
$ tar xvf InsightToolkit-4.13.2.tar.gz
$ cd InsightToolkit-4.13.2
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=../../itk ..
$ make -j4
$ make install
```

Alternatively, you can check out these [ITK install instructions](https://itk.org/Wiki/ITK/Getting_Started/Build/Linux).

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

### Command line arguments

#### General arguments

`-h [--help]`     print list of arguments

*drop2* currently supports the following image file formats: uncompressed and compressed NIfTI (default), JPEG, PNG, DICOM, Meta Image, Analyze, NRRD.

`-s [--source]`   filename of source image
`-t [--target]`   filename of target image
`-o [--output]`   filename of output image (basename is also used for transformation files)

Binary image masks can be provided which restrict the calculation of similarity measures to the non-zero mask points.

`--smask`   filename of source mask
`--tmask`   filename of target mask

The affine transformation can be initialized by providing a text file with a 4-by-4 matrix.

`--transform`   filename of transformation file

Same for displacement fields, which can be initialized by providing three images containing the displacement components.

`--fx`  filename of displacement field for X-component
`--fy`  filename of displacement field for Y-component
`--fz`  filename of displacement field for Z-component

The same tool can be used for 2D and 3D images. For 2D images, an argument needs to be added. Note, the resulting transformations will be still in 3D format, with the third dimension being set to identiy transformations (or zero displacements).

`--mode2d`

The interpolation method being used to generate the resulting warped output image can be selected.

`--ointerp`   0=NEAREST, 1=LINEAR (default)

The fill value for out of bounds interpolation can be set.

`--ofill`   default is 0.0

If no output image is required, this can be disabled.

`--onoimage`

The affine transformation can be composed into the resulting displacement field. This may be useful if the displacement fields are used to warp landmarks or similar after registration. The default is to have the affine component of the overall transformation separately in a text file, and displacement field only contains the non-linear component from the Free Form Deformation.

`--ocompose`

#### Center-of-mass alignment

*drop2* provides an optional initialization method using intensity center-of-mass alignment. This is particularly effective for brain image registration.

`-c [--com]`

#### Linear registration

The linear registration of *drop2* is based on gradient-free optimization of similarity measures using Downhill Simplex (aka [Neader-Mead Simplex](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)).

Linear registration is enabled by adding the argument

`-l [--linear]`

There are three types of linear transformations, rigid, simiarity (rigid plus isotropic scaling), and fully affine.

`--ltype` 0=RIGID (default), 1=SIMILARITY, 2=AFFINE

*drop2* comes with three different intensity-based similarity measures, mean absolute differences (MAD), correlation coefficient (CC), and entropy correlation coefficient (ECC). The later is a variant of normalized mutual information and may work for multi-modal image registration. ECC is computed from joint histograms with 16 bins.

`--lsim`  0=MAD (default), 1=CC, 2=ECC

An important set of parameters is concerning the image resolution and multi-scale pyramid on which registration is performed. The pyramid is configured with triplets of scaling factors. For example, a three-level registration with image downscaling factors of 4, 2, and 1 is done with

`--llevels 4 4 4 2 2 2 1 1 1`

Note, a scaling factor of 1 means full, native image resolution. Factors can be different for each image dimension. For 2D, the last number in the triplet should be 1. For example `--llevels 8 8 1 4 4 1` would run a two-level registration on images downsampled by a factor of 8, and then 4.

The number iterations per level of the Downhill Simplex optimizer can be set.

`--liters`  default is 100

The interpolation method for computing the similarity measure can be set. Using nearest neighbors will be speed up registration, but may be less accurate.

`--linterp` 0=NEAREST, 1=LINEAR (default)

The linear registration will use random subsampling for calculating the similarity measure. The sampling rate controls the number of image points being used and directly influences the efficiency of the registration.

`--lsampling` default is 0.05 (5% of the total image points)

#### Non-linear registration

The non-linear registration in *drop2* is based on discrete MRF optimization of control point displacements in a cubic B-spline FFD transformation model.

Non-linear registration is enabled by adding the argument

`-n [--nonlinear]`

To be done...

### Examples

Example configurations for common registration problems are coming soon...

## Acknowledgements

The original idea for using discrete MRF optimization for image registration was developed back in 2006 by [Ben Glocker](http://www.doc.ic.ac.uk/~bglocker/) and [Nikos Paragios](https://en.wikipedia.org/wiki/Nikos_Paragios). The core algorithm has been patented (Pub. No. [WO/2009/010860](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2009010860)).

We are very grateful for the support and contribution by Jiří Borovec and Jan Kybic. Their efforts on the [Automatic Non-rigid Histological Image Registration (ANHIR) challenge](https://anhir.grand-challenge.org/) made us push for the release of this updated version of the original *drop* registration method which had been lying around in a private repository for too long.

Special thanks go to [Hauke Heibel](https://github.com/hauke76) who has contributed significantly to the re-implementation of the C++ image processing backend.
