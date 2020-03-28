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

IMPORTANT: *drop2* uses the image position and orientation information from the image headers. This information is used to create the initial (identity) transformation. The estimated transformation is thus the residual transformation that is required to align the images after taking the internal image transformations into account. This is essential for registering images acquired in the same session (e.g., different MR sequences, pre- and post-contrast, etc.) where the internal transformations typically already pre-align the images very well (and remaining misalignment may come from patient motion or breathing).

When images from different sessions or different patients are to be registered, the internal transformations stored in the image headers are meaningless and an initialization such as center of mass alignment (see below) is required to obtain a rough pre-alignment.

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

*drop2* provides an optional initialization method using intensity center-of-mass alignment. This is particularly effective for inter-subject or multi-modal brain image registration.

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

Note, a scaling factor of 1 means full, native image resolution. Factors can be different for each image dimension. For 2D, the last number in the triplet should be 1. For example `--llevels 8 8 1 4 4 1` would run a two-level registration on images downsampled by a factor of 8, and then 4. There is hard coded limit on the number of image points, and an image won't be resampled if the number of voxels would fall below 32.

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

There is an option to choose between first-order and second-order MRFs. First-order MRFs use pairwise potentials for the regularization term, while second-order MRFs use triple cliques. The latter is experimental and the default is to use first-order MRFs.

`--ntype`   0=FIRST_ORDER (default), 1=SECOND_ORDER

The regularization mode can be set to fluid or elastic. Fluid regularization penalizes only updates of the displacement field, while elastic penalizes the overall displacements. For first-order MRFs, pairwise potentials penalize the L2 norm of the difference vector between neighboring control point displacements (an approximation to first order derivates). In second-rder MRFs, triple clique potentials penalize an approximation to the second order derivates of the displacement field.

`--nreg`  0=FULID (default), 2=ELASTIC

An important parameter is the initial control point spacing of the Free Form Deformation grid which is set in millimeters.

`--nffd`  default is 80

There are the same three similarity measures available for non-linear registration as for linear registration.

`--nsim`  0=MAD (default), 1=CC, 2=ECC

See comments for linear registration on how to set the image levels. In addition, for the case of non-linear registration the image levels are connected to an increase in resolution of the FFD control point grid. Starting from the initial FFD, the control point spacing is halved on each subsequent image level. Note, it is possible to run several levels of FFD spacings on the same image resolution, for example, by setting `--nlevels 2 2 2 2 2 2 2 2 2`. Assuming an initial control points spacing of 80, this would run a three level registration on two times downsampled images for 80, 40 and 20 mm FFD.

`--nlevels`

The number of iterative MRF optimizations per level can be set.

`--niters`  default is 10

The interpolation method for computing the similarity measure can be set. Using nearest neighbors will be speed up registration, but may be less accurate.

`--ninterp` 0=NEAREST, 1=LINEAR (default)

The non-linear registration can use random subsampling for calculating the similarity measure, but this will degrade the quality of the local similarity measures.

`--nsampling` default is 1.0 (100% of the total image points)

The weigthing of the regularization term needs to be set, and this is sensitive to the used similarity measure. For CC and ECC, try values between [0.1,1.0], for MAD the weighting may need to be in the hundreds or more.

`--nlambda`   default is 0.0 (no regularization)

In some applications it may be beneficial to explicitly constraint the deformations by pinning a virtual set of boundary control points. Pinning means to enfource the control point displacements outside the image domain to be zero. This has a strong regularization effect on the resulting deformation. This can be enabled by adding the argument

`--npin`

### Examples

Example configurations for common registration problems:

#### Inter-subject linear registration of 3D brain MRI

Let's say we want to rigidly register a subject's brain MRI (source) to an MNI atlas (target). The following setting should provide a good starting point:

`./dropreg -s <source_fname> -t <target_fname> -o <output_fname> -c -l --lsim 1 --ltype 0 --llevels 4 4 4 2 2 2 --lsampling 0.1`

Note, we are using `-c` to initialize the transformation with center of mass alignment. If we want to do affine registration, it may help to add another level to the image pyramid and maybe increase the sampling rate for the similarity measure calculation (due to the larger set of transformation parameters):

`./dropreg -s <source_fname> -t <target_fname> -o <output_fname> -c -l --lsim 1 --ltype 2 --llevels 6 6 6 4 4 4 2 2 2 --lsampling 0.2`

#### Intra-subject non-rigid registration of 3D abdominal CT scans

Say we want to register abdominal images at different breathing cycles (assuming the images are acquired in the same session):

`./dropreg -s <source_fname> -t <target_fname> -o <output_fname> -n --nsim 1 --nffd 80 --nlevels 4 4 4 4 4 4 2 2 2 --nlambda 0.5`

This is using an FFD control grid pyrmaid with 80, 40 and 20mm spacing. We may need to fine-tune the `nlambda` parameter.


#### Mono-modal non-registration of anisotropic 3D data

In cases where the image data is highly anistropic, we may want to apply different downsampling factors for different dimensions. The following is an example for registration data from the [POPI dataset](https://www.creatis.insa-lyon.fr/rio/dir_validation_data):

`./dropreg -s <source_fname> -t <target_fname> -o <output_fname> -n --nsim 0 --nffd 40 --nlambda 100 --nlevels 4 4 2 2 2 1`

#### Non-rigid 2D registration of histopathology images

This is assuming we have very large images of say 16k by 16k pixels:

`./dropreg --mode2d -s <source_fname> -t <target_fname> -o <output_fname> -l --llevels 4 4 1 2 2 1 --lsim 1 --lsampling 0.2 -n --nsim 1 --nffd 1000 --nlevels 4 4 1 2 2 1 2 2 1 --nlambda 0.5 --mode2d --ocompose`

The argument `ocompose` will result in the affine transformation being composed with the non-rigid deformation and the resulting displacement field will reflect the entire transformation (not just the non-rigid part).

## Acknowledgements

The original idea for using discrete MRF optimization for image registration was developed back in 2006 by [Ben Glocker](http://www.doc.ic.ac.uk/~bglocker/) and [Nikos Paragios](https://en.wikipedia.org/wiki/Nikos_Paragios). The core algorithm has been patented (Pub. No. [WO/2009/010860](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2009010860)).

We are very grateful for the support and contribution by Jiří Borovec and Jan Kybic. Their efforts on the [Automatic Non-rigid Histological Image Registration (ANHIR) challenge](https://anhir.grand-challenge.org/) made us push for the release of this updated version of the original *drop* registration method which had been lying around in a private repository for too long.

Special thanks go to [Hauke Heibel](https://github.com/hauke76) who has contributed significantly to the re-implementation of the C++ image processing backend.
