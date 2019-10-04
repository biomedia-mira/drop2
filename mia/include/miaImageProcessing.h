/*
* mia - a lightweight, cross-platform C++ image processing library.
*
* Copyright 2016 Ben Glocker <b.glocker@imperial.ac.uk>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#pragma once

#include "miaImage.h"
#include "miaHistogram.h"

namespace mia
{
  /**
   * \brief Sets all image points to zero.
   * \param image The image.
   **/
  void zeros(Image& image);

  /**
   * \brief Sets all image points to one.
   * \param image The image.
   **/
  void ones(Image& image);

  /**
   * \brief Sets all image points to a specified value.
   * \param image Image.
   * \param value The value.
   **/
  void fill(Image& image, float value);

  /**
   * \brief Sets all image boundary points to specified value.
   * \param image Image.
   * \param value The value.
   **/
  void fill_boundaries(Image& image, float value);

  /**
   * \brief Copies image values from one image to another.
   * \param src The source image.
   * \param dst The destination image.
   **/
  void copy(const Image& src, Image& dst);

  /**
   * \brief Copies image values from a specified region from one image to another.
   * \param src The source image.
   * \param dst The destination image.
   * \param srcX Minimum x-coordinate of the region in the source image.
   * \param srcY Minimum y-coordinate of the region in the source image.
   * \param srcZ Minimum z-coordinate of the region in the source image.
   * \param dstX Minimum x-coordinate of the region in the destination image.
   * \param dstY Minimum y-coordinate of the region in the destination image.
   * \param dstZ Minimum z-coordinate of the region in the destination image.
   * \param sizeX Region size along x-axis.
   * \param sizeY Region size along y-axis.
   * \param sizeZ Region size along z-axis.
   **/
  void copy_region(const Image& src, Image& dst, int srcX, int srcY, int srcZ, int dstX, int dstY, int dstZ, int sizeX, int sizeY, int sizeZ);

  /**
   * \brief Returns the minimum image value.
   * \param image The image.
   * \return The minimum value.
   **/
  float min(const Image& image);

  /**
   * \brief Returns the maximum image value.
   * \param image The image.
   * \return The maximum value.
   **/
  float max(const Image& image);

  /**
   * \brief Returns the mean image value.
   * \param image The image.
   * \return The mean value.
   **/
  float mean(const Image& image);

  /**
  * \brief Returns the variance.
  * \param image The image.
  * \return The variance.
  **/
  float var(const Image& image);

  /**
   * \brief Adds a specified value to all image points.
   * \param image The image.
   * \param value The value.
   **/
  void add(Image& image, float value);

  /**
  * \brief Subtracts a specified value to all image points.
  * \param image The image.
  * \param value The value.
  **/
  void sub(Image& image, float value);

  /**
   * \brief Multiplies all image points with a specified value.
   * \param image The image.
   * \param value The value.
   **/
  void mul(Image& image, float value);

  /**
   * \brief Divides all image points with a specified value.
   * \param image The image.
   * \param value The value.
   **/
  void div(Image& image, float value);

  /**
   * \brief Adds two images pointwise, C=A+B.
   * \param imageA The input image A.
   * \param imageB The input image B.
   * \param imageC The output image C.
   **/
  void add(const Image& imageA, const Image& imageB, Image& imageC);

  /**
   * \brief Subtracts two images pointwise, C=A-B.
   * \param imageA The input image A.
   * \param imageB The input image B.
   * \param imageC The output image C.
   **/
  void sub(const Image& imageA, const Image& imageB, Image& imageC);

  /**
   * \brief Multiplies two images pointwise, C=A*B.
   * \param imageA The input image A.
   * \param imageB The input image B.
   * \param imageC The output image C.
   **/
  void mul(const Image& imageA, const Image& imageB, Image& imageC);

  /**
   * \brief Divides two images pointwise, C=A/B.
   * \param imageA The input image A.
   * \param imageB The input image B.
   * \param imageC The output image C.
   **/
  void div(const Image& imageA, const Image& imageB, Image& imageC);

  /**
   * \brief Computes an intensity histogram of an image for a given intensity range.
   * \param image The image.
   * \param numBins Number of histogram bins.
   * \param minValue Minimum intensity value considered for the histogram counts.
   * \param maxValue Maximum intensity value considered for the histogram counts.
   * \return The histogram.
   **/
  Histogram histogram(const Image& image, int numBins, float minValue, float maxValue);

  /**
   * \brief Returns a subimage of a given image (without copying the image data, and with correct world coordinates).
   * \param image The image.
   * \param x The first x-coordinate (in image space) of the subimage.
   * \param y The first y-coordinate (in image space) of the subimage.
   * \param z The first z-coordinate (in image space) of the subimage.
   * \param sizeX Size of the subimage along the x-axis.
   * \param sizeY Size of the subimage along the y-axis.
   * \param sizeZ Size of the subimage along the z-axis.
   * \return The subimage.
   **/
  Image subimage(Image& image, int x, int y, int z, int sizeX, int sizeY, int sizeZ);

  /**
   * \brief Resamples a given image with the specified element spacing.
   * \param image The image.
   * \param spacing The new element spacing.
   * \param interpolation The interpolation mode.
   * \param outsideValue The value to be used for out of bounds image points.
   * \param minimumSize The minimum size along each dimension of the new image.
   * \return The resampled image.
   **/
  Image resample(const Image& image, const Eigen::Vector3d& spacing, Interpolation interpolation = Interpolation::LINEAR, float outsideValue = 0.0f, int minimumSize = 1);

  /**
   * \brief Resamples a given input image with the specifications of the given output image.
   * \param input The input image.
   * \param output The output image.
   * \param interpolation The interpolation mode.
   * \param outsideValue The value to be used for out of bounds image points.
   **/
  void resample(const Image& input, Image& output, Interpolation interpolation = Interpolation::LINEAR, float outsideValue = 0.0f);

  /**
   * \brief Warps a given image with the specified transformation. Uses backward warping.
   * \param input The image.
   * \param output The warped output image.
   * \param transform The affine transformation matrix.
   * \param interpolation The interpolation mode.
   * \param outsideValue The value to be used for out of bounds image points.
   **/
  void warp(const Image& input, Image& output, const Eigen::Matrix4d& transform, Interpolation interpolation = Interpolation::LINEAR, float outsideValue = 0.0f);

  /**
   * \brief Warps a given image with the specified displacement field. Uses backward warping.
   * \param input The image.
   * \param output The warped output image.
   * \param fieldX The x-component of the displacement field.
   * \param fieldY The y-component of the displacement field.
   * \param fieldZ The z-component of the displacement field.
   * \param interpolation The interpolation mode.
   * \param outsideValue The value to be used for out of bounds image points.
   **/
  void warp(const Image& input, Image& output, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation = Interpolation::LINEAR, float outsideValue = 0.0f);

  /**
   * \brief Warps a given image with the specified transformation and displacement field. Uses backward warping.
   * \param input The image.
   * \param output The warped output image.
   * \param transform The affine transformation matrix.
   * \param fieldX The x-component of the displacement field.
   * \param fieldY The y-component of the displacement field.
   * \param fieldZ The z-component of the displacement field.
   * \param interpolation The interpolation mode.
   * \param outsideValue The value to be used for out of bounds image points.
   **/
  void warp(const Image& input, Image& output, const Eigen::Matrix4d& transform, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation = Interpolation::LINEAR, float outsideValue = 0.0f);

  /**
  * \brief Composes two displacement fields. Result is stored in field B.
  * \param aX The x-component of the displacement field A.
  * \param aY The y-component of the displacement field A.
  * \param aZ The z-component of the displacement field A.
  * \param bX The x-component of the displacement field B.
  * \param bY The y-component of the displacement field B.
  * \param bZ The z-component of the displacement field B.
  **/
  void compose(const Image& aX, const Image& aY, const Image& aZ, Image& bX, Image& bY, Image& bZ);

  /**
  * \brief Composes the specified transformation with the given displacement field.
  * \param transform The affine transformation matrix.
  * \param fieldX The x-component of the displacement field.
  * \param fieldY The y-component of the displacement field.
  * \param fieldZ The z-component of the displacement field.
  **/
  void compose(const Eigen::Matrix4d& transform, Image& fieldX, Image& fieldY, Image& fieldZ);

  /**
   * \brief Applies a Gaussian smoothing.
   * \param input The input image.
   * \param output The output image.
   * \param sigmaX Sigma (in element spacing units) of the Gaussian kernel along the x-axis.
   * \param sigmaY Sigma (in element spacing units) of the Gaussian kernel along the y-axis.
   * \param sigmaZ Sigma (in element spacing units) of the Gaussian kernel along the z-axis.
   **/
  void gauss(const Image& input, Image& output, double sigmaX, double sigmaY, double sigmaZ);

  /**
   * \brief Applies a median filter. Radius of one equals a 3x3x3 filter.
   * \param input The input image.
   * \param output The output image.
   * \param filterRadiusX Filter radius (in pixels) along the x-axis.
   * \param filterRadiusY Filter radius (in pixels) along the y-axis.
   * \param filterRadiusZ Filter radius (in pixels) along the z-axis.
   **/
  void median(const Image& input, Image& output, int filterRadiusX = 1, int filterRadiusY = 1, int filterRadiusZ = 1);

  /**
   * \brief Applies a 1-D filter on the x-axis.
   * \param input The input image.
   * \param output The output image.
   * \param filter The filter.
   **/
  void filter_x(const Image& input, Image& output, const std::vector<double>& filter);

  /**
  * \brief Applies a 1-D filter on the y-axis.
  * \param input The input image.
  * \param output The output image.
  * \param filter The filter.
  **/
  void filter_y(const Image& input, Image& output, const std::vector<double>& filter);

  /**
   * \brief Applies a 1-D filter on the z-axis.
   * \param input The input image.
   * \param output The output image.
   * \param filter The filter.
   **/
  void filter_z(const Image& input, Image& output, const std::vector<double>& filter);

  /**
   * \brief Checks whether image A and image B are defined on the same image domain in world coordinate space.
   * \param imageA The image A.
   * \param imageB The image B.
   * \return True if image domain is the same.
   **/
  bool same_domain(const Image& imageA, const Image& imageB);

  /**
   * \brief Checks the orientation of the image and whether any of the three axes should be flipped when visualizing the image.
   * \param image The image.
   * \return A vector with three boolean flags, one per image dimension.
   **/
  std::vector<bool> flip_flags(const Image& image);

  /**
  * \brief Reorients the image to the standard coordinate system.
  * \return The reoriented image.
  **/
  Image reorient(const Image& image);

  /**
   * \brief Computes the center of intensity mass (in world coordinates) for a given image.
   * \param image The image.
   * \return The center of mass in world coordinates.
   **/
  Eigen::Vector3d center_of_mass(const Image& image);

  /**
   * \brief Pads an image by a given number of pixels on each side. Copies the image data.
   * \param image The image.
   * \param padX1 Padding on the 'left' side of the x-axis.
   * \param padX2 Padding on the 'right' side of the x-axis.
   * \param padY1 Padding on the 'top' side of the y-axis.
   * \param padY2 Padding on the 'bottom' side of the y-axis.
   * \param padZ1 Padding on the 'front' side of the z-axis.
   * \param padZ2 Padding on the 'back' side of the z-axis.
   * \return The padded image.
   **/
  Image pad(const Image& image, int padX1, int padX2, int padY1, int padY2, int padZ1, int padZ2, float value = 0.0f);

  /**
   * \brief Crops an image by a given number of pixels on each side. Copies the image data.
   * \param image The image.
   * \param cropX1 Cropping on the 'left' side of the x-axis.
   * \param cropX2 Cropping on the 'right' side of the x-axis.
   * \param cropY1 Cropping on the 'top' side of the y-axis.
   * \param cropY2 Cropping on the 'bottom' side of the y-axis.
   * \param cropZ1 Cropping on the 'front' side of the z-axis.
   * \param cropZ2 Cropping on the 'back' side of the z-axis.
   * \return The cropped image.
   **/
  Image crop(const Image& image, int cropX1, int cropX2, int cropY1, int cropY2, int cropZ1, int cropZ2);

  /**
   * \brief Computes an integral image for given intensity image.
   * \param image The image.
   * \return The integral image.
   **/
  Image integral_image(const Image& image);

  /**
   * \brief Alternative method for computing an integral image for given intensity image. Uses temporary images.
   * \param image The image.
   * \return The integral image.
   **/
  Image integral_image_alternative(const Image& image);

  /**
   * \brief Reconstructs an intensity image from a given integral image.
   * \param integral The integral image.
   * \param image Reconstructed intensity image.
   **/
  void reverse_integral_image(const Image& integral, Image& image);

  /**
   * \brief Computes the sum of intensities in a given cuboid using an integral image.
   * \param integral The integral image.
   * \param minX Minimum x-coordinate (in image space) of cuboid.
   * \param minY Minimum y-coordinate (in image space) of cuboid.
   * \param minZ Minimum z-coordinate (in image space) of cuboid.
   * \param maxX Maximum x-coordinate (in image space) of cuboid.
   * \param maxY Maximum y-coordinate (in image space) of cuboid.
   * \param maxZ Maximum z-coordinate (in image space) of cuboid.
   * \return Intensity image.
   **/
  float evaluate_integral_image(const Image& integral, int minX, int minY, int minZ, int maxX, int maxY, int maxZ);

  /**
   * \brief Computes an integral histogram for given intensity image.
   * \param image The image.
   * \param bins Number of histogram bins.
   * \param minValue Minimum intensity value considered for the histogram counts.
   * \param maxValue Maximum intensity value considered for the histogram counts.
   * \return The integral histogram.
   **/
  IntegralHistogram integral_histogram(const Image& image, int bins, float minValue, float maxValue);

  /**
   * \brief Alternative method for Computing an integral histogram for given intensity image. Uses temporary images.
   * \param image The image.
   * \param bins Number of histogram bins.
   * \param minValue Minimum intensity value considered for the histogram counts.
   * \param maxValue Maximum intensity value considered for the histogram counts.
   * \return The integral histogram.
   **/
  IntegralHistogram integral_histogram_alternative(const Image& image, int bins, float minValue, float maxValue);

  /**
   * \brief Reconstructs an intensity image from a given integral histogram. Uses the bin index with maximum count as reconstructed value.
   * \param integral The integral image.
   * \param image Reconstructed intensity image.
   **/
  void reverse_integral_histogram(const IntegralHistogram& integral, Image& image);

  /**
   * \brief Computes the histogram for a given cuboid using an integral histogram.
   * \param integral The integral histogram.
   * \param minX Minimum x-coordinate (in image space) of cuboid.
   * \param minY Minimum y-coordinate (in image space) of cuboid.
   * \param minZ Minimum z-coordinate (in image space) of cuboid.
   * \param maxX Maximum x-coordinate (in image space) of cuboid.
   * \param maxY Maximum y-coordinate (in image space) of cuboid.
   * \param maxZ Maximum z-coordinate (in image space) of cuboid.
   * \return Intensity histogram.
   **/
  Histogram evaluate_integral_histogram(const IntegralHistogram& integral, int minX, int minY, int minZ, int maxX, int maxY, int maxZ);

  /**
   * \brief Thresholds an intensity image and returns the binary mask. Points in the binary mask are set according to B(x) = I(x) >= minValue && I(x) <= maxValue.
   * \param image The image.
   * \param binary The binary mask.
   * \param minValue The minimum intensity value to be included.
   * \param maxValue The maximum intensity value to be included.
   **/
  void threshold(const Image& image, Image& binary, float minValue, float maxValue);

  /**
   * \brief Generates a random binary mask given a specified probability for an image point to be set to one.
   * \param image The binary image.
   * \param probability The probability for an image point to be set to one.
   **/
  void random_binary(Image& image, float probability);

  /**
   * \brief Inverts a binary image.
   * \param input The input binary image.
   * \param output The inverted binary image.
   **/
  void invert_binary(const Image& input, Image& output);

  /**
   * \brief Dilates a binary image with 3x3x3 box filter.
   * \param input The input binary image.
   * \param output The dilated binary image.
   **/
  void dilate_binary(const Image& input, Image& output);

  /**
   * \brief Erodes a binary image with 3x3x3 box filter.
   * \param input The input binary image.
   * \param output The eroded binary image.
   **/
  void erode_binary(const Image& input, Image& output);

  /**
   * \brief Computes an Euclidean distance map for a given binary image.
   * \param binary The binary image.
   * \param distmap The Euclidean distance map.
   * \param numPasses The number of forward/backward passes.
   **/
  void euclidean_distmap(const Image& binary, Image& distmap, int numPasses);

  /**
   * \brief More efficient computation of an Euclidean distance map for a given binary image. Uses temporary padded images to avoid boundary checks.
   * \param binary The binary image.
   * \param distmap The Euclidean distance map.
   * \param numPasses The number of forward/backward passes.
   **/
  void euclidean_distmap_fast(const Image& binary, Image& distmap, int numPasses);

  /**
   * \brief Computes a geodesic distance map for a given intensity and binary image.
   * \param image The intensity image.
   * \param binary The binary image.
   * \param distmap The geodesic distance map.
   * \param intensityWeight The weighting for the intensity component.
   * \param numPasses The number of forward/backward passes.
   **/
  void geodesic_distmap(const Image& image, const Image& binary, Image& distmap, float intensityWeight, int numPasses);

  /**
   * \brief More efficient computation of a geodesic distance map for a intensity and given binary image. Uses temporary padded images to avoid boundary checks.
   * \param image The intensity image.
   * \param binary The binary image.
   * \param distmap The geodesic distance map.
   * \param intensityWeight The weighting for the intensity component.
   * \param numPasses The number of forward/backward passes.
   **/
  void geodesic_distmap_fast(const Image& image, const Image& binary, Image& distmap, float intensityWeight, int numPasses);

  /**
   * \brief Extracts a 2-D slice from a 3-D volume. Copies image data. Slice does not have correct image to world transformations.
   * \param image The 3-D image.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   * \return The 2-D slice.
   **/
  Image extract_slice(const Image& image, ImagePlane plane, int sliceNumber);

  /**
   * \brief Copies the values of a 2-D slice into a 3-D volume at the specified position.
   * \param slice The 2-D slice.
   * \param image The 3-D image.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   **/
  void insert_slice(const Image& slice, Image& image, ImagePlane plane, int sliceNumber);

  /**
   * \brief Generates a 8-bit grayscale 2D texture for image rendering. Intensity values are mapped to [0,255].
   * \param image The image.
   * \param texture Pointer to the texture buffer.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   * \param window The intensity window.
   * \param level The intensity level.
   **/
  void texture_8bit_gray(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level);

  /**
   * \brief Generates a 32-bit (RGBA) grayscale 2-D texture for image rendering. Intensity values are mapped to [0,255].
   * \param image The image.
   * \param texture Pointer to the texture buffer.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   * \param window The intensity window.
   * \param level The intensity level.
   **/
  void texture_32bit_gray(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level);

  /**
   * \brief Generates a 32-bit (RGBA) pseudo-color 2-D texture for image rendering. Uses the JET64 colormap. Intensity values are mapped to [0,255].
   * \param image The image.
   * \param texture Pointer to the texture buffer.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   * \param window The intensity window.
   * \param level The intensity level.
   * \param alpha The alpha value used for blending.
   **/
  void texture_32bit_color_jet64(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level, float alpha);

  /**
   * \brief Generates a 32-bit (RGBA) pseudo-color 2-D texture for image rendering. Uses the INDEXED32 colormap Intensity values are mapped to [0,255].
   * \param image The image.
   * \param texture Pointer to the texture buffer.
   * \param plane The plane orientation.
   * \param sliceNumber The slice number.
   * \param window The intensity window.
   * \param level The intensity level.
   * \param alpha The alpha value used for blending.
   **/
  void texture_32bit_color_indexed32(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float alpha);

  /**
   * \brief Generates a binary checkerboard image.
   * \param image The image.
   * \param patchSizeX Patch size (in pixels) along the x-axis.
   * \param patchSizeY Patch size (in pixels) along the y-axis.
   * \param patchSizeZ Patch size (in pixels) along the z-axis.
   **/
  void checkerboard(Image& image, int patchSizeX, int patchSizeY, int patchSizeZ);

  /**
   * \brief Generates a checkerboard image for two given intensity images. The images must have equal size.
   * \param inputA The input image A.
   * \param inputB The input image B.
   * \param patchSizeX Patch size (in pixels) along the x-axis.
   * \param patchSizeY Patch size (in pixels) along the y-axis.
   * \param patchSizeZ Patch size (in pixels) along the z-axis.
   * \return The checkerboard image.
   **/
  Image checkerboard(const Image& inputA, const Image& inputB, int patchSizeX, int patchSizeY, int patchSizeZ);

  namespace impl
  {
    /**
     * \brief Internal implementation of image transformation. Uses backward warping.
     * \param input The image.
     * \param output The warped output image.
     * \param transformInImageSpace The affine transformation matrix mapping from image to image space.
     * \param interpolation The interpolation mode.
     * \param outsideValue The value to be used for out of bounds image points.
     **/
    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, Interpolation interpolation, float outsideValue);

    /**
     * \brief Internal implementation of image transformation including a discplacement field. Uses backward warping.
     * \param input The image.
     * \param output The warped output image.
     * \param transformInImageSpace The affine transformation matrix mapping from image to image space.
     * \param fieldX The x-component of the displacement field.
     * \param fieldY The y-component of the displacement field.
     * \param fieldZ The z-component of the displacement field.
     * \param interpolation The interpolation mode.
     * \param outsideValue The value to be used for out of bounds image points.
     **/
    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation, float outsideValue);

    /**
    * \brief Generates a 1-D Gauss filter kernel with a given sigma.
    * \param sigmaInPixels Sigma in pixels.
    * \return The filter kernel.
    **/
    std::vector<double> generate_gauss_filter(double sigmaInPixels);

    /**
    * \brief Forward propagation pass for Euclidean distance map computation.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    **/
    void forward_pass(Image& distmap, const Image& point_distances);

    /**
    * \brief Backward propagation pass for Euclidean distance map computation.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    **/
    void backward_pass(Image& distmap, const Image& point_distances);

    /**
    * \brief Fast forward propagation pass for Euclidean distance map computation. No boundary checks.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    **/
    void forward_pass_fast(Image& distmap, const Image& point_distances);

    /**
    * \brief Fast backward propagation pass for Euclidean distance map computation. No boundary checks.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    **/
    void backward_pass_fast(Image& distmap, const Image& point_distances);

    /**
    * \brief Forward propagation pass for geodesic distance map computation.
    * \param image The intensity image.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    * \param intensityWeight The weighting of the intensity component.
    **/
    void forward_pass(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight);

    /**
    * \brief Backward propagation pass for geodesic distance map computation.
    * \param image The intensity image.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    * \param intensityWeight The weighting of the intensity component.
    **/
    void backward_pass(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight);

    /**
    * \brief Fast forward propagation pass for geodesic distance map computation. No boundary checks.
    * \param image The intensity image.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    * \param intensityWeight The weighting of the intensity component.
    **/
    void forward_pass_fast(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight);

    /**
    * \brief Fast backward propagation pass for geodesic distance map computation. No boundary checks.
    * \param image The intensity image.
    * \param distmap The distance map.
    * \param point_distances A 3x3x3 image containing the 26-neighbourhood distances from the center point in element spacing units.
    * \param intensityWeight The weighting of the intensity component.
    **/
    void backward_pass_fast(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight);
  }

  namespace colormaps
  {
    /**
    * \brief JET64 RGBA colormap.
    **/
    static double jet64[][4] =
    {
      { 0.0000000e+000,  0.0000000e+000,  5.6250000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  6.2500000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  6.8750000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  7.5000000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  8.1250000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  8.7500000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  9.3750000e-001,   1.0000},
      { 0.0000000e+000,  0.0000000e+000,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  6.2500000e-002,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  1.2500000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  1.8750000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  2.5000000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  3.1250000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  3.7500000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  4.3750000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  5.0000000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  5.6250000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  6.2500000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  6.8750000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  7.5000000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  8.1250000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  8.7500000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  9.3750000e-001,  1.0000000e+000,   1.0000},
      { 0.0000000e+000,  1.0000000e+000,  1.0000000e+000,   1.0000},
      { 6.2500000e-002,  1.0000000e+000,  9.3750000e-001,   1.0000},
      { 1.2500000e-001,  1.0000000e+000,  8.7500000e-001,   1.0000},
      { 1.8750000e-001,  1.0000000e+000,  8.1250000e-001,   1.0000},
      { 2.5000000e-001,  1.0000000e+000,  7.5000000e-001,   1.0000},
      { 3.1250000e-001,  1.0000000e+000,  6.8750000e-001,   1.0000},
      { 3.7500000e-001,  1.0000000e+000,  6.2500000e-001,   1.0000},
      { 4.3750000e-001,  1.0000000e+000,  5.6250000e-001,   1.0000},
      { 5.0000000e-001,  1.0000000e+000,  5.0000000e-001,   1.0000},
      { 5.6250000e-001,  1.0000000e+000,  4.3750000e-001,   1.0000},
      { 6.2500000e-001,  1.0000000e+000,  3.7500000e-001,   1.0000},
      { 6.8750000e-001,  1.0000000e+000,  3.1250000e-001,   1.0000},
      { 7.5000000e-001,  1.0000000e+000,  2.5000000e-001,   1.0000},
      { 8.1250000e-001,  1.0000000e+000,  1.8750000e-001,   1.0000},
      { 8.7500000e-001,  1.0000000e+000,  1.2500000e-001,   1.0000},
      { 9.3750000e-001,  1.0000000e+000,  6.2500000e-002,   1.0000},
      { 1.0000000e+000,  1.0000000e+000,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  9.3750000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  8.7500000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  8.1250000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  7.5000000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  6.8750000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  6.2500000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  5.6250000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  5.0000000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  4.3750000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  3.7500000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  3.1250000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  2.5000000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  1.8750000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  1.2500000e-001,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  6.2500000e-002,  0.0000000e+000,   1.0000},
      { 1.0000000e+000,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 9.3750000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 8.7500000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 8.1250000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 7.5000000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 6.8750000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 6.2500000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 5.6250000e-001,  0.0000000e+000,  0.0000000e+000,   1.0000},
      { 5.0000000e-001,  0.0000000e+000,  0.0000000e+000,   0.0000}
    };

    /**
    * \brief INDEXED32 RGBA colormap.
    **/
    static double indexed32[][4] =
    {
      {       0,        0,        0,        0 }, //BLACK
      {       0,        0,   1.0000,   1.0000 },
      {       0,   1.0000,        0,   1.0000 },
      {  1.0000,        0,        0,   1.0000 },
      {  1.0000,        0,   0.7241,   1.0000 },
      {  1.0000,   0.8276,        0,   1.0000 },
      {       0,   0.4483,   0.9310,   1.0000 },
      {       0,   1.0000,   0.7586,   1.0000 },
      {  0.3103,   0.5172,        0,   1.0000 },
      {  0.6207,   0.0345,   0.2069,   1.0000 },
      {  0.3103,        0,   0.3793,   1.0000 },
      {  1.0000,   0.6207,   0.4138,   1.0000 },
      {  0.9655,   0.6552,   1.0000,   1.0000 },
      {       0,   0.7586,   1.0000,   1.0000 },
      {  0.6897,   0.3103,   0.9655,   1.0000 },
      {  0.8621,   1.0000,   0.6207,   1.0000 },
      {  0.4483,   0.3103,        0,   1.0000 },
      {  0.7931,   1.0000,        0,   1.0000 },
      {       0,   0.5517,   0.4138,   1.0000 },
      {       0,        0,   0.6207,   1.0000 },
      {  1.0000,   0.5517,   0.6552,   1.0000 },
      {       0,   0.2414,   0.4828,   1.0000 },
      {  0.3448,   0.8621,   0.3448,   1.0000 },
      {  0.5172,   1.0000,   1.0000,   1.0000 },
      {       0,   0.2414,        0,   1.0000 },
      {  0.6897,   0.2414,        0,   1.0000 },
      {  1.0000,        0,   0.3103,   1.0000 },
      {  0.2759,        0,   0.1379,   1.0000 },
      {  0.6207,   0.3448,   0.5517,   1.0000 },
      {  0.5517,   0.5862,   0.8276,   1.0000 },
      {  0.7241,   0.6552,   0.2414,   1.0000 },
      {  1.0000,        0,   1.0000,   1.0000 },
    };
  }
}
