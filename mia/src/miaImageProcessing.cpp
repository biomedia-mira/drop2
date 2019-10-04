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

#include "miaImageProcessing.h"
#include "miaImageIterator.h"

#include <limits>
#include <random>



namespace mia
{
  const double MIA_PI = 3.141592653589793;

  void zeros(Image& image)
  {
    fill(image, 0.0f);
  }

  void ones(Image& image)
  {
    fill(image, 1.0f);
  }

  void fill(Image& image, float value)
  {
    for (auto& p : image)
      p = value;
  }

  void copy(const Image& src, Image& dst)
  {
    if (src.contiguous() && dst.contiguous())
    {
      const float* src_start = &src(0,0,0);
      float* dst_start = &dst(0,0,0);
      std::copy(src_start, src_start+src.size(), dst_start);
    }
    else
    {
      if (src.stepX() == 1 && dst.stepX() == 1)
      {
        for (int z = 0; z < src.sizeZ(); z++)
        {
          for (int y = 0; y < src.sizeY(); y++)
          {
            const auto* src_start = &src(0,y,z);
            auto* dst_start = &dst(0,y,z);
            std::copy(src_start, src_start+src.sizeX(), dst_start);
          }
        }
      }
      else
      {
        for (int z = 0; z < src.sizeZ(); z++)
        {
          for (int y = 0; y < src.sizeY(); y++)
          {
            for (int x = 0; x < src.sizeX(); x++)
            {
              dst(x,y,z) = src(x,y,z);
            }
          }
        }
      }
    }
  }

  void copy_region(const Image& src, Image& dst, int srcX, int srcY, int srcZ, int dstX, int dstY, int dstZ, int regionSizeX, int regionSizeY, int regionSizeZ)
  {
    for (int z = 0; z < regionSizeZ; z++)
    {
      for (int y = 0; y < regionSizeY; y++)
      {
        for (int x = 0; x < regionSizeX; x++)
        {
          dst(dstX+x,dstY+y,dstZ+z) = src(srcX+x,srcY+y,srcZ+z);
        }
      }
    }
  }

  float min(const Image& image)
  {
    float value = std::numeric_limits<float>::max();
    for (auto& p : image)
    {
      if (p < value) value = p;
    }
    return value;
  }

  float max(const Image& image)
  {
    float value = std::numeric_limits<float>::lowest();
    for (auto& p : image)
    {
      if (p > value) value = p;
    }
    return value;
  }

  float mean(const Image& image)
  {
    float value = 0.0f;
    for (auto& p : image)
    {
      value += p;
    }
    return value / static_cast<float>(image.size());
  }

  float var(const Image& image)
  {
    float sum = 0.0f;
    float sumsq = 0.0f;
    for (auto& p : image)
    {
      sum += p;
      sumsq += p * p;
    }
    float mean = sum / static_cast<float>(image.size());

    return sumsq / static_cast<float>(image.size()) - mean;
  }

  void add(Image& image, float value)
  {
    for (auto& p : image)
      p += value;
  }

  void sub(Image& image, float value)
  {
    for (auto& p : image)
      p -= value;
  }

  void mul(Image& image, float value)
  {
    for (auto& p : image)
      p *= value;
  }

  void div(Image& image, float value)
  {
    for (auto& p : image)
      p /= value;
  }

  void add(const Image& imageA, const Image& imageB, Image& imageC)
  {
    for (int z = 0; z < imageA.sizeZ(); z++)
    {
      for (int y = 0; y < imageA.sizeY(); y++)
      {
        for (int x = 0; x < imageA.sizeX(); x++)
        {
          imageC(x,y,z) = imageA(x,y,z) + imageB(x,y,z);
        }
      }
    }
  }

  void sub(const Image& imageA, const Image& imageB, Image& imageC)
  {
    for (int z = 0; z < imageA.sizeZ(); z++)
    {
      for (int y = 0; y < imageA.sizeY(); y++)
      {
        for (int x = 0; x < imageA.sizeX(); x++)
        {
          imageC(x,y,z) = imageA(x,y,z) - imageB(x,y,z);
        }
      }
    }
  }

  void mul(const Image& imageA, const Image& imageB, Image& imageC)
  {
    for (int z = 0; z < imageA.sizeZ(); z++)
    {
      for (int y = 0; y < imageA.sizeY(); y++)
      {
        for (int x = 0; x < imageA.sizeX(); x++)
        {
          imageC(x,y,z) = imageA(x,y,z) * imageB(x,y,z);
        }
      }
    }
  }

  void div(const Image& imageA, const Image& imageB, Image& imageC)
  {
    for (int z = 0; z < imageA.sizeZ(); z++)
    {
      for (int y = 0; y < imageA.sizeY(); y++)
      {
        for (int x = 0; x < imageA.sizeX(); x++)
        {
          imageC(x,y,z) = imageA(x,y,z) / imageB(x,y,z);
        }
      }
    }
  }

  Histogram histogram(const Image& image, int numBins, float minValue, float maxValue)
  {
    Histogram h(numBins, minValue, maxValue);
    for (auto& p : image)
      h.add(p);
    return h;
  }

  Image subimage(Image& image, int x, int y, int z, int sizeX, int sizeY, int sizeZ)
  {
    float* data = image.data() + x*image.stepX() + y*image.stepY() + z*image.stepZ();
    Image subimage(data, sizeX, sizeY, sizeZ, image.stepX(), image.stepY(), image.stepZ());

    // set image datatype
    subimage.dataType(image.dataType());

    // set element spacing
    subimage.spacing(image.spacing());

    // set image-to-world orientation
    subimage.imageToWorld(image.imageToWorld());

    // adjust image origin
    Eigen::Vector3d offset(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z));
    subimage.origin(image.origin() + (image.imageToWorld() * offset.cwiseProduct(image.spacing())));

    return subimage;
  }

  Image resample(const Image& image, const Eigen::Vector3d& spacing, Interpolation interpolation, float outsideValue, int minimumSize)
  {
    Eigen::Vector3d spacing_output = spacing;

    int sizeX = static_cast<int>(image.sizeX() * image.spacing()[0] / spacing[0]);
    int sizeY = static_cast<int>(image.sizeY() * image.spacing()[1] / spacing[1]);
    int sizeZ = static_cast<int>(image.sizeZ() * image.spacing()[2] / spacing[2]);

    if (sizeX < minimumSize)
    {
      sizeX = minimumSize;
      spacing_output[0] = static_cast<double>(image.sizeX()) / static_cast<double>(minimumSize) * image.spacing()[0];
    }
    if (sizeY < minimumSize)
    {
      sizeY = minimumSize;
      spacing_output[1] = static_cast<double>(image.sizeY()) / static_cast<double>(minimumSize) * image.spacing()[1];
    }
    if (sizeZ < minimumSize)
    {
      sizeZ = minimumSize;
      spacing_output[2] = static_cast<double>(image.sizeZ()) / static_cast<double>(minimumSize) * image.spacing()[2];
    }

    Image output(sizeX, sizeY, sizeZ);

    // set image datatype
    output.dataType(image.dataType());

    // set element spacing
    output.spacing(spacing_output);

    // set image-to-world orientation
    output.imageToWorld(image.imageToWorld());

    // set image origin
    output.origin(image.origin());
    output.origin(image.origin() + (image.imageCenterToWorld() - output.imageCenterToWorld()));

    resample(image, output, interpolation, outsideValue);

    return output;
  }
  void resample(const Image& input, Image& output, Interpolation interpolation, float outsideValue)
  {
    Eigen::Matrix4d transformInImageSpace = input.worldToImageTransform() * output.imageToWorldTransform();
    impl::transform(input, output, transformInImageSpace, interpolation, outsideValue);
  }

  void warp(const Image& input, Image& output, const Eigen::Matrix4d& transform, Interpolation interpolation, float outsideValue)
  {
    Eigen::Matrix4d transformInImageSpace = input.worldToImageTransform() * transform * output.imageToWorldTransform();
    impl::transform(input, output, transformInImageSpace, interpolation, outsideValue);
  }
  void warp(const Image& input, Image& output, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation, float outsideValue)
  {
    Eigen::Matrix4d transformInImageSpace = input.worldToImageTransform() * output.imageToWorldTransform();
    impl::transform(input, output, transformInImageSpace, fieldX, fieldY, fieldZ, interpolation, outsideValue);
  }
  void warp(const Image& input, Image& output, const Eigen::Matrix4d& transform, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation, float outsideValue)
  {
    Eigen::Matrix4d transformInImageSpace = input.worldToImageTransform() * transform * output.imageToWorldTransform();
    impl::transform(input, output, transformInImageSpace, fieldX, fieldY, fieldZ, interpolation, outsideValue);
  }

  void compose(const Image& aX, const Image& aY, const Image& aZ, Image& bX, Image& bY, Image& bZ)
  {
    int sizeX = bX.sizeX();
    int sizeY = bX.sizeY();
    int sizeZ = bX.sizeZ();

    Eigen::Vector3d spacingB = bX.spacing();

    for (int z = 0; z < sizeZ; z++)
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          auto dx = x + bX(x, y, z) / spacingB[0];
          auto dy = y + bY(x, y, z) / spacingB[1];
          auto dz = z + bZ(x, y, z) / spacingB[2];

          if (dx < 0) dx = 0;
          else if (dx > sizeX - 1) dx = sizeX - 1;
          if (dy < 0) dy = 0;
          else if (dy > sizeY - 1) dy = sizeY - 1;
          if (dz < 0) dz = 0;
          else if (dz > sizeZ - 1) dz = sizeZ - 1;

          Eigen::Vector3d d(dx, dy, dz);

          bX(x, y, z) += aX.linear(d);
          bY(x, y, z) += aY.linear(d);
          bZ(x, y, z) += aZ.linear(d);
        }
      }
    }
  }

  void compose(const Eigen::Matrix4d& transform, Image& fieldX, Image& fieldY, Image& fieldZ)
  {
    for (int z = 0; z < fieldX.sizeZ(); z++)
    {
      for (int y = 0; y < fieldX.sizeY(); y++)
      {
        for (int x = 0; x < fieldX.sizeX(); x++)
        {
          Eigen::Vector3d s = fieldX.mapImageToWorld(Eigen::Vector3d(x, y, z));
          Eigen::Vector4d t = transform * Eigen::Vector4d(s[0] + fieldX(x,y,z), s[1] + fieldY(x, y, z), s[2] + fieldZ(x, y, z), 1);
          Eigen::Vector3d d = t.segment(0, 3) - s;
          fieldX(x, y, z) = static_cast<float>(d[0]);
          fieldY(x, y, z) = static_cast<float>(d[1]);
          fieldZ(x, y, z) = static_cast<float>(d[2]);
        }
      }
    }
  }

  void gauss(const Image& input, Image& output, double sigmaX, double sigmaY, double sigmaZ)
  {
    double sigmaXinPixels = sigmaX / input.spacing()[0];
    double sigmaYinPixels = sigmaY / input.spacing()[1];
    double sigmaZinPixels = sigmaZ / input.spacing()[2];

    auto filterX = impl::generate_gauss_filter(sigmaXinPixels);
    auto filterY = impl::generate_gauss_filter(sigmaYinPixels);
    auto filterZ = impl::generate_gauss_filter(sigmaZinPixels);

    if (static_cast<int>(filterX.size()) > input.sizeX() * 2)
    {
      filterX.resize(1);
      filterX[0] = 1.0;
    }
    if (static_cast<int>(filterY.size()) > input.sizeY() * 2)
    {
      filterY.resize(1);
      filterY[0] = 1.0;
    }
    if (static_cast<int>(filterZ.size()) > input.sizeZ() * 2)
    {
      filterZ.resize(1);
      filterZ[0] = 1.0;
    }

    auto interm = input.clone();

    filter_x(interm, output, filterX);
    filter_y(output, interm, filterY);
    filter_z(interm, output, filterZ);
  }

  void median(const Image& input, Image& output, int filterRadiusX, int filterRadiusY, int filterRadiusZ)
  {
    copy(input, output);

    int numFilterElements = (filterRadiusX*2+1)*(filterRadiusY*2+1)*(filterRadiusZ*2+1);

    for (int z = filterRadiusZ; z < input.sizeZ()-filterRadiusZ; z++)
    {
      for (int y = filterRadiusY; y < input.sizeY()-filterRadiusY; y++)
      {
        for (int x = filterRadiusX; x < input.sizeX()-filterRadiusX; x++)
        {
          std::vector<float> values;
          for (int dz = -filterRadiusZ; dz <= filterRadiusZ; dz++)
          {
            for (int dy = -filterRadiusY; dy <= filterRadiusY; dy++)
            {
              for (int dx = -filterRadiusX; dx <= filterRadiusX; dx++)
              {
                values.push_back(input(x+dx,y+dy,z+dz));
              }
            }
          }
          std::sort(values.begin(), values.end());
          output(x,y,z) = values[(numFilterElements-1)/2];
        }
      }
    }
  }

  void filter_x(const Image& input, Image& output, const std::vector<double>& filter)
  {
    int filterLength = static_cast<int>(filter.size());
    int filterOffset = (filterLength - 1) / 2;

    //inner part
    for (int z = 0; z < input.sizeZ(); z++)
    {
      for (int y = 0; y < input.sizeY(); y++)
      {
        for (int x = filterOffset; x < input.sizeX() - filterOffset; x++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            filterResponse += static_cast<float>(input(x+k-filterOffset,y,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }

    // borders
    int limit1 = filterOffset;
    if (limit1 > input.sizeX() - 1) limit1 = input.sizeX() - 1;
    int limit2 = input.sizeX() - filterOffset;
    if (limit2 < 0) limit2 = 0;
    for (int z = 0; z < input.sizeZ(); z++)
    {
      for (int y = 0; y < input.sizeY(); y++)
      {
        for (int x = 0; x < limit1; x++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = x+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeX()) inputOffset = 2*input.sizeX()-inputOffset-2;
            filterResponse += static_cast<float>(input(inputOffset,y,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
        for (int x = limit2; x < input.sizeX(); x++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = x+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeX()) inputOffset = 2*input.sizeX()-inputOffset-2;
            filterResponse += static_cast<float>(input(inputOffset,y,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }
  }
  void filter_y(const Image& input, Image& output, const std::vector<double>& filter)
  {
    int filterLength = static_cast<int>(filter.size());
    int filterOffset = (filterLength - 1) / 2;

    //inner part
    for (int z = 0; z < input.sizeZ(); z++)
    {
      for (int x = 0; x < input.sizeX(); x++)
      {
        for (int y = filterOffset; y < input.sizeY() - filterOffset; y++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            filterResponse += static_cast<float>(input(x,y+k-filterOffset,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }

    // borders
    int limit1 = filterOffset;
    if (limit1 > input.sizeY() - 1) limit1 = input.sizeY() - 1;
    int limit2 = input.sizeY() - filterOffset;
    if (limit2 < 0) limit2 = 0;
    for (int z = 0; z < input.sizeZ(); z++)
    {
      for (int x = 0; x < input.sizeX(); x++)
      {
        for (int y = 0; y < limit1; y++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = y+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeY()) inputOffset = 2*input.sizeY()-inputOffset-2;
            filterResponse += static_cast<float>(input(x,inputOffset,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
        for (int y = limit2; y < input.sizeY(); y++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = y+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeY()) inputOffset = 2*input.sizeY()-inputOffset-2;
            filterResponse += static_cast<float>(input(x,inputOffset,z) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }
  }
  void filter_z(const Image& input, Image& output, const std::vector<double>& filter)
  {
    int filterLength = static_cast<int>(filter.size());
    int filterOffset = (filterLength - 1) / 2;

    //inner part
    for (int y = 0; y < input.sizeY(); y++)
    {
      for (int x = 0; x < input.sizeX(); x++)
      {
        for (int z = filterOffset; z < input.sizeZ() - filterOffset; z++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            filterResponse += static_cast<float>(input(x,y,z+k-filterOffset) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }

    // borders
    int limit1 = filterOffset;
    if (limit1 > input.sizeZ() - 1) limit1 = input.sizeZ() - 1;
    int limit2 = input.sizeZ() - filterOffset;
    if (limit2 < 0) limit2 = 0;
    for (int y = 0; y < input.sizeY(); y++)
    {
      for (int x = 0; x < input.sizeX(); x++)
      {
        for (int z = 0; z < limit1; z++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = z+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeZ()) inputOffset = 2*input.sizeZ()-inputOffset-2;
            filterResponse += static_cast<float>(input(x,y,inputOffset) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
        for (int z = limit2; z < input.sizeZ(); z++)
        {
          float filterResponse = 0.0f;
          for (int k = 0; k < filterLength; k++)
          {
            int inputOffset = z+k-filterOffset;
            if (inputOffset < 0) inputOffset = -inputOffset;
            else if (inputOffset >= input.sizeZ()) inputOffset = 2*input.sizeZ()-inputOffset-2;
            filterResponse += static_cast<float>(input(x,y,inputOffset) * filter[k]);
          }
          output(x,y,z) = filterResponse;
        }
      }
    }
  }

  void fill_boundaries(Image& image, float value)
  {
    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    for (int z = 0; z < sizeZ; z++)
    {
      for (int x = 0; x < sizeX; x++)
      {
        image(x, 0, z) = value;
        image(x, sizeY - 1, z) = value;
      }
    }
    for (int z = 0; z < sizeZ; z++)
    {
      for (int y = 1; y < sizeY - 1; y++)
      {
        image(0, y, z) = value;
        image(sizeX - 1, y, z) = value;
      }
    }
    if (sizeZ > 2)
    {
      for (int y = 1; y < sizeY - 1; y++)
      {
        for (int x = 1; x < sizeX - 1; x++)
        {
          image(x, y, 0) = value;
          image(x, y, sizeZ - 1) = value;
        }
      }
    }
  }

  bool same_domain(const Image& imageA, const Image& imageB)
  {
    bool sameSize = imageA.sizeX() == imageB.sizeX() && imageA.sizeY() == imageB.sizeY() && imageA.sizeZ() == imageB.sizeZ();
    bool sameSpacing = imageA.spacing() == imageB.spacing();
    bool sameOrigin = imageA.origin() == imageB.origin();
    bool sameOrientation = imageA.imageToWorld() == imageB.imageToWorld();

    return sameSize && sameSpacing && sameOrigin && sameOrientation;
  }

  std::vector<bool> flip_flags(const Image& image)
  {
    std::vector<bool> flags(3);

    flags[0] = image.imageToWorld()(0,0) < 0;
    flags[1] = image.imageToWorld()(1,1) < 0;
    flags[2] = image.imageToWorld()(2,2) < 0;

    return flags;
  }

  Image reorient(const Image& image)
  {
    std::vector<double> max_values(3);
    std::vector<int> max_value_indices(3);
    std::fill(max_values.begin(), max_values.end(), 0.0);
    std::fill(max_value_indices.begin(), max_value_indices.end(), 0);

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        auto value = abs(image.imageToWorld()(i,j));
        if (value > max_values[i])
        {
          max_values[i] = value;
          max_value_indices[i] = j;
        }
      }
    }

    std::vector<int> old_size(3);
    old_size[0] = image.sizeX();
    old_size[1] = image.sizeY();
    old_size[2] = image.sizeZ();

    std::vector<int> new_size(3);
    Eigen::Vector3d new_spacing;
    Eigen::Vector3d new_extent;
    Eigen::Matrix3d new_imageToWorld;

    for (int i = 0; i < 3; i++)
    {
      new_size[i] = old_size[max_value_indices[i]];
      new_spacing[i] = image.spacing()[max_value_indices[i]];
      new_extent[i] = static_cast<double>(new_size[i]) * new_spacing[i];
      for (int j = 0; j < 3; j++)
      {
        new_imageToWorld(i,j) = image.imageToWorld()(i, max_value_indices[j]);
      }
    }

    Eigen::Vector3d flip_vector;
    Eigen::Matrix3d flip_matrix = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; i++)
    {
      bool flip = new_imageToWorld(i,i) < 0;
      flip_vector[i] = flip ? 1.0 : 0.0;
      flip_matrix(i, i) = flip ? -1.0 : 1.0;
    }

    Image reoriented(new_size[0], new_size[1], new_size[2]);
    reoriented.spacing(new_spacing);
    reoriented.origin(image.origin() + new_imageToWorld * new_extent.cwiseProduct(flip_vector));
    reoriented.imageToWorld(new_imageToWorld * flip_matrix);
    reoriented.dataType(image.dataType());

    resample(image, reoriented);

    return reoriented;
  }

  Eigen::Vector3d center_of_mass(const Image& image)
  {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();

    float sum_weights = 0.0f;
    float min_value = min(image);
    for (int z = 0; z < image.sizeZ(); z++)
    {
      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          float weight = image(x,y,z) - min_value;
          center += image.mapImageToWorld(Eigen::Vector3d(x,y,z)) * weight;
          sum_weights += weight;
        }
      }
    }

    if (sum_weights > 0.0f) center /= sum_weights;

    return center;
  }

  Image pad(const Image& image, int padX1, int padX2, int padY1, int padY2, int padZ1, int padZ2, float value)
  {
    Image padimage(image.sizeX()+padX1+padX2,image.sizeY()+padY1+padY2,image.sizeZ()+padZ1+padZ2);

    // set image datatype
    padimage.dataType(image.dataType());

    // set element spacing
    padimage.spacing(image.spacing());

    // set image-to-world orientation
    padimage.imageToWorld(image.imageToWorld());

    // adjust image origin
    Eigen::Vector3d offset(static_cast<double>(padX1), static_cast<double>(padY1), static_cast<double>(padZ1));
    padimage.origin(image.origin() - (image.imageToWorld() * offset.cwiseProduct(image.spacing())));

    fill(padimage, value);

    copy_region(image, padimage, 0, 0, 0, padX1, padY1, padZ1, image.sizeX(), image.sizeY(), image.sizeZ());

    return padimage;
  }

  Image crop(const Image& image, int cropX1, int cropX2, int cropY1, int cropY2, int cropZ1, int cropZ2)
  {
    return subimage(const_cast<Image&>(image), cropX1, cropY1, cropZ1, image.sizeX()-cropX1-cropX2, image.sizeY()-cropY1-cropY2, image.sizeZ()-cropZ1-cropZ2).clone();
  }

  Image integral_image(const Image& image)
  {
    Image integral = pad(image,1,0,1,0,1,0);
    integral.dataType(ImageDataType::FLOAT);

    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    for (int x = 1; x <= sizeX; x++)
    {
      integral(x,0,0) = image(x-1,0,0) + integral(x-1,0,0);
    }

    for (int y = 1; y <= sizeY; y++)
    {
      integral(0,y,0) = image(0,y-1,0) + integral(0,y-1,0);
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      integral(0,0,z) = image(0,0,z-1) + integral(0,0,z-1);
    }

    for (int y = 1; y <= sizeY; y++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        integral(x,y,0) = image(x-1,y-1,0) + integral(x-1,y,0) + integral(x,y-1,0) - integral(x-1,y-1,0);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        integral(x,0,z) = image(x-1,0,z-1) + integral(x-1,0,z) + integral(x,0,z-1) - integral(x-1,0,z-1);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        integral(0,y,z) = image(0,y-1,z-1) + integral(0,y-1,z) + integral(0,y,z-1) - integral(0,y-1,z-1);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        for (int x = 1; x <= sizeX; x++)
        {
          integral(x,y,z) = image(x-1,y-1,z-1) + integral(x,y,z-1) + integral(x,y-1,z) - integral(x,y-1,z-1) + integral(x-1,y,z) - integral(x-1,y,z-1) - integral(x-1,y-1,z) + integral(x-1,y-1,z-1);
        }
      }
    }

    return integral;
  }

  Image integral_image_alternative(const Image& image)
  {
    Image integral = pad(image,1,0,1,0,1,0);
    Image temp1 = integral.clone();
    Image temp2 = integral.clone();

    integral.dataType(ImageDataType::FLOAT);

    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        for (int x = 1; x <= sizeX; x++)
        {
          temp1(x,y,z) = temp1(x,y,z-1) + image(x-1,y-1,z-1);
          temp2(x,y,z) = temp2(x,y-1,z) + temp1(x,y,z);
          integral(x,y,z) = integral(x-1,y,z) + temp2(x,y,z);
        }
      }
    }

    return integral;
  }

  void reverse_integral_image(const Image& integral, Image& image)
  {
    for (int z = 0; z < image.sizeZ(); z++)
    {
      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          image(x,y,z) = evaluate_integral_image(integral, x, y, z, x+1, y+1, z+1);
        }
      }
    }
  }

  float evaluate_integral_image(const Image& integral, int minX, int minY, int minZ, int maxX, int maxY, int maxZ)
  {
    float A = integral(minX,minY,minZ);
    float B = integral(maxX,minY,minZ);
    float C = integral(maxX,minY,maxZ);
    float D = integral(minX,minY,maxZ);
    float E = integral(minX,maxY,minZ);
    float F = integral(maxX,maxY,minZ);
    float G = integral(maxX,maxY,maxZ);
    float H = integral(minX,maxY,maxZ);

    return (-A + B - C + D + E - F + G - H);
  }

  IntegralHistogram integral_histogram(const Image& image, int bins, float minValue, float maxValue)
  {
    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    IntegralHistogram histograms(sizeX+1,sizeY+1,sizeZ+1);

    histograms(0,0,0) = Histogram(bins, minValue, maxValue);

    for (int x = 1; x <= sizeX; x++)
    {
      histograms(x,0,0) = histograms(x-1,0,0);
      histograms(x,0,0).add(image(x-1,0,0));
    }

    for (int y = 1; y <= sizeY; y++)
    {
      histograms(0,y,0) = histograms(0,y-1,0);
      histograms(0,y,0).add(image(0,y-1,0));
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      histograms(0,0,z) = histograms(0,0,z-1);
      histograms(0,0,z).add(image(0,0,z-1));
    }

    for (int y = 1; y <= sizeY; y++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        histograms(x,y,0) = histograms(x-1,y,0);
        histograms(x,y,0) += histograms(x,y-1,0);
        histograms(x,y,0) -= histograms(x-1,y-1,0);
        histograms(x,y,0).add(image(x-1,y-1,0));
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        histograms(x,0,z) = histograms(x-1,0,z);
        histograms(x,0,z) += histograms(x,0,z-1);
        histograms(x,0,z) -= histograms(x-1,0,z-1);
        histograms(x,0,z).add(image(x-1,0,z-1));
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        histograms(0,y,z) = histograms(0,y-1,z);
        histograms(0,y,z) += histograms(0,y,z-1);
        histograms(0,y,z) -= histograms(0,y-1,z-1);
        histograms(0,y,z).add(image(0,y-1,z-1));
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        for (int x = 1; x <= sizeX; x++)
        {
          histograms(x,y,z) = histograms(x,y,z-1);
          histograms(x,y,z) += histograms(x,y-1,z);
          histograms(x,y,z) -= histograms(x,y-1,z-1);
          histograms(x,y,z) += histograms(x-1,y,z);
          histograms(x,y,z) -= histograms(x-1,y,z-1);
          histograms(x,y,z) -= histograms(x-1,y-1,z);
          histograms(x,y,z) += histograms(x-1,y-1,z-1);
          histograms(x,y,z).add(image(x-1,y-1,z-1));
        }
      }
    }

    return histograms;
  }

  IntegralHistogram integral_histogram_alternative(const Image& image, int bins, float minValue, float maxValue)
  {
    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    IntegralHistogram histograms(sizeX+1,sizeY+1,sizeZ+1);
    IntegralHistogram temp1(sizeX+1,sizeY+1,sizeZ+1);
    IntegralHistogram temp2(sizeX+1,sizeY+1,sizeZ+1);

    histograms(0,0,0) = Histogram(bins, minValue, maxValue);

    for (int x = 1; x <= sizeX; x++)
    {
      histograms(x,0,0) = Histogram(bins, minValue, maxValue);
    }

    for (int y = 1; y <= sizeY; y++)
    {
      histograms(0,y,0) = Histogram(bins, minValue, maxValue);
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      histograms(0,0,z) = Histogram(bins, minValue, maxValue);
    }

    for (int y = 1; y <= sizeY; y++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        temp1(x,y,0) = Histogram(bins, minValue, maxValue);
        histograms(x,y,0) = Histogram(bins, minValue, maxValue);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int x = 1; x <= sizeX; x++)
      {
        temp2(x,0,z) = Histogram(bins, minValue, maxValue);
        histograms(x,0,z) = Histogram(bins, minValue, maxValue);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        histograms(0,y,z) = Histogram(bins, minValue, maxValue);
      }
    }

    for (int z = 1; z <= sizeZ; z++)
    {
      for (int y = 1; y <= sizeY; y++)
      {
        for (int x = 1; x <= sizeX; x++)
        {
          temp1(x,y,z) = temp1(x,y,z-1);
          temp1(x,y,z).add(image(x-1,y-1,z-1));
          temp2(x,y,z) = temp1(x,y,z);
          temp2(x,y,z) += temp2(x,y-1,z);
          histograms(x,y,z) = temp2(x,y,z);
          histograms(x,y,z) += histograms(x-1,y,z);
        }
      }
    }

    return histograms;
  }

  void reverse_integral_histogram(const IntegralHistogram& integral, Image& image)
  {
    for (int z = 0; z < image.sizeZ(); z++)
    {
      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          Histogram h = evaluate_integral_histogram(integral, x, y, z, x+1, y+1, z+1);
          image(x,y,z) = static_cast<float>(maxCountIndex(h.counts()));
        }
      }
    }
  }

  Histogram evaluate_integral_histogram(const IntegralHistogram& integral, int minX, int minY, int minZ, int maxX, int maxY, int maxZ)
  {
    Histogram A = integral(minX,minY,minZ);
    Histogram B = integral(maxX,minY,minZ);
    Histogram C = integral(maxX,minY,maxZ);
    Histogram D = integral(minX,minY,maxZ);
    Histogram E = integral(minX,maxY,minZ);
    Histogram F = integral(maxX,maxY,minZ);
    Histogram G = integral(maxX,maxY,maxZ);
    Histogram H = integral(minX,maxY,maxZ);

    return (B - C + D + E - F + G - H - A);
  }

  void threshold(const Image& image, Image& binary, float minValue, float maxValue)
  {
    for (int z = 0; z < image.sizeZ(); z++)
    {
      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          auto value = image(x,y,z);
          auto result = (value >= minValue && value <= maxValue) ? 1.0f : 0.0f;
          binary(x,y,z) = result;
        }
      }
    }
  }

  void random_binary(Image& image, float probability)
  {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    for (auto& p : image)
      p = dist(mt) < probability ? 1.0f : 0.0f;
  }

  void invert_binary(const Image& input, Image& output)
  {
    for (int z = 0; z < input.sizeZ(); z++)
    {
      for (int y = 0; y < input.sizeY(); y++)
      {
        for (int x = 0; x < input.sizeX(); x++)
        {
          output(x,y,z) = input(x,y,z) != 0.0f ? 0.0f : 1.0f;
        }
      }
    }
  }

  void dilate_binary(const Image& input, Image& output)
  {
    copy(input, output);

    for (int z = 1; z < input.sizeZ()-1; z++)
    {
      for (int y = 1; y < input.sizeY()-1; y++)
      {
        for (int x = 1; x < input.sizeX()-1; x++)
        {
          if (input(x,y,z) != 0.0f)
          {
            for (int dz = -1; dz <= 1; dz++)
            {
              for (int dy = -1; dy <= 1; dy++)
              {
                for (int dx = -1; dx <= 1; dx++)
                {
                  output(x+dx,y+dy,z+dz) = 1.0f;
                }
              }
            }
          }
        }
      }
    }
  }

  void erode_binary(const Image& input, Image& output)
  {
    copy(input, output);

    for (int z = 1; z < input.sizeZ()-1; z++)
    {
      for (int y = 1; y < input.sizeY()-1; y++)
      {
        for (int x = 1; x < input.sizeX()-1; x++)
        {
          if (input(x,y,z) != 0.0f)
          {
            bool preserve = true;
            for (int dz = -1; dz <= 1; dz++)
            {
              for (int dy = -1; dy <= 1; dy++)
              {
                for (int dx = -1; dx <= 1; dx++)
                {
                  if (input(x+dx,y+dy,z+dz) == 0.0f)
                  {
                    preserve = false;
                  }
                  if (!preserve) break;
                }
                if (!preserve) break;
              }
              if (!preserve) break;
            }
            output(x,y,z) = preserve ? 1.0f : 0.0f;
          }
        }
      }
    }
  }

  void euclidean_distmap(const Image& binary, Image& distmap, int numPasses)
  {
    Image point_distances(3,3,3);
    for (int z = 0; z < 3; z++)
    {
      for (int y = 0; y < 3; y++)
      {
        for (int x = 0; x < 3; x++)
        {
          Eigen::Vector3d dist = Eigen::Vector3d(x-1,y-1,z-1).cwiseProduct(binary.spacing());
          point_distances(x,y,z) = static_cast<float>(dist.norm());
        }
      }
    }

    invert_binary(binary, distmap);
    mul(distmap, std::numeric_limits<float>::max());

    for (int i = 0; i < numPasses; i++)
    {
      impl::forward_pass(distmap, point_distances);
      impl::backward_pass(distmap, point_distances);
    }
  }

  void euclidean_distmap_fast(const Image& binary, Image& distmap, int numPasses)
  {
    Image point_distances(3,3,3);
    for (int z = 0; z < 3; z++)
    {
      for (int y = 0; y < 3; y++)
      {
        for (int x = 0; x < 3; x++)
        {
          Eigen::Vector3d dist = Eigen::Vector3d(x-1,y-1,z-1).cwiseProduct(binary.spacing());
          point_distances(x,y,z) = static_cast<float>(dist.norm());
        }
      }
    }

    invert_binary(binary, distmap);
    Image padded_distmap = pad(distmap, 1, 1, 1, 1, 1, 1, 1.0f);
    mul(padded_distmap, std::numeric_limits<float>::max());

    for (int i = 0; i < numPasses; i++)
    {
      impl::forward_pass_fast(padded_distmap, point_distances);
      impl::backward_pass_fast(padded_distmap, point_distances);
    }

    copy_region(padded_distmap, distmap, 1, 1, 1, 0, 0, 0, distmap.sizeX(), distmap.sizeY(), distmap.sizeZ());
  }

  void geodesic_distmap(const Image& image, const Image& binary, Image& distmap, float intensityWeight, int numPasses)
  {
    Image point_distances(3,3,3);
    for (int z = 0; z < 3; z++)
    {
      for (int y = 0; y < 3; y++)
      {
        for (int x = 0; x < 3; x++)
        {
          Eigen::Vector3d dist = Eigen::Vector3d(x-1,y-1,z-1).cwiseProduct(binary.spacing());
          point_distances(x,y,z) = static_cast<float>(dist.norm() * dist.norm());
        }
      }
    }

    invert_binary(binary, distmap);
    mul(distmap, std::numeric_limits<float>::max());

    for (int i = 0; i < numPasses; i++)
    {
      impl::forward_pass(image, distmap, point_distances, intensityWeight);
      impl::backward_pass(image, distmap, point_distances, intensityWeight);
    }
  }

  void geodesic_distmap_fast(const Image& image, const Image& binary, Image& distmap, float intensityWeight, int numPasses)
  {
    Image point_distances(3,3,3);
    for (int z = 0; z < 3; z++)
    {
      for (int y = 0; y < 3; y++)
      {
        for (int x = 0; x < 3; x++)
        {
          Eigen::Vector3d dist = Eigen::Vector3d(x-1,y-1,z-1).cwiseProduct(binary.spacing());
          point_distances(x,y,z) = static_cast<float>(dist.norm() * dist.norm());
        }
      }
    }

    invert_binary(binary, distmap);
    Image padded_image = pad(image, 1, 1, 1, 1, 1, 1, 0.0f);
    Image padded_distmap = pad(distmap, 1, 1, 1, 1, 1, 1, 1.0f);
    mul(padded_distmap, std::numeric_limits<float>::max());

    for (int i = 0; i < numPasses; i++)
    {
      impl::forward_pass_fast(padded_image, padded_distmap, point_distances, intensityWeight);
      impl::backward_pass_fast(padded_image, padded_distmap, point_distances, intensityWeight);
    }

    copy_region(padded_distmap, distmap, 1, 1, 1, 0, 0, 0, distmap.sizeX(), distmap.sizeY(), distmap.sizeZ());
  }

  Image extract_slice(const Image& image, ImagePlane plane, int sliceNumber)
  {
    if (plane == ImagePlane::XY) //axial
    {
      Image slice(image.sizeX(), image.sizeY(), 1);
      slice.dataType(image.dataType());

      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          slice(x,y,0) = image(x,y,sliceNumber);
        }
      }
      return slice;
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      Image slice(image.sizeX(), image.sizeZ(), 1);
      slice.dataType(image.dataType());

      for (int z = 0; z < image.sizeZ(); z++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          slice(x,z,0) = image(x,sliceNumber,image.sizeZ()-z-1);
        }
      }
      return slice;
    }
    else //sagittal
    {
      Image slice(image.sizeY(), image.sizeZ(), 1);
      slice.dataType(image.dataType());

      for (int z = 0; z < image.sizeZ(); z++)
      {
        for (int y = 0; y < image.sizeY(); y++)
        {
          slice(y,z,0) = image(sliceNumber,y,image.sizeZ()-z-1);
        }
      }
      return slice;
    }
  }

  void insert_slice(const Image& slice, Image& image, ImagePlane plane, int sliceNumber)
  {
    if (plane == ImagePlane::XY) //axial
    {
      for (int y = 0; y < image.sizeY(); y++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          image(x,y,sliceNumber) = slice(x,y,0);
        }
      }
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      for (int z = 0; z < image.sizeZ(); z++)
      {
        for (int x = 0; x < image.sizeX(); x++)
        {
          image(x,sliceNumber,image.sizeZ()-z-1) = slice(x,z,0);
        }
      }
    }
    else //sagittal
    {
      for (int z = 0; z < image.sizeZ(); z++)
      {
        for (int y = 0; y < image.sizeY(); y++)
        {
          image(sliceNumber,y,image.sizeZ()-z-1) = slice(y,z,0);
        }
      }
    }
  }

  void texture_8bit_gray(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level)
  {
    float wlMin = level - window / 2.0f;

    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    if (plane == ImagePlane::XY) //axial
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          float value = (image(x,y,sliceNumber) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[x + y*sizeX] = static_cast<unsigned char>(value);
        }
      }
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          float value = (image(x,sliceNumber,sizeZ-z-1) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[x + z*sizeX] = static_cast<unsigned char>(value);
        }
      }
    }
    else //sagittal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          float value = (image(sliceNumber,y,sizeZ-z-1) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[y + z*sizeY] = static_cast<unsigned char>(value);
        }
      }
    }
  }

  void texture_32bit_gray(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level)
  {
    float wlMin = level - window / 2.0f;

    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    if (plane == ImagePlane::XY) //axial
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          float value = (image(x,y,sliceNumber) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[x*4 + y*sizeX*4 + 0] = static_cast<unsigned char>(value);
          texture[x*4 + y*sizeX*4 + 1] = static_cast<unsigned char>(value);
          texture[x*4 + y*sizeX*4 + 2] = static_cast<unsigned char>(value);
          texture[x*4 + y*sizeX*4 + 3] = 255;
        }
      }
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          float value = (image(x,sliceNumber,sizeZ-z-1) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[x*4 + z*sizeX*4 + 0] = static_cast<unsigned char>(value);
          texture[x*4 + z*sizeX*4 + 1] = static_cast<unsigned char>(value);
          texture[x*4 + z*sizeX*4 + 2] = static_cast<unsigned char>(value);
          texture[x*4 + z*sizeX*4 + 3] = 255;
        }
      }
    }
    else //sagittal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          float value = (image(sliceNumber,y,sizeZ-z-1) - wlMin) / window * 255.0f;
          if (value > 255.0f) value = 255.0f;
          else if (value < 0.0f) value = 0.0f;
          texture[y*4 + z*sizeY*4 + 0] = static_cast<unsigned char>(value);
          texture[y*4 + z*sizeY*4 + 1] = static_cast<unsigned char>(value);
          texture[y*4 + z*sizeY*4 + 2] = static_cast<unsigned char>(value);
          texture[y*4 + z*sizeY*4 + 3] = 255;
        }
      }
    }
  }

  void texture_32bit_color_jet64(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float window, float level, float alpha)
  {
    float wlMin = level - window / 2.0f;

    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    if (plane == ImagePlane::XY) //axial
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          int value = static_cast<int>((image(x,y,sliceNumber) - wlMin) / window * 63);
          if (value > 63) value = 63;
          else if (value < 0) value = 0;
          texture[x*4 + y*sizeX*4 + 0] = static_cast<unsigned char>(colormaps::jet64[63-value][0] * 255);
          texture[x*4 + y*sizeX*4 + 1] = static_cast<unsigned char>(colormaps::jet64[63-value][1] * 255);
          texture[x*4 + y*sizeX*4 + 2] = static_cast<unsigned char>(colormaps::jet64[63-value][2] * 255);
          texture[x*4 + y*sizeX*4 + 3] = static_cast<unsigned char>(colormaps::jet64[63-value][3] * 255 * alpha);
        }
      }
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          int value = static_cast<int>((image(x,sliceNumber,sizeZ-z-1) - wlMin) / window * 63);
          if (value > 63) value = 63;
          else if (value < 0) value = 0;
          texture[x*4 + z*sizeX*4 + 0] = static_cast<unsigned char>(colormaps::jet64[63-value][0] * 255);
          texture[x*4 + z*sizeX*4 + 1] = static_cast<unsigned char>(colormaps::jet64[63-value][1] * 255);
          texture[x*4 + z*sizeX*4 + 2] = static_cast<unsigned char>(colormaps::jet64[63-value][2] * 255);
          texture[x*4 + z*sizeX*4 + 3] = static_cast<unsigned char>(colormaps::jet64[63-value][3] * 255 * alpha);
        }
      }
    }
    else //sagittal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          int value = static_cast<int>((image(sliceNumber,y,sizeZ-z-1) - wlMin) / window * 63);
          if (value > 63) value = 63;
          else if (value < 0) value = 0;
          texture[y*4 + z*sizeY*4 + 0] = static_cast<unsigned char>(colormaps::jet64[63-value][0] * 255);
          texture[y*4 + z*sizeY*4 + 1] = static_cast<unsigned char>(colormaps::jet64[63-value][1] * 255);
          texture[y*4 + z*sizeY*4 + 2] = static_cast<unsigned char>(colormaps::jet64[63-value][2] * 255);
          texture[y*4 + z*sizeY*4 + 3] = static_cast<unsigned char>(colormaps::jet64[63-value][3] * 255 * alpha);
        }
      }
    }
  }

  void texture_32bit_color_indexed32(const Image& image, unsigned char* texture, ImagePlane plane, int sliceNumber, float alpha)
  {
    int sizeX = image.sizeX();
    int sizeY = image.sizeY();
    int sizeZ = image.sizeZ();

    if (plane == ImagePlane::XY) //axial
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          int value = static_cast<int>(image(x,y,sliceNumber)) % 32;
          value = (value + 32) % 32;
          texture[x*4 + y*sizeX*4 + 0] = static_cast<unsigned char>(colormaps::indexed32[value][0] * 255);
          texture[x*4 + y*sizeX*4 + 1] = static_cast<unsigned char>(colormaps::indexed32[value][1] * 255);
          texture[x*4 + y*sizeX*4 + 2] = static_cast<unsigned char>(colormaps::indexed32[value][2] * 255);
          texture[x*4 + y*sizeX*4 + 3] = static_cast<unsigned char>(colormaps::indexed32[value][3] * 255 * alpha);
        }
      }
    }
    else if (plane == ImagePlane::XZ) //coronal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          int value = static_cast<int>(image(x,sliceNumber,sizeZ-z-1)) % 32;
          value = (value + 32) % 32;
          texture[x*4 + z*sizeX*4 + 0] = static_cast<unsigned char>(colormaps::indexed32[value][0] * 255);
          texture[x*4 + z*sizeX*4 + 1] = static_cast<unsigned char>(colormaps::indexed32[value][1] * 255);
          texture[x*4 + z*sizeX*4 + 2] = static_cast<unsigned char>(colormaps::indexed32[value][2] * 255);
          texture[x*4 + z*sizeX*4 + 3] = static_cast<unsigned char>(colormaps::indexed32[value][3] * 255 * alpha);
        }
      }
    }
    else //sagittal
    {
      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          int value = static_cast<int>(image(sliceNumber,y,sizeZ-z-1)) % 32;
          value = (value + 32) % 32;
          texture[y*4 + z*sizeY*4 + 0] = static_cast<unsigned char>(colormaps::indexed32[value][0] * 255);
          texture[y*4 + z*sizeY*4 + 1] = static_cast<unsigned char>(colormaps::indexed32[value][1] * 255);
          texture[y*4 + z*sizeY*4 + 2] = static_cast<unsigned char>(colormaps::indexed32[value][2] * 255);
          texture[y*4 + z*sizeY*4 + 3] = static_cast<unsigned char>(colormaps::indexed32[value][3] * 255 * alpha);
        }
      }
    }
  }

  void checkerboard(Image& image, int patchSizeX, int patchSizeY, int patchSizeZ)
  {
    for (int i=0; i<image.sizeX()/patchSizeX; ++i)
    {
      for (int j=0; j<image.sizeY()/patchSizeY; ++j)
      {
        for (int k=0; k<image.sizeZ()/patchSizeZ; ++k)
        {
          int x = i*patchSizeX;
          int y = j*patchSizeY;
          int z = k*patchSizeZ;
          auto patch = subimage(image, x, y, z, patchSizeX, patchSizeY, patchSizeZ);
          if ((i%2 == 0 && j%2==1 && k%2==0) || (i%2==1 && j%2==0 && k%2==0) || (i%2 == 1 && j%2==1 && k%2==1) || (i%2==0 && j%2==0 && k%2==1))
            zeros(patch);
          else
            ones(patch);
        }
      }
    }
  }

  Image checkerboard(const Image& inputA, const Image& inputB, int patchSizeX, int patchSizeY, int patchSizeZ)
  {
    Image image = inputA.clone();
    for (int i=0; i<image.sizeX()/patchSizeX; ++i)
    {
      for (int j=0; j<image.sizeY()/patchSizeY; ++j)
      {
        for (int k=0; k<image.sizeZ()/patchSizeZ; ++k)
        {
          int x = i*patchSizeX;
          int y = j*patchSizeY;
          int z = k*patchSizeZ;
          auto patch = subimage(image, x, y, z, patchSizeX, patchSizeY, patchSizeZ);
          if ((i%2 == 0 && j%2==1 && k%2==0) || (i%2==1 && j%2==0 && k%2==0) || (i%2 == 1 && j%2==1 && k%2==1) || (i%2==0 && j%2==0 && k%2==1))
          {
            auto patchA = subimage(const_cast<Image&>(inputA), x, y, z, patchSizeX, patchSizeY, patchSizeZ);
            copy(patchA, patch);
          }
          else
          {
            auto patchB = subimage(const_cast<Image&>(inputB), x, y, z, patchSizeX, patchSizeY, patchSizeZ);
            copy(patchB, patch);
          }
        }
      }
    }

    return image;
  }

  namespace impl
  {
    template <typename Sampler>
    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, float outsideValue, Sampler sampler)
    {
      Eigen::Vector3d deltaX = Eigen::Vector3d(transformInImageSpace(0,0), transformInImageSpace(1,0), transformInImageSpace(2,0));
      Eigen::Vector3d deltaY = Eigen::Vector3d(transformInImageSpace(0,1), transformInImageSpace(1,1), transformInImageSpace(2,1));

      for (int z = 0; z < output.sizeZ(); z++)
      {
        Eigen::Vector4d p_h = transformInImageSpace * Eigen::Vector4d(0, 0, z, 1);
        Eigen::Vector3d p = p_h.segment(0,3);
        for (int y = 0; y < output.sizeY(); y++)
        {
          Eigen::Vector3d pixel = p;
          for (int x = 0; x < output.sizeX(); x++)
          {
            output(x,y,z) = sampler(pixel, outsideValue);
            pixel += deltaX;
          }
          p += deltaY;
        }
      }
    }

    template<typename Sampler>
    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, const Image& fieldX, const Image& fieldY, const Image& fieldZ, float outsideValue, Sampler sampler)
    {
      for (int z = 0; z < output.sizeZ(); z++)
      {
        for (int y = 0; y < output.sizeY(); y++)
        {
          for (int x = 0; x < output.sizeX(); x++)
          {
            Eigen::Vector4d p = transformInImageSpace * Eigen::Vector4d(x + fieldX(x,y,z) / output.spacing()[0], y + fieldY(x,y,z) / output.spacing()[1], z + fieldZ(x,y,z) / output.spacing()[2], 1);
            Eigen::Vector3d pixel = p.segment(0,3);
            output(x,y,z) = sampler(pixel, outsideValue);
          }
        }
      }
    }

    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, Interpolation interpolation, float outsideValue)
    {
      if (interpolation == LINEAR)
      {
        transform(input, output, transformInImageSpace, outsideValue, [&input](const Eigen::Vector3d& p, float f)
        {
          return input.linear(p,f);
        });
      }
      else
      {
        transform(input, output, transformInImageSpace, outsideValue, [&input](const Eigen::Vector3d& p, float f)
        {
          return input.nearest(p,f);
        });
      }
    }
    void transform(const Image& input, Image& output, const Eigen::Matrix4d& transformInImageSpace, const Image& fieldX, const Image& fieldY, const Image& fieldZ, Interpolation interpolation, float outsideValue)
    {
      if (interpolation == LINEAR)
      {
        transform(input, output, transformInImageSpace, fieldX, fieldY, fieldZ, outsideValue, [&input](const Eigen::Vector3d& p, float f)
        {
          return input.linear(p,f);
        });
      }
      else
      {
        transform(input, output, transformInImageSpace, fieldX, fieldY, fieldZ, outsideValue, [&input](const Eigen::Vector3d& p, float f)
        {
          return input.nearest(p,f);
        });
      }
    }

    std::vector<double> generate_gauss_filter(double sigmaInPixels)
    {
      int filterRadius = static_cast<int>(3 * sigmaInPixels);
      int filterSize = filterRadius * 2 + 1;

      double filterSum = 0.0;
      std::vector<double> gaussFilter(filterSize);
      double sigmaSq = sigmaInPixels * sigmaInPixels;

      if (filterSize == 1)
      {
        gaussFilter[0] = 1.0;
      }
      else
      {
        for (int i = -filterRadius; i <= filterRadius; i++)
        {
          gaussFilter[i + filterRadius] = 1.0 / (sigmaInPixels * sqrt(2 * MIA_PI)) * exp(-static_cast<double>(i * i) / (2.0 * sigmaSq));
          filterSum += gaussFilter[i + filterRadius];
        }

        for (int i = 0; i < filterSize; i++)
        {
          gaussFilter[i] /= filterSum;
        }
      }

      return gaussFilter;
    }

    void forward_pass(Image& distmap, const Image& point_distances)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          for (int x = 0; x < sizeX; x++)
          {
            float minDist = distmap(x,y,z);

            for (int dy = -1; dy <= 1; dy++)
            {
              for (int dx = -1; dx <= 1; dx++)
              {
                if (x + dx >= 0 && x + dx < sizeX && y + dy >= 0 && y + dy < sizeY && z - 1 >= 0 && z - 1 < sizeZ)
                {
                  float dist = distmap(x+dx,y+dy,z-1) + point_distances(dx+1,dy+1,0);
                  if (dist < minDist) minDist = dist;
                }
              }
            }

            for (int dx = -1; dx <= 1; dx++)
            {
              if (x + dx >= 0 && x + dx < sizeX && y - 1 >= 0 && y - 1 < sizeY)
              {
                float dist = distmap(x+dx,y-1,z) + point_distances(dx+1,0,1);
                if (dist < minDist) minDist = dist;
              }
            }

            {
              if (x - 1 >= 0 && x - 1 < sizeX)
              {
                float dist = distmap(x-1,y,z) + point_distances(0,1,1);
                if (dist < minDist) minDist = dist;
              }
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }
    void backward_pass(Image& distmap, const Image& point_distances)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = sizeZ-1; z >= 0; z--)
      {
        for (int y = sizeY-1; y >= 0; y--)
        {
          for (int x = sizeX-1; x >= 0; x--)
          {
            float minDist = distmap(x,y,z);

            for (int dy = 1; dy >= -1; dy--)
            {
              for (int dx = 1; dx >= -1; dx--)
              {
                if (x + dx >= 0 && x + dx < sizeX && y + dy >= 0 && y + dy < sizeY && z + 1 >= 0 && z + 1 < sizeZ)
                {
                  float dist = distmap(x+dx,y+dy,z+1) + point_distances(dx+1,dy+1,2);
                  if (dist < minDist) minDist = dist;
                }
              }
            }

            for (int dx = 1; dx >= -1; dx--)
            {
              if (x + dx >= 0 && x + dx < sizeX && y + 1 >= 0 && y + 1 < sizeY)
              {
                float dist = distmap(x+dx,y+1,z) + point_distances(dx+1,2,1);
                if (dist < minDist) minDist = dist;
              }
            }

            {
              if (x - 1 >= 0 && x + 1 < sizeX)
              {
                float dist = distmap(x+1,y,z) + point_distances(2,1,1);
                if (dist < minDist) minDist = dist;
              }
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }

    void forward_pass_fast(Image& distmap, const Image& point_distances)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = 1; z < sizeZ-1; z++)
      {
        for (int y = 1; y < sizeY-1; y++)
        {
          for (int x = 1; x < sizeX-1; x++)
          {
            float minDist = distmap(x,y,z);

            for (int dy = -1; dy <= 1; dy++)
            {
              for (int dx = -1; dx <= 1; dx++)
              {
                float dist = distmap(x+dx,y+dy,z-1) + point_distances(dx+1,dy+1,0);
                if (dist < minDist) minDist = dist;
              }
            }

            for (int dx = -1; dx <= 1; dx++)
            {
              float dist = distmap(x+dx,y-1,z) + point_distances(dx+1,0,1);
              if (dist < minDist) minDist = dist;
            }

            {
              float dist = distmap(x-1,y,z) + point_distances(0,1,1);
              if (dist < minDist) minDist = dist;
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }
    void backward_pass_fast(Image& distmap, const Image& point_distances)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = sizeZ-2; z >= 1; z--)
      {
        for (int y = sizeY-2; y >= 1; y--)
        {
          for (int x = sizeX-2; x >= 1; x--)
          {
            float minDist = distmap(x,y,z);

            for (int dy = 1; dy >= -1; dy--)
            {
              for (int dx = 1; dx >= -1; dx--)
              {
                float dist = distmap(x+dx,y+dy,z+1) + point_distances(dx+1,dy+1,2);
                if (dist < minDist) minDist = dist;
              }
            }

            for (int dx = 1; dx >= -1; dx--)
            {
              float dist = distmap(x+dx,y+1,z) + point_distances(dx+1,2,1);
              if (dist < minDist) minDist = dist;
            }

            {
              float dist = distmap(x+1,y,z) + point_distances(2,1,1);
              if (dist < minDist) minDist = dist;
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }

    void forward_pass(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = 0; z < sizeZ; z++)
      {
        for (int y = 0; y < sizeY; y++)
        {
          for (int x = 0; x < sizeX; x++)
          {
            float intensity = image(x,y,z);
            float minDist = distmap(x,y,z);

            for (int dy = -1; dy <= 1; dy++)
            {
              for (int dx = -1; dx <= 1; dx++)
              {
                if (x + dx >= 0 && x + dx < sizeX && y + dy >= 0 && y + dy < sizeY && z - 1 >= 0 && z - 1 < sizeZ)
                {
                  float grad = (intensity - image(x+dx,y+dy,z-1)) * intensityWeight;
                  float dist = distmap(x+dx,y+dy,z-1) + sqrt(point_distances(dx+1,dy+1,0) + grad*grad);
                  if (dist < minDist) minDist = dist;
                }
              }
            }

            for (int dx = -1; dx <= 1; dx++)
            {
              if (x + dx >= 0 && x + dx < sizeX && y - 1 >= 0 && y - 1 < sizeY)
              {
                float grad = (intensity - image(x+dx,y-1,z)) * intensityWeight;
                float dist = distmap(x+dx,y-1,z) + sqrt(point_distances(dx+1,0,1) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            {
              if (x - 1 >= 0 && x - 1 < sizeX)
              {
                float grad = (intensity - image(x-1,y,z)) * intensityWeight;
                float dist = distmap(x-1,y,z) + sqrt(point_distances(0,1,1) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }
    void backward_pass(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = sizeZ-1; z >= 0; z--)
      {
        for (int y = sizeY-1; y >= 0; y--)
        {
          for (int x = sizeX-1; x >= 0; x--)
          {
            float intensity = image(x,y,z);
            float minDist = distmap(x,y,z);

            for (int dy = 1; dy >= -1; dy--)
            {
              for (int dx = 1; dx >= -1; dx--)
              {
                if (x + dx >= 0 && x + dx < sizeX && y + dy >= 0 && y + dy < sizeY && z + 1 >= 0 && z + 1 < sizeZ)
                {
                  float grad = (intensity - image(x+dx,y+dy,z+1)) * intensityWeight;
                  float dist = distmap(x+dx,y+dy,z+1) + sqrt(point_distances(dx+1,dy+1,2) + grad*grad);
                  if (dist < minDist) minDist = dist;
                }
              }
            }

            for (int dx = 1; dx >= -1; dx--)
            {
              if (x + dx >= 0 && x + dx < sizeX && y + 1 >= 0 && y + 1 < sizeY)
              {
                float grad = (intensity - image(x+dx,y+1,z)) * intensityWeight;
                float dist = distmap(x+dx,y+1,z) + sqrt(point_distances(dx+1,2,1) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            {
              if (x - 1 >= 0 && x + 1 < sizeX)
              {
                float grad = (intensity - image(x+1,y,z)) * intensityWeight;
                float dist = distmap(x+1,y,z) + sqrt(point_distances(2,1,1) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }

    void forward_pass_fast(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = 1; z < sizeZ-1; z++)
      {
        for (int y = 1; y < sizeY-1; y++)
        {
          for (int x = 1; x < sizeX-1; x++)
          {
            float intensity = image(x,y,z);
            float minDist = distmap(x,y,z);

            for (int dy = -1; dy <= 1; dy++)
            {
              for (int dx = -1; dx <= 1; dx++)
              {
                float grad = (intensity - image(x+dx,y+dy,z-1)) * intensityWeight;
                float dist = distmap(x+dx,y+dy,z-1) + sqrt(point_distances(dx+1,dy+1,0) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            for (int dx = -1; dx <= 1; dx++)
            {
              float grad = (intensity - image(x+dx,y-1,z)) * intensityWeight;
              float dist = distmap(x+dx,y-1,z) + sqrt(point_distances(dx+1,0,1) + grad*grad);
              if (dist < minDist) minDist = dist;
            }

            {
              float grad = (intensity - image(x-1,y,z)) * intensityWeight;
              float dist = distmap(x-1,y,z) + sqrt(point_distances(0,1,1) + grad*grad);
              if (dist < minDist) minDist = dist;
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }
    void backward_pass_fast(const Image& image, Image& distmap, const Image& point_distances, float intensityWeight)
    {
      int sizeX = distmap.sizeX();
      int sizeY = distmap.sizeY();
      int sizeZ = distmap.sizeZ();

      for (int z = sizeZ-2; z >= 1; z--)
      {
        for (int y = sizeY-2; y >= 1; y--)
        {
          for (int x = sizeX-2; x >= 1; x--)
          {
            float intensity = image(x,y,z);
            float minDist = distmap(x,y,z);

            for (int dy = 1; dy >= -1; dy--)
            {
              for (int dx = 1; dx >= -1; dx--)
              {
                float grad = (intensity - image(x+dx,y+dy,z+1)) * intensityWeight;
                float dist = distmap(x+dx,y+dy,z+1) + sqrt(point_distances(dx+1,dy+1,2) + grad*grad);
                if (dist < minDist) minDist = dist;
              }
            }

            for (int dx = 1; dx >= -1; dx--)
            {
              float grad = (intensity - image(x+dx,y+1,z)) * intensityWeight;
              float dist = distmap(x+dx,y+1,z) + sqrt(point_distances(dx+1,2,1) + grad*grad);
              if (dist < minDist) minDist = dist;
            }

            {
              float grad = (intensity - image(x+1,y,z)) * intensityWeight;
              float dist = distmap(x+1,y,z) + sqrt(point_distances(2,1,1) + grad*grad);
              if (dist < minDist) minDist = dist;
            }

            distmap(x,y,z) = minDist;
          }
        }
      }
    }
  }
}
