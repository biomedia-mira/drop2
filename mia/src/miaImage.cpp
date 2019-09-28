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

#include "miaImage.h"

namespace mia
{
  Image::Image()
    : m_datatype(FLOAT)
    , m_size_x(0)
    , m_size_y(0)
    , m_size_z(0)
    , m_step_x(0)
    , m_step_y(0)
    , m_step_z(0)
    , m_spacing()
    , m_origin()
    , m_image_to_world()
    , m_world_to_image()
    , m_data( nullptr )
  {}

  Image::Image(int sizeX, int sizeY, int sizeZ)
    : m_datatype(FLOAT)
    , m_size_x(sizeX)
    , m_size_y(sizeY)
    , m_size_z(sizeZ)
    , m_step_x(1)
    , m_step_y(sizeX)
    , m_step_z(sizeX*sizeY)
    , m_spacing(1.0, 1.0, 1.0)
    , m_origin(0.0, 0.0, 0.0)
    , m_image_to_world(Eigen::Matrix3d::Identity())
    , m_world_to_image(Eigen::Matrix3d::Identity())
    , m_data( new float[sizeX*sizeY*sizeZ], [](float* d)
  {
    delete [] d;
  } )
  {}

  Image::Image(float* data, int sizeX, int sizeY, int sizeZ)
    : m_datatype(FLOAT)
    , m_size_x(sizeX)
    , m_size_y(sizeY)
    , m_size_z(sizeZ)
    , m_step_x(1)
    , m_step_y(sizeX)
    , m_step_z(sizeX*sizeY)
    , m_spacing(1.0, 1.0, 1.0)
    , m_origin(0.0, 0.0, 0.0)
    , m_image_to_world(Eigen::Matrix3d::Identity())
    , m_world_to_image(Eigen::Matrix3d::Identity())
    , m_data( data, [](float*) {} )
  {}

  Image::Image(float* data, int sizeX, int sizeY, int sizeZ, int stepX, int stepY, int stepZ)
    : m_datatype(FLOAT)
    , m_size_x(sizeX)
    , m_size_y(sizeY)
    , m_size_z(sizeZ)
    , m_step_x(stepX)
    , m_step_y(stepY)
    , m_step_z(stepZ)
    , m_spacing(1.0, 1.0, 1.0)
    , m_origin(0.0, 0.0, 0.0)
    , m_image_to_world(Eigen::Matrix3d::Identity())
    , m_world_to_image(Eigen::Matrix3d::Identity())
    , m_data( data, [](float*) {} )
  {}

  Image::Image(const Image& other)
    : m_datatype(other.m_datatype)
    , m_size_x(other.m_size_x)
    , m_size_y(other.m_size_y)
    , m_size_z(other.m_size_z)
    , m_step_x(other.m_step_x)
    , m_step_y(other.m_step_y)
    , m_step_z(other.m_step_z)
    , m_spacing(other.m_spacing)
    , m_origin(other.m_origin)
    , m_image_to_world(other.m_image_to_world)
    , m_world_to_image(other.m_world_to_image)
    , m_data(other.m_data)
  {}

  Image::Image(Image&& other)
    : m_datatype(std::move(other.m_datatype))
    , m_size_x(std::move(other.m_size_x))
    , m_size_y(std::move(other.m_size_y))
    , m_size_z(std::move(other.m_size_z))
    , m_step_x(std::move(other.m_step_x))
    , m_step_y(std::move(other.m_step_y))
    , m_step_z(std::move(other.m_step_z))
    , m_spacing(std::move(other.m_spacing))
    , m_origin(std::move(other.m_origin))
    , m_image_to_world(std::move(other.m_image_to_world))
    , m_world_to_image(std::move(other.m_world_to_image))
    , m_data(std::move(other.m_data))
  {}

  Image& Image::operator=(const Image& other)
  {
    if (this != &other)
    {
      *this = Image(other);
    }
    return *this;
  }

  Image& Image::operator=(Image&& other)
  {
    //if (this != &other)
    //{
    //  *this = Image(other);
    //}
    m_datatype = std::move(other.m_datatype);
    m_size_x = std::move(other.m_size_x);
    m_size_y = std::move(other.m_size_y);
    m_size_z = std::move(other.m_size_z);
    m_step_x = std::move(other.m_step_x);
    m_step_y = std::move(other.m_step_y);
    m_step_z = std::move(other.m_step_z);
    m_spacing = std::move(other.m_spacing);
    m_origin = std::move(other.m_origin);
    m_image_to_world = std::move(other.m_image_to_world);
    m_world_to_image = std::move(other.m_world_to_image);
    m_data = std::move(other.m_data);
    return *this;
  }

  Image Image::clone() const
  {
    Image image(m_size_x, m_size_y, m_size_z);
    image.m_datatype = m_datatype;
    image.m_spacing = m_spacing;
    image.m_origin = m_origin;
    image.m_image_to_world = m_image_to_world;
    image.m_world_to_image = m_world_to_image;

    if (contiguous())
    {
      const float* src_start = &(*this)(0,0,0);
      float* dst_start = &image(0,0,0);
      std::copy(src_start, src_start+size(), dst_start);
    }
    else
    {
      for (int z = 0; z < sizeZ(); z++)
      {
        for (int y = 0; y < sizeY(); y++)
        {
          if (m_step_x == 1)
          {
            const auto* src_start = &(*this)(0,y,z);
            auto* dst_start = &image(0,y,z);
            std::copy(src_start, src_start+sizeX(), dst_start);
          }
          else
          {
            for (int x = 0; x < sizeX(); x++)
            {
              image(x,y,z) = (*this)(x,y,z);
            }
          }
        }
      }
    }
    return image;
  }

  void Image::imageToWorld(const Eigen::Matrix3d& imageToWorld)
  {
    m_image_to_world = imageToWorld;
    m_world_to_image = imageToWorld.transpose();
  }

  const Eigen::Matrix3d& Image::imageToWorld() const
  {
    return m_image_to_world;
  }

  void Image::worldToImage(const Eigen::Matrix3d& worldToImage)
  {
    m_world_to_image = worldToImage;
    m_image_to_world = worldToImage.transpose();
  }

  const Eigen::Matrix3d& Image::worldToImage() const
  {
    return m_world_to_image;
  }

  Eigen::Vector3d Image::mapImageToWorld(const Eigen::Vector3d& point) const
  {
    return m_image_to_world * (point.cwiseProduct(m_spacing)) + m_origin;
  }

  Eigen::Vector3d Image::mapWorldToImage(const Eigen::Vector3d& point) const
  {
    return (m_world_to_image * (point - m_origin)).cwiseQuotient(m_spacing);
  }

  Eigen::Vector3d Image::imageCenterToWorld() const
  {
    Eigen::Vector3d imageCenter(static_cast<double>(m_size_x - 1) / 2.0, static_cast<double>(m_size_y - 1) / 2.0, static_cast<double>(m_size_z - 1) / 2.0);
    return mapImageToWorld(imageCenter);
  }

  Eigen::Matrix4d Image::imageToWorldTransform() const
  {
    Eigen::Matrix3d scaling = m_spacing.asDiagonal();
    Eigen::Matrix3d transform = m_image_to_world * scaling;

    Eigen::Matrix4d fullTransform = Eigen::Matrix4d::Identity();
    fullTransform.block(0, 0, 3, 3) = transform.block(0, 0, 3, 3);
    fullTransform.block(0, 3, 3, 1) = m_origin;

    return fullTransform;
  }

  Eigen::Matrix4d Image::worldToImageTransform() const
  {
    Eigen::Matrix3d scaling = m_spacing.cwiseInverse().asDiagonal();
    Eigen::Matrix3d transform = scaling * m_world_to_image;

    Eigen::Matrix4d fullTransform = Eigen::Matrix4d::Identity();
    fullTransform.block(0, 0, 3, 3) = transform;

    Eigen::Matrix4d translation = Eigen::Matrix4d::Identity();
    translation.block(0, 3, 3, 1) = -m_origin;

    return fullTransform * translation;
  }

  float Image::nearest(const Eigen::Vector3d& point, float outsideValue) const
  {
    if (point[0] < -0.5 || point[1] < -0.5 || point[2] < -0.5 || point[0] >= static_cast<double>(m_size_x) - 0.5 || point[1] >= static_cast<double>(m_size_y) - 0.5 || point[2] >= static_cast<double>(m_size_z) - 0.5)
    {
      return outsideValue;
    }

    int x = static_cast<int>(point[0] + 0.5);
    int y = static_cast<int>(point[1] + 0.5);
    int z = static_cast<int>(point[2] + 0.5);

    return (*this)(x,y,z);
  }

  float Image::linear(const Eigen::Vector3d& point, float outsideValue) const
  {
    int xx = static_cast<int>(point[0]);
    int yy = static_cast<int>(point[1]);
    int zz = static_cast<int>(point[2]);

    auto data = m_data.get();

    double dx1 = point[0] - xx, dy1 = point[1] - yy, dz1 = point[2] - zz;

    // boundary check
    if (point[0] < 0 || point[1] < 0 || point[2] < 0 || point[0] >= m_size_x - 1 || point[1] >= m_size_y - 1 || point[2] >= m_size_z - 1)
    {
      if (point[0] < -0.5 || point[1] < -0.5 || point[2] < -0.5 || point[0] >= static_cast<double>(m_size_x) - 0.5 || point[1] >= static_cast<double>(m_size_y) - 0.5 || point[2] >= static_cast<double>(m_size_z) - 0.5)
      {
        return outsideValue;
      }
      else
      {
        if (point[0] < 0)
        {
          dx1 = 0;
          xx = 0;
        }
        else if (point[0] > m_size_x - 1)
        {
          dx1 = 0;
        }

        if (point[1] < 0)
        {
          dy1 = 0;
          yy = 0;
        }
        else if (point[1] > m_size_y - 1)
        {
          dy1 = 0;
        }

        if (point[2] < 0)
        {
          dz1 = 0;
          zz = 0;
        }
        else if (point[2] > m_size_z - 1)
        {
          dz1 = 0;
        }

        double dx2 = 1.0 - dx1, dy2 = 1.0 - dy1, dz2 = 1.0 - dz1;
        int index = xx * m_step_x + yy * m_step_y + zz * m_step_z;

        // this is always valid, return only this if x==_dimX and y==_dimY and z==_dimZ
        double value = data[index] * dx2 * dy2 * dz2;

        // check all cases where at least one of the dx1,dy1,dz1 values is 0, but not all are 0.

        if (dx1 == 0 && dy1 == 0 && dz1 != 0)
        {
          // 0 0 1
          value += data[index + m_step_z] * dz1;
        }
        else if (dx1 == 0 && dy1 != 0 && dz1 == 0)
        {
          // 0 1 0
          value += data[index + m_step_y] * dy1;
        }
        else if (dx1 != 0 && dy1 == 0 && dz1 == 0)
        {
          // 1 0 0
          value += data[index + m_step_x] * dx1;
        }
        else if (dx1 == 0 && dy1 != 0 && dz1 != 0)
        {
          // 0 1 1
          value += data[index + m_step_y + m_step_z] * dy1 * dz1;
          value += data[index + m_step_y] * dy1 * dz2;
          value += data[index + m_step_z] * dy2 * dz1;
        }
        else if (dx1 != 0 && dy1 == 0 && dz1 != 0)
        {
          // 1 0 1
          value += data[index + m_step_x + m_step_z] * dx1 * dz1;
          value += data[index + m_step_x] * dx1 * dz2;
          value += data[index + m_step_z] * dx2 * dz1;
        }
        else if (dx1 != 0 && dy1 != 0 && dz1 == 0)
        {
          // 1 1 0
          value += data[index + m_step_x + m_step_y] * dx1 * dy1;
          value += data[index + m_step_y] * dx2 * dy1;
          value += data[index + m_step_x] * dx1 * dy2;
        }

        return static_cast<float>(value);
      }
    }

    // everything is inside
    double dx2 = 1.0 - dx1, dy2 = 1.0 - dy1, dz2 = 1.0 - dz1;
    int ind = xx * m_step_x + yy * m_step_y + zz * m_step_z;

    return static_cast<float>(
             data[ind] * dx2 * dy2 * dz2
             + data[ind + m_step_x] * dx1 * dy2 * dz2
             + data[ind + m_step_y] * dx2 * dy1 * dz2
             + data[ind + m_step_z] * dx2 * dy2 * dz1
             + data[ind + m_step_x + m_step_z] * dx1 * dy2 * dz1
             + data[ind + m_step_y + m_step_z] * dx2 * dy1 * dz1
             + data[ind + m_step_x + m_step_y] * dx1 * dy1 * dz2
             + data[ind + m_step_x + m_step_y + m_step_z] * dx1 * dy1 * dz1);
  }
}
