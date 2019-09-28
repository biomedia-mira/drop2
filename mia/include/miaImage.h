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

#include <vector>
#include <memory>
#include <Eigen/Dense>

namespace mia
{
  enum Interpolation
  {
    NEAREST,
    LINEAR,
  };

  enum ImageDataType
  {
    BYTE,
    SHORT,
    USHORT,
    INT,
    UINT,
    FLOAT,
    DOUBLE,
  };

  enum ImagePlane
  {
    XY,
    XZ,
    YZ,
  };

  class Image
  {
  public:
    Image();

    Image(int sizeX, int sizeY, int sizeZ);

    Image(float* data, int sizeX, int sizeY, int sizeZ);

    Image(float* data, int sizeX, int sizeY, int sizeZ, int stepX, int stepY, int stepZ);

    Image(const Image& other);

    Image(Image&& other);

    Image& operator=(const Image& other);

    Image& operator=(Image&& other);

    /**
     * \brief Returns a deep copy of the image.
     * \return The copy.
     **/
    Image clone() const;

    /**
     * \brief Returns a pointer to the image data.
     * \return The data pointer.
     **/
    float* data()
    {
      return m_data.get();
    }

    /**
     * \brief Returns a const pointer to the image data.
     * \return The data pointer.
     **/
    const float* data() const
    {
      return m_data.get();
    }

    /**
     * \brief Returns a reference to the image data at specified image coordinate.
     * \param x Image x-coordinate.
     * \param y Image y-coordinate.
     * \param z Image z-coordinate.
     * \return The data reference.
     **/
    float& operator()(int x, int y, int z)
    {
      return const_cast<float&>(static_cast<const Image&>(*this)(x,y,z));
    }

    /**
     * \brief Returns a const reference to the image data at specified image coordinate.
     * \param x Image x-coordinate.
     * \param y Image y-coordinate.
     * \param z Image z-coordinate.
     * \return The data reference.
     **/
    const float& operator()(int x, int y, int z) const
    {
      return *(m_data.get() + x*m_step_x + y*m_step_y + z*m_step_z);
    }

    /**
     * \brief Sets the external data type of the image (used when writing the image to a file).
     * \param dataType The data type.
     **/
    void dataType(const ImageDataType& dataType)
    {
      m_datatype = dataType;
    }

    /**
     * \brief Gets the external data type of the image (used when writing the image to a file).
     * \return The data type.
     **/
    const ImageDataType& dataType() const
    {
      return m_datatype;
    }

    /**
     * \brief Returns the size of the image, i.e. number of image points.
     * \return The image size.
     **/
    size_t size() const
    {
      return static_cast<size_t>(m_size_x) * static_cast<size_t>(m_size_y) * static_cast<size_t>(m_size_z);
    }

    /**
     * \brief Returns the size of the image along the x-axis.
     * \return The size along x-axis.
     **/
    int sizeX() const
    {
      return m_size_x;
    }

    /**
     * \brief Returns the size of the image along the y-axis.
     * \return The size along y-axis.
     **/
    int sizeY() const
    {
      return m_size_y;
    }

    /**
     * \brief Returns the size of the image along the z-axis.
     * \return The size along z-axis.
     **/
    int sizeZ() const
    {
      return m_size_z;
    }

    /**
     * \brief Returns the step size between data points along the x-axis.
     * \return The step size along x-axis.
     **/
    int stepX() const
    {
      return m_step_x;
    }

    /**
     * \brief Returns the step size between data points along the y-axis.
     * \return The step size along y-axis.
     **/
    int stepY() const
    {
      return m_step_y;
    }

    /**
     * \brief Returns the step size between data points along the z-axis.
     * \return The step size along z-axis.
     **/
    int stepZ() const
    {
      return m_step_z;
    }

    /**
     * \brief Checks whether the image data is located in a contiguous block.
     * \return True if data is contiguous.
     **/
    bool contiguous() const
    {
      return m_step_x==1 && m_step_y==m_size_x && m_step_z==m_size_x*m_size_y;
    }

    /**
     * \brief Sets the element spacing.
     * \param spacing A 3-vector containing the element spacing.
     **/
    void spacing(const Eigen::Vector3d& spacing)
    {
      m_spacing = spacing;
    }

    /**
     * \brief Gets the element spacing.
     * \return A 3-vector containing the element spacing.
     **/
    const Eigen::Vector3d& spacing() const
    {
      return m_spacing;
    }

    /**
     * \brief Sets the image origin in world coordinates.
     * \param origin A 3-vector containing the image origin.
     **/
    void origin(const Eigen::Vector3d& origin)
    {
      m_origin = origin;
    }

    /**
     * \brief Gets the image origin in world coordinates.
     * \return A 3-vector containing the image origin.
     **/
    const Eigen::Vector3d& origin() const
    {
      return m_origin;
    }

    /**
     * \brief Sets the image to world orientation.
     * \param imageToWorld A 3x3-matrix containing the image to world orientation.
     **/
    void imageToWorld(const Eigen::Matrix3d& imageToWorld);

    /**
     * \brief Gets the image to world orientation.
     * \return A 3x3-matrix containing the image to world orientation.
     **/
    const Eigen::Matrix3d& imageToWorld() const;

    /**
     * \brief Sets the world to image orientation.
     * \param worldToImage A 3x3-matrix containing the world to image orientation.
     **/
    void worldToImage(const Eigen::Matrix3d& worldToImage);

    /**
     * \brief Gets the world to image orientation.
     * \return A 3x3-matrix containing the world to image orientation.
     **/
    const Eigen::Matrix3d& worldToImage() const;

    /**
     * \brief Transforms a point from image to world coordinates.
     * \param point Image coordinates.
     * \return World coordinates.
     **/
    Eigen::Vector3d mapImageToWorld(const Eigen::Vector3d& point) const;

    /**
     * \brief Transforms a point from world to image coordinates.
     * \param point World coordinates.
     * \return Image coordinates.
     **/
    Eigen::Vector3d mapWorldToImage(const Eigen::Vector3d& point) const;

    /**
     * \brief Returns the image center in world coordinates.
     * \return World coordinates.
     **/
    Eigen::Vector3d imageCenterToWorld() const;

    /**
     * \brief Returns the full image to world transformation.
     * \return Image to world transformation.
     **/
    Eigen::Matrix4d imageToWorldTransform() const;

    /**
     * \brief Returns the full world to image transformation.
     * \return World to image transformation.
     **/
    Eigen::Matrix4d worldToImageTransform() const;

    /**
     * \brief Nearest neighbor interpolation at specified image coordinate.
     *\ param point Image coordinates.
     *\ param outsideValue Intensity value used for out of bounds access.
     * \return The interpolated intensity value.
     **/
    float nearest(const Eigen::Vector3d& point, float outsideValue = 0.0f) const;

    /**
     * \brief Nearest neighbor interpolation at specified image coordinate.
     *\ param point Image coordinates.
     *\ param outsideValue Intensity value used for out of bounds access.
     * \return The interpolated intensity value.
     **/
    float linear(const Eigen::Vector3d& point, float outsideValue = 0.0f) const;

  private:
    std::shared_ptr<float> m_data; // raw image data

    ImageDataType m_datatype; // image data type

    int m_size_x; // image size along x
    int m_size_y; // image size along y
    int m_size_z; // image size along z

    int m_step_x; // the number of floats between two elements along x
    int m_step_y; // the number of floats between two elements along y
    int m_step_z; // the number of floats between two elements along z

    Eigen::Vector3d m_spacing; // element spacing
    Eigen::Vector3d m_origin; // image origin in world coordinates
    Eigen::Matrix3d m_image_to_world; // image-to-world orientation
    Eigen::Matrix3d m_world_to_image; // world-to-image orientation
  };
}
