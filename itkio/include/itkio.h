/*
* itkio - a cross-platform C++ library for simplified IO of image data using the ITK library.
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

#include "itkImage.h"
#include "miaImage.h"

namespace itkio
{
  enum ComponentType
  {
    UNKNOWNCOMPONENTTYPE,
    UCHAR,
    CHAR,
    USHORT,
    SHORT,
    UINT,
    INT,
    ULONG,
    LONG,
    FLOAT,
    DOUBLE
  };

  /**
   * \brief Saves an image to a file with a particular data type.
   * \param image The image to be saved.
   * \param filename The filename.
   **/
  template <typename PixelType>
  void save(const mia::Image& image, const std::string& filename);

  /**
   * \brief Saves an image to a file.
   * \param image The image to be saved.
   * \param filename The filename.
   **/
  void save(const mia::Image& image, const std::string& filename);

  /**
   * \brief Loads an image from a file.
   * \param filename The filename.
   * \return The image.
   **/
  mia::Image load(const std::string& filename);

  /**
   * \brief Loads an image from a DICOM sequence corresponding to a given DICOM file.
   * \param filename The filename.
   * \return The image.
   **/
  mia::Image load_dicom(const std::string& filename);

  /**
   * \brief Gets all DICOM sequence UIDs found in the directory of a given DICOM file.
   * \param filename The filename.
   * \return Array of UIDs.
   **/
  std::vector<std::string> get_uids(const std::string& filename);

  /**
  * \brief Converts a mia::Image to itk::Image.
  * \param image The mia::Image.
  * \return itk::Image.
  **/
  itk::Image<float, 3>::Pointer miaToItk(const mia::Image& image);
}
