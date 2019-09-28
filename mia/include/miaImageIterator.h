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

#include <iterator>

namespace mia
{
  template <class T>
  class ImageIterator
  {
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;

  public:
    ImageIterator()
      : m_data(nullptr)
      , m_x(0)
      , m_y(0)
      , m_z(0)
      , m_size_x(0)
      , m_size_y(0)
      , m_step_x(0)
      , m_step_y(0)
      , m_step_z(0)
    {}

    explicit ImageIterator(Image& image, bool make_begin = true)
      : m_data( make_begin ? image.data() : image.data() + image.sizeZ()*image.stepZ() )
      , m_x(0)
      , m_y(0)
      , m_z(make_begin ? 0 : image.sizeZ())
      , m_size_x(image.sizeX())
      , m_size_y(image.sizeY())
      , m_step_x(image.stepX())
      , m_step_y(image.stepY())
      , m_step_z(image.stepZ())
    {}

    explicit ImageIterator(const Image& image, bool make_begin = true)
      : m_data( make_begin ? const_cast<float*>(image.data()) : const_cast<float*>(image.data()) + image.sizeZ()*image.stepZ() )
      , m_x(0)
      , m_y(0)
      , m_z(make_begin ? 0 : image.sizeZ())
      , m_size_x(image.sizeX())
      , m_size_y(image.sizeY())
      , m_step_x(image.stepX())
      , m_step_y(image.stepY())
      , m_step_z(image.stepZ())
    {}

    ImageIterator(const ImageIterator& i)
      : m_data(i.m_data)
      , m_x(i.m_x)
      , m_y(i.m_y)
      , m_z(i.m_z)
      , m_size_x(i.m_size_x)
      , m_size_y(i.m_size_y)
      , m_step_x(i.m_step_x)
      , m_step_y(i.m_step_y)
      , m_step_z(i.m_step_z)
    {}

    bool operator==( const ImageIterator& other ) const
    {
      return m_data == other.m_data;
    }

    /// Inequality comparison.
    bool operator!=( const ImageIterator& other ) const
    {
      return !(*this == other);
    }

    /// Pre-increment.
    ImageIterator& operator++()
    {
      m_data += m_step_x;
      ++m_x;
      if ( m_x >= m_size_x )
      {
        m_data = m_data - m_step_x*m_size_x + m_step_y;
        m_x = 0;
        ++m_y;

        if ( m_y >= m_size_y )
        {
          m_data = m_data - m_step_y*m_size_y + m_step_z;
          m_y = 0;
          ++m_z;
        }
      }
      return *this;
    }

    /// Post-increment.
    ImageIterator operator++( int )
    {
      ImageIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    /// Pre-decrement.
    ImageIterator& operator--()
    {
      m_data -= m_step_x;
      --m_x;
      if ( m_x < 0 )
      {
        m_data = m_data + m_step_x*m_size_x - m_step_y;
        m_x = m_size_x - 1;
        --m_y;

        if ( m_y < 0 )
        {
          m_data = m_data + m_step_y*m_size_y - m_step_z;
          m_y = m_size_y - 1;
          --m_z;
        }
      }
      return *this;
    }

    /// Post-decrement.
    ImageIterator operator--( int )
    {
      ImageIterator tmp(*this);
      --(*this);
      return tmp;
    }

    /// Dereference operator.
    reference operator*() const
    {
      return *m_data;
    }

    /// Pointer access operator.
    pointer operator->() const
    {
      return m_data;
    }

  private:
    value_type* m_data;

    int m_x;
    int m_y;
    int m_z;

    int m_size_x;
    int m_size_y;

    int m_step_x;
    int m_step_y;
    int m_step_z;
  };

  inline ImageIterator<float> begin(Image& image)
  {
    return ImageIterator<float>(image);
  }

  inline ImageIterator<float> end(Image& image)
  {
    return ImageIterator<float>(image, false);
  }

  inline ImageIterator<const float> begin(const Image& image)
  {
    return ImageIterator<const float>(image);
  }

  inline ImageIterator<const float> end(const Image& image)
  {
    return ImageIterator<const float>(image, false);
  }
}
