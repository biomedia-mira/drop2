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
#include <Eigen/Core>

namespace mia
{
  class Histogram
  {
  public:
    Histogram();

    Histogram(int numBins, float minValue, float maxValue);

    Histogram(const Histogram& other);

    Histogram(Histogram&& other);

    Histogram& operator=(const Histogram& other);

    Histogram& operator=(Histogram&& other);

    Histogram& operator+=(const Histogram& rhs)
    {
      if (m_counts.size() == rhs.m_counts.size())
      {
        m_total_count += rhs.m_total_count;
        m_counts += rhs.m_counts;
      }
      return *this;
    }

    Histogram operator+(const Histogram& rhs) const
    {
      Histogram res(*this);
      res += rhs;
      return res;
    }

    Histogram& operator-=(const Histogram& rhs)
    {
      if (m_counts.size() == rhs.m_counts.size())
      {
        m_total_count -= rhs.m_total_count;
        m_counts -= rhs.m_counts;
      }
      return *this;
    }

    Histogram operator-(const Histogram& rhs) const
    {
      Histogram res(*this);
      res -= rhs;
      return res;
    }

    /**
     * \brief Sets all histogram counts to zero.
     **/
    void reset()
    {
      m_counts.setZero();
    }

    /**
     * \brief Returns the sum of all counts.
     * \return The total count.
     **/
    size_t totalCount() const
    {
      return m_total_count;
    }

    /**
     * \brief Adds a new value to the histogram by increasing the count of the corresponding bin.
     * \param value The value to be added.
     **/
    void add(float value);

    /**
     * \brief Computes the empirical probability for a given value based on the histogram counts.
     * \param value The value.
     * \return The probability.
     **/
    float probability(float value) const;

    /**
     * \brief Returns the number of bins.
     * \return The number of bins.
     **/
    int numBins() const
    {
      return static_cast<int>(m_counts.size());
    }

    /**
     * \brief Returns the bin array with counts.
     * \return The bin array with counts.
     **/
    const Eigen::ArrayXd& counts()
    {
      return m_counts;
    }

    /**
     * \brief Returns the minimum value considered for the histogram.
     * \return The minimum value.
     **/
    float minValue() const
    {
      return m_min_value;
    }

    /**
     * \brief Returns the maximum value considered for the histogram.
     * \return The maximum value.
     **/
    float maxValue() const
    {
      return m_max_value;
    }

  private:
    size_t m_total_count;

    float m_min_value;

    float m_max_value;

    Eigen::ArrayXd m_counts;
  };

  class JointHistogram
  {
  public:
    JointHistogram();

    JointHistogram(int numBinsX, int numBinsY, float minValueX, float maxValueX, float minValueY, float maxValueY);

    JointHistogram(const JointHistogram& other);

    JointHistogram(JointHistogram&& other);

    JointHistogram& operator=(const JointHistogram& other);

    JointHistogram& operator=(JointHistogram&& other);

    JointHistogram& operator+=(const JointHistogram& rhs)
    {
      if (m_counts.size() == rhs.m_counts.size())
      {
        m_total_count += rhs.m_total_count;
        m_counts += rhs.m_counts;
      }
      return *this;
    }

    JointHistogram operator+(const JointHistogram& rhs) const
    {
      JointHistogram res(*this);
      res += rhs;
      return res;
    }

    JointHistogram& operator-=(const JointHistogram& rhs)
    {
      if (m_counts.size() == rhs.m_counts.size())
      {
        m_total_count -= rhs.m_total_count;
        m_counts -= rhs.m_counts;
      }
      return *this;
    }

    JointHistogram operator-(const JointHistogram& rhs) const
    {
      JointHistogram res(*this);
      res -= rhs;
      return res;
    }

    /**
    * \brief Sets all histogram counts to zero.
    **/
    void reset()
    {
      m_counts.setZero();
    }

    /**
    * \brief Returns the sum of all counts.
    * \return The total count.
    **/
    size_t totalCount() const
    {
      return m_total_count;
    }

    /**
    * \brief Adds a new pair of values to the histogram by increasing the count of the corresponding bin.
    * \param valueX The value for X.
    * \param valueY The value for Y.
    **/
    void add(float valueX, float valueY);

    /**
    * \brief Computes the empirical probability for a given pair of values based on the histogram counts.
    * \param valueX The value for X.
    * \param valueY The value for Y.
    * \return The probability.
    **/
    float probability(float valueX, float valueY) const;

    /**
    * \brief Returns the number of bins for variable X.
    * \return The number of bins.
    **/
    int numBinsX() const
    {
      return static_cast<int>(m_counts.rows());
    }

    /**
    * \brief Returns the number of bins for variable Y.
    * \return The number of bins.
    **/
    int numBinsY() const
    {
      return static_cast<int>(m_counts.cols());
    }

    /**
    * \brief Returns the bin matrix with counts.
    * \return The bin matrix with counts.
    **/
    const Eigen::MatrixXd& counts()
    {
      return m_counts;
    }

    /**
    * \brief Returns the minimum value considered for the histogram for variable X.
    * \return The minimum value.
    **/
    float minValueX() const
    {
      return m_min_value_x;
    }

    /**
    * \brief Returns the maximum value considered for the histogram for variable X.
    * \return The maximum value.
    **/
    float maxValueX() const
    {
      return m_max_value_x;
    }

    /**
    * \brief Returns the minimum value considered for the histogram for variable Y.
    * \return The minimum value.
    **/
    float minValueY() const
    {
      return m_min_value_y;
    }

    /**
    * \brief Returns the maximum value considered for the histogram for variable Y.
    * \return The maximum value.
    **/
    float maxValueY() const
    {
      return m_max_value_y;
    }

  private:
    size_t m_total_count;

    float m_min_value_x;

    float m_max_value_x;

    float m_min_value_y;

    float m_max_value_y;

    Eigen::MatrixXd m_counts;
  };

  class IntegralHistogram
  {
  public:
    IntegralHistogram()
      : m_size_x(0)
      , m_size_y(0)
      , m_size_z(0)
      , m_histograms(0)
    {}

    IntegralHistogram(int sizeX, int sizeY, int sizeZ)
      : m_size_x(sizeX)
      , m_size_y(sizeY)
      , m_size_z(sizeZ)
      , m_histograms(sizeX*sizeY*sizeZ)
    {}

    Histogram& operator()(int x, int y, int z)
    {
      return const_cast<Histogram&>(static_cast<const IntegralHistogram&>(*this)(x,y,z));
    }
    const Histogram& operator()(int x, int y, int z) const
    {
      return m_histograms[x+y*m_size_x+z*m_size_x*m_size_y];
    }

  private:
    int m_size_x;
    int m_size_y;
    int m_size_z;

    std::vector<Histogram> m_histograms;
  };

  int maxCountIndex(const Eigen::ArrayXd& counts);

  double entropy(const Eigen::ArrayXd& counts, size_t totalCount);

  double entropy(const Eigen::MatrixXd& counts, size_t totalCount);
}
