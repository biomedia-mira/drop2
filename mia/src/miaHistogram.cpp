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

#include "miaHistogram.h"

namespace mia
{
  Histogram::Histogram()
    : m_total_count(0)
    , m_min_value(0.0f)
    , m_max_value(0.0f)
    , m_counts()
  {}

  Histogram::Histogram(int numBins, float minValue, float maxValue)
    : m_total_count(0)
    , m_min_value(minValue)
    , m_max_value(maxValue)
    , m_counts(numBins)
  {
    reset();
  }

  Histogram::Histogram(const Histogram& other)
    : m_total_count(other.m_total_count)
    , m_min_value(other.m_min_value)
    , m_max_value(other.m_max_value)
    , m_counts(other.m_counts)
  {}

  Histogram::Histogram(Histogram&& other)
    : m_total_count(std::move(other.m_total_count))
    , m_min_value(std::move(other.m_min_value))
    , m_max_value(std::move(other.m_max_value))
    , m_counts(std::move(other.m_counts))
  {}

  Histogram& Histogram::operator=(const Histogram& other)
  {
    if (this != &other)
    {
      *this = Histogram(other);
    }
    return *this;
  }

  Histogram& Histogram::operator=(Histogram&& other)
  {
    //if (this != &other)
    //{
    //  *this = Histogram(other);
    //}
    m_total_count = std::move(other.m_total_count);
    m_min_value = std::move(other.m_min_value);
    m_max_value = std::move(other.m_max_value);
    m_counts = std::move(other.m_counts);
    return *this;
  }

  void Histogram::add(float value)
  {
    if ((value - m_min_value) >= 0)
    {
      //should be equivalent to Matlab's hist function
      int bin = static_cast<int>((value - m_min_value) / (m_max_value - m_min_value + 1.0f) * static_cast<float>(m_counts.size()));

      //maximum values are included in last bin
      if (bin == m_counts.size()) bin--;

      if (bin < m_counts.size())
      {
        m_counts[bin]++;
        m_total_count++;
      }
    }
  }

  float Histogram::probability(float value) const
  {
    if ((value - m_min_value) < 0) return -1.0f;

    int bin = static_cast<int>((value - m_min_value) / (m_max_value - m_min_value) * m_counts.size());

    //maximum values are included in last bin
    if (bin == m_counts.size()) bin--;

    if (bin < m_counts.size())
      return static_cast<float>(m_counts[bin]) / static_cast<float>(m_total_count);
    else
      return -1.0f;
  }

  JointHistogram::JointHistogram()
    : m_total_count(0)
    , m_min_value_x(0.0f)
    , m_max_value_x(0.0f)
    , m_min_value_y(0.0f)
    , m_max_value_y(0.0f)
    , m_counts()
  {}

  JointHistogram::JointHistogram(int numBinsX, int numBinsY, float minValueX, float maxValueX, float minValueY, float maxValueY)
    : m_total_count(0)
    , m_min_value_x(minValueX)
    , m_max_value_x(maxValueX)
    , m_min_value_y(minValueY)
    , m_max_value_y(maxValueY)
    , m_counts(numBinsX, numBinsY)
  {
    reset();
  }

  JointHistogram::JointHistogram(const JointHistogram& other)
    : m_total_count(other.m_total_count)
    , m_min_value_x(other.m_min_value_x)
    , m_max_value_x(other.m_max_value_x)
    , m_min_value_y(other.m_min_value_y)
    , m_max_value_y(other.m_max_value_y)
    , m_counts(other.m_counts)
  {}

  JointHistogram::JointHistogram(JointHistogram&& other)
    : m_total_count(std::move(other.m_total_count))
    , m_min_value_x(std::move(other.m_min_value_x))
    , m_max_value_x(std::move(other.m_max_value_x))
    , m_min_value_y(std::move(other.m_min_value_y))
    , m_max_value_y(std::move(other.m_max_value_y))
    , m_counts(std::move(other.m_counts))
  {}

  JointHistogram& JointHistogram::operator=(const JointHistogram& other)
  {
    if (this != &other)
    {
      *this = JointHistogram(other);
    }
    return *this;
  }

  JointHistogram& JointHistogram::operator=(JointHistogram&& other)
  {
    m_total_count = std::move(other.m_total_count);
    m_min_value_x = std::move(other.m_min_value_x);
    m_max_value_x = std::move(other.m_max_value_x);
    m_min_value_y = std::move(other.m_min_value_y);
    m_max_value_y = std::move(other.m_max_value_y);
    m_counts = std::move(other.m_counts);
    return *this;
  }

  void JointHistogram::add(float valueX, float valueY)
  {
    if ((valueX - m_min_value_x) >= 0 && (valueY - m_min_value_y) >= 0)
    {
      //should be equivalent to Matlab's hist function
      int bin_x = static_cast<int>((valueX - m_min_value_x) / (m_max_value_x - m_min_value_x + 1.0f) * static_cast<float>(m_counts.rows()));
      int bin_y = static_cast<int>((valueY - m_min_value_y) / (m_max_value_y - m_min_value_y + 1.0f) * static_cast<float>(m_counts.cols()));

      //maximum values are included in last bin
      if (bin_x == m_counts.rows()) bin_x--;
      if (bin_y == m_counts.cols()) bin_y--;

      if (bin_x < m_counts.rows() && bin_y < m_counts.cols())
      {
        m_counts(bin_x, bin_y)++;
        m_total_count++;
      }
    }
  }

  float JointHistogram::probability(float valueX, float valueY) const
  {
    if ((valueX - m_min_value_x) < 0 && (valueY - m_min_value_y) < 0) return -1.0f;

    int bin_x = static_cast<int>((valueX - m_min_value_x) / (m_max_value_x - m_min_value_x) * m_counts.rows());
    int bin_y = static_cast<int>((valueY - m_min_value_y) / (m_max_value_y - m_min_value_y) * m_counts.cols());

    //maximum values are included in last bin
    if (bin_x == m_counts.rows()) bin_x--;
    if (bin_y == m_counts.cols()) bin_y--;

    if (bin_x < m_counts.rows() && bin_y < m_counts.cols())
      return static_cast<float>(m_counts(bin_x, bin_y)) / static_cast<float>(m_total_count);
    else
      return -1.0f;
  }

  int maxCountIndex(const Eigen::ArrayXd& counts)
  {
    int index = 0;
    double count = counts[0];
    for (int i = 1; i < counts.size(); i++)
    {
      if (counts[i] > count)
      {
        count = counts[i];
        index = i;
      }
    }

    if (count > 0)
      return index;
    else
      return -1;
  }

  double entropy(const Eigen::ArrayXd& counts, size_t totalCount)
  {
    double entropy = 0;
    for (int i = 0; i < counts.size(); i++)
    {
      double prob = counts[i] / static_cast<double>(totalCount);
      if (prob > 0) entropy -= prob * log(prob);
    }
    return entropy;
  }

  double entropy(const Eigen::MatrixXd& counts, size_t totalCount)
  {
    double entropy = 0;
    for (int i = 0; i < counts.size(); i++)
    {
      double prob = counts(i) / static_cast<double>(totalCount);
      if (prob > 0) entropy -= prob * log(prob);
    }
    return entropy;
  }
}
