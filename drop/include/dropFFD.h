#pragma once

#include <Eigen/Dense>
#include "miaImage.h"

namespace drop
{
  class FFD
  {
  public:

    FFD(double spacing, double domainX, double domainY, double domainZ);

    void compute_weights(int dimX, int dimY, int dimZ);

    void evaluate(std::vector<mia::Image> &field, bool pinBoundary);

    void project(const std::vector<mia::Image> &field);

    void project_linear(const std::vector<mia::Image> &field);

    void reset();

    size_t size() const
    {
      return m_pts_x * m_pts_y * m_pts_z;
    }

    int sizeX() const
    {
      return m_pts_x;
    }

    int sizeY() const
    {
      return m_pts_y;
    }

    int sizeZ() const
    {
      return m_pts_z;
    }

    const std::vector<Eigen::Vector4d>& weightsX() const
    {
      return m_weights_x;
    }

    const std::vector<Eigen::Vector4d>& weightsY() const
    {
      return m_weights_y;
    }

    const std::vector<Eigen::Vector4d>& weightsZ() const
    {
      return m_weights_z;
    }

    const std::vector<Eigen::Vector2d>& weightsLinearX() const
    {
      return m_weights_linear_x;
    }

    const std::vector<Eigen::Vector2d>& weightsLinearY() const
    {
      return m_weights_linear_y;
    }

    const std::vector<Eigen::Vector2d>& weightsLinearZ() const
    {
      return m_weights_linear_z;
    }

    const Eigen::Vector3d& spacing() const
    {
      return m_spacing;
    }

    std::vector<mia::Image>& positions()
    {
      return m_pos;
    }

    std::vector<mia::Image>& displacements()
    {
      return m_disp;
    }

  private:

    void compute_weights(int dim, int pts, std::vector<Eigen::Vector4d>& weights) const;

    void compute_weights_linear(int dim, int pts, std::vector<Eigen::Vector2d>& weights) const;

    void extrapolate();

    int m_pts_x;
    int m_pts_y;
    int m_pts_z;

    Eigen::Vector3d m_domain;
    Eigen::Vector3d m_spacing;

    std::vector<mia::Image> m_pos;
    std::vector<mia::Image> m_disp;
    std::vector<mia::Image> m_disp_padded;

    std::vector<Eigen::Vector4d> m_weights_x;
    std::vector<Eigen::Vector4d> m_weights_y;
    std::vector<Eigen::Vector4d> m_weights_z;
    std::vector<Eigen::Vector2d> m_weights_linear_x;
    std::vector<Eigen::Vector2d> m_weights_linear_y;
    std::vector<Eigen::Vector2d> m_weights_linear_z;
  };
}
