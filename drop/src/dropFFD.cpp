#include "dropFFD.h"
#include "miaImageProcessing.h"

namespace drop
{
  FFD::FFD(double spacing, double domainX, double domainY, double domainZ)
  {
    m_domain = Eigen::Vector3d(domainX, domainY, domainZ);
    m_pts_x = std::max(static_cast<int>(round(m_domain[0] / spacing + 1)), 2);
    m_pts_y = std::max(static_cast<int>(round(m_domain[1] / spacing + 1)), 2);
    m_pts_z = std::max(static_cast<int>(round(m_domain[2] / spacing + 1)), 2);
    m_spacing = m_domain.cwiseQuotient(Eigen::Vector3d(m_pts_x - 1, m_pts_y - 1, m_pts_z - 1));

    m_pos = std::vector<mia::Image>(3);
    m_disp = std::vector<mia::Image>(3);
    m_disp_padded = std::vector<mia::Image>(3);
    for (int i = 0; i < 3; i++)
    {
      m_pos[i] = mia::Image(m_pts_x, m_pts_y, m_pts_z);

      m_disp_padded[i] = mia::Image(m_pts_x + 2, m_pts_y + 2, m_pts_z + 2);
      zeros(m_disp_padded[i]);

      m_disp[i] = mia::subimage(m_disp_padded[i], 1, 1, 1, m_pts_x, m_pts_y, m_pts_z);
    }

    for (int z = 0; z < m_pts_z; z++)
    {
      for (int y = 0; y < m_pts_y; y++)
      {
        for (int x = 0; x < m_pts_x; x++)
        {
          m_pos[0](x, y, z) = static_cast<float>(x * m_spacing[0]);
          m_pos[1](x, y, z) = static_cast<float>(y * m_spacing[1]);
          m_pos[2](x, y, z) = static_cast<float>(z * m_spacing[2]);
        }
      }
    }
  }

  void FFD::compute_weights(int dimX, int dimY, int dimZ)
  {
    m_weights_x = std::vector<Eigen::Vector4d>(dimX);
    m_weights_y = std::vector<Eigen::Vector4d>(dimY);
    m_weights_z = std::vector<Eigen::Vector4d>(dimZ);

    compute_weights(dimX, m_pts_x, m_weights_x);
    compute_weights(dimY, m_pts_y, m_weights_y);
    compute_weights(dimZ, m_pts_z, m_weights_z);

    m_weights_linear_x = std::vector<Eigen::Vector2d>(dimX);
    m_weights_linear_y = std::vector<Eigen::Vector2d>(dimY);
    m_weights_linear_z = std::vector<Eigen::Vector2d>(dimZ);

    compute_weights_linear(dimX, m_pts_x, m_weights_linear_x);
    compute_weights_linear(dimY, m_pts_y, m_weights_linear_y);
    compute_weights_linear(dimZ, m_pts_z, m_weights_linear_z);
  }

  void FFD::compute_weights(int dim, int pts, std::vector<Eigen::Vector4d>& weights) const
  {
    double scale = static_cast<double>(pts - 1) / static_cast<double>(dim);

    for (int i = 0; i < dim; i++)
    {
      double pos = (static_cast<double>(i)+0.5) * scale;
      double h = pos - static_cast<int>(pos);
      weights[i][0] = ((1.0 - h) * (1.0 - h) * (1.0 - h)) / 6.0;
      weights[i][1] = (((3.0 * h - 6.0) * h) * h + 4.0) / 6.0;
      weights[i][2] = ((((-3.0 * h + 3.0) * h + 3.0) * h) + 1.0) / 6.0;
      weights[i][3] = (h * h * h) / 6.0;
    }
  }

  void FFD::compute_weights_linear(int dim, int pts, std::vector<Eigen::Vector2d>& weights) const
  {
    double scale = static_cast<double>(pts - 1) / static_cast<double>(dim);

    for (int i = 0; i < dim; i++)
    {
      double pos = (static_cast<double>(i)+0.5) * scale;
      double h = pos - static_cast<int>(pos);
      weights[i][0] = 1.0 - h;
      weights[i][1] = h;
    }
  }

  void FFD::evaluate(std::vector<mia::Image> &field, bool pinBoundary)
  {
    if (!pinBoundary)
		{
			extrapolate();
		}
		else
		{
      fill_boundaries(m_disp[0], 0.0f);
      fill_boundaries(m_disp[1], 0.0f);
      fill_boundaries(m_disp[2], 0.0f);
		}

    auto dim_x = field[0].sizeX();
    auto dim_y = field[0].sizeY();
    auto dim_z = field[0].sizeZ();

    double scale_x = static_cast<double>(m_pts_x - 1) / static_cast<double>(dim_x);
    double scale_y = static_cast<double>(m_pts_y - 1) / static_cast<double>(dim_y);
    double scale_z = static_cast<double>(m_pts_z - 1) / static_cast<double>(dim_z);

    for (int z = 0; z < dim_z; z++)
    {
      for (int y = 0; y < dim_y; y++)
      {
        for (int x = 0; x < dim_x; x++)
        {
          int ind_x = static_cast<int>((x + 0.5) * scale_x);
          int ind_y = static_cast<int>((y + 0.5) * scale_y);
          int ind_z = static_cast<int>((z + 0.5) * scale_z);

          for (int dz = 0; dz < 4; dz++)
          {
            for (int dy = 0; dy < 4; dy++)
            {
              for (int dx = 0; dx < 4; dx++)
              {
                int pt_x = ind_x + dx;
                int pt_y = ind_y + dy;
                int pt_z = ind_z + dz;

                double weight = m_weights_x[x][dx] * m_weights_y[y][dy] * m_weights_z[z][dz];
                field[0](x, y, z) += static_cast<float>(m_disp_padded[0](pt_x, pt_y, pt_z) * weight);
                field[1](x, y, z) += static_cast<float>(m_disp_padded[1](pt_x, pt_y, pt_z) * weight);
                field[2](x, y, z) += static_cast<float>(m_disp_padded[2](pt_x, pt_y, pt_z) * weight);
              }
            }
          }
        }
      }
    }
  }

  void FFD::project(const std::vector<mia::Image> &field)
  {
    auto dim_x = field[0].sizeX();
    auto dim_y = field[0].sizeY();
    auto dim_z = field[0].sizeZ();

    double scale_x = static_cast<double>(m_pts_x - 1) / static_cast<double>(dim_x);
    double scale_y = static_cast<double>(m_pts_y - 1) / static_cast<double>(dim_y);
    double scale_z = static_cast<double>(m_pts_z - 1) / static_cast<double>(dim_z);

    zeros(m_disp_padded[0]);
    zeros(m_disp_padded[1]);
    zeros(m_disp_padded[2]);
    mia::Image W = m_disp_padded[0].clone();

    for (int z = 0; z < dim_z; z++)
    {
      for (int y = 0; y < dim_y; y++)
      {
        for (int x = 0; x < dim_x; x++)
        {
          int ind_x = static_cast<int>((x + 0.5) * scale_x);
          int ind_y = static_cast<int>((y + 0.5) * scale_y);
          int ind_z = static_cast<int>((z + 0.5) * scale_z);

          for (int dz = 0; dz < 4; dz++)
          {
            for (int dy = 0; dy < 4; dy++)
            {
              for (int dx = 0; dx < 4; dx++)
              {
                int pt_x = ind_x + dx;
                int pt_y = ind_y + dy;
                int pt_z = ind_z + dz;

                double weight = m_weights_x[x][dx] * m_weights_y[y][dy] * m_weights_z[z][dz];

                m_disp_padded[0](pt_x, pt_y, pt_z) += static_cast<float>(field[0](x, y, z) * weight);
                m_disp_padded[1](pt_x, pt_y, pt_z) += static_cast<float>(field[1](x, y, z) * weight);
                m_disp_padded[2](pt_x, pt_y, pt_z) += static_cast<float>(field[2](x, y, z) * weight);
                W(pt_x, pt_y, pt_z) += static_cast<float>(weight);
              }
            }
          }
        }
      }
    }

    div(m_disp_padded[0], W, m_disp_padded[0]);
    div(m_disp_padded[1], W, m_disp_padded[1]);
    div(m_disp_padded[2], W, m_disp_padded[2]);

  }

  void FFD::project_linear(const std::vector<mia::Image> &field)
  {
    auto dim_x = field[0].sizeX();
    auto dim_y = field[0].sizeY();
    auto dim_z = field[0].sizeZ();

    double scale_x = static_cast<double>(m_pts_x - 1) / static_cast<double>(dim_x);
    double scale_y = static_cast<double>(m_pts_y - 1) / static_cast<double>(dim_y);
    double scale_z = static_cast<double>(m_pts_z - 1) / static_cast<double>(dim_z);

    zeros(m_disp[0]);
    zeros(m_disp[1]);
    zeros(m_disp[2]);
    mia::Image W = m_disp[0].clone();

    for (int z = 0; z < dim_z; z++)
    {
      for (int y = 0; y < dim_y; y++)
      {
        for (int x = 0; x < dim_x; x++)
        {
          int ind_x = static_cast<int>((x + 0.5) * scale_x);
          int ind_y = static_cast<int>((y + 0.5) * scale_y);
          int ind_z = static_cast<int>((z + 0.5) * scale_z);

          for (int dz = 0; dz < 2; dz++)
          {
            for (int dy = 0; dy < 2; dy++)
            {
              for (int dx = 0; dx < 2; dx++)
              {
                int pt_x = ind_x + dx;
                int pt_y = ind_y + dy;
                int pt_z = ind_z + dz;

                double weight = m_weights_linear_x[x][dx] * m_weights_linear_y[y][dy] * m_weights_linear_z[z][dz];

                m_disp[0](pt_x, pt_y, pt_z) += static_cast<float>(field[0](x, y, z) * weight);
                m_disp[1](pt_x, pt_y, pt_z) += static_cast<float>(field[1](x, y, z) * weight);
                m_disp[2](pt_x, pt_y, pt_z) += static_cast<float>(field[2](x, y, z) * weight);
                W(pt_x, pt_y, pt_z) += static_cast<float>(weight);
              }
            }
          }
        }
      }
    }

    div(m_disp[0], W, m_disp[0]);
    div(m_disp[1], W, m_disp[1]);
    div(m_disp[2], W, m_disp[2]);

  }

  void FFD::reset()
  {
    zeros(m_disp_padded[0]);
    zeros(m_disp_padded[1]);
    zeros(m_disp_padded[2]);
  }

  void FFD::extrapolate()
  {
    int pad_pts_x = m_pts_x + 2;
    int pad_pts_y = m_pts_y + 2;
    int pad_pts_z = m_pts_z + 2;

    for (int x = 1; x < pad_pts_x - 1; x++)
    {
      for (int y = 1; y < pad_pts_y - 1; y++)
      {
        m_disp_padded[0](x, y, 0) = 2 * m_disp_padded[0](x, y, 1) - m_disp_padded[0](x, y, 2);
        m_disp_padded[1](x, y, 0) = 2 * m_disp_padded[1](x, y, 1) - m_disp_padded[1](x, y, 2);
        m_disp_padded[2](x, y, 0) = 2 * m_disp_padded[2](x, y, 1) - m_disp_padded[2](x, y, 2);
        m_disp_padded[0](x, y, pad_pts_z - 1) = 2 * m_disp_padded[0](x, y, pad_pts_z - 2) - m_disp_padded[0](x, y, pad_pts_z - 3);
        m_disp_padded[1](x, y, pad_pts_z - 1) = 2 * m_disp_padded[1](x, y, pad_pts_z - 2) - m_disp_padded[1](x, y, pad_pts_z - 3);
        m_disp_padded[2](x, y, pad_pts_z - 1) = 2 * m_disp_padded[2](x, y, pad_pts_z - 2) - m_disp_padded[2](x, y, pad_pts_z - 3);
      }
    }

    for (int z = 1; z < pad_pts_z - 1; z++)
    {
      for (int x = 1; x < pad_pts_x - 1; x++)
      {
        m_disp_padded[0](x, 0, z) = 2 * m_disp_padded[0](x, 1, z) - m_disp_padded[0](x, 2, z);
        m_disp_padded[1](x, 0, z) = 2 * m_disp_padded[1](x, 1, z) - m_disp_padded[1](x, 2, z);
        m_disp_padded[2](x, 0, z) = 2 * m_disp_padded[2](x, 1, z) - m_disp_padded[2](x, 2, z);
        m_disp_padded[0](x, pad_pts_y - 1, z) = 2 * m_disp_padded[0](x, pad_pts_y - 2, z) - m_disp_padded[0](x, pad_pts_y - 3, z);
        m_disp_padded[1](x, pad_pts_y - 1, z) = 2 * m_disp_padded[1](x, pad_pts_y - 2, z) - m_disp_padded[1](x, pad_pts_y - 3, z);
        m_disp_padded[2](x, pad_pts_y - 1, z) = 2 * m_disp_padded[2](x, pad_pts_y - 2, z) - m_disp_padded[2](x, pad_pts_y - 3, z);
      }
    }

    for (int z = 1; z < pad_pts_z - 1; z++)
    {
      for (int y = 1; y < pad_pts_y - 1; y++)
      {
        m_disp_padded[0](0, y, z) = 2 * m_disp_padded[0](1, y, z) - m_disp_padded[0](2, y, z);
        m_disp_padded[1](0, y, z) = 2 * m_disp_padded[1](1, y, z) - m_disp_padded[1](2, y, z);
        m_disp_padded[2](0, y, z) = 2 * m_disp_padded[2](1, y, z) - m_disp_padded[2](2, y, z);
        m_disp_padded[0](pad_pts_x - 1, y, z) = 2 * m_disp_padded[0](pad_pts_x - 2, y, z) - m_disp_padded[0](pad_pts_x - 3, y, z);
        m_disp_padded[1](pad_pts_x - 1, y, z) = 2 * m_disp_padded[1](pad_pts_x - 2, y, z) - m_disp_padded[1](pad_pts_x - 3, y, z);
        m_disp_padded[2](pad_pts_x - 1, y, z) = 2 * m_disp_padded[2](pad_pts_x - 2, y, z) - m_disp_padded[2](pad_pts_x - 3, y, z);
      }
    }

    for (int x = 1; x < pad_pts_x - 1; x++)
    {
      m_disp_padded[0](x, 0, 0) = ((2 * m_disp_padded[0](x, 1, 0) - m_disp_padded[0](x, 2, 0)) + (2 * m_disp_padded[0](x, 0, 1) - m_disp_padded[0](x, 0, 2))) / 2;
      m_disp_padded[0](x, pad_pts_y - 1, 0) = ((2 * m_disp_padded[0](x, pad_pts_y - 2, 0) - m_disp_padded[0](x, pad_pts_y - 3, 0)) + (2 * m_disp_padded[0](x, pad_pts_y - 1, 1) - m_disp_padded[0](x, pad_pts_y - 1, 2))) / 2;
      m_disp_padded[0](x, 0, pad_pts_z - 1) = ((2 * m_disp_padded[0](x, 1, pad_pts_z - 1) - m_disp_padded[0](x, 2, pad_pts_z - 1)) + (2 * m_disp_padded[0](x, 0, pad_pts_z - 2) - m_disp_padded[0](x, 0, pad_pts_z - 3))) / 2;
      m_disp_padded[0](x, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[0](x, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[0](x, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[0](x, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[0](x, pad_pts_y - 1, pad_pts_z - 3))) / 2;
      m_disp_padded[1](x, 0, 0) = ((2 * m_disp_padded[1](x, 1, 0) - m_disp_padded[1](x, 2, 0)) + (2 * m_disp_padded[1](x, 0, 1) - m_disp_padded[1](x, 0, 2))) / 2;
      m_disp_padded[1](x, pad_pts_y - 1, 0) = ((2 * m_disp_padded[1](x, pad_pts_y - 2, 0) - m_disp_padded[1](x, pad_pts_y - 3, 0)) + (2 * m_disp_padded[1](x, pad_pts_y - 1, 1) - m_disp_padded[1](x, pad_pts_y - 1, 2))) / 2;
      m_disp_padded[1](x, 0, pad_pts_z - 1) = ((2 * m_disp_padded[1](x, 1, pad_pts_z - 1) - m_disp_padded[1](x, 2, pad_pts_z - 1)) + (2 * m_disp_padded[1](x, 0, pad_pts_z - 2) - m_disp_padded[1](x, 0, pad_pts_z - 3))) / 2;
      m_disp_padded[1](x, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[1](x, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[1](x, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[1](x, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[1](x, pad_pts_y - 1, pad_pts_z - 3))) / 2;
      m_disp_padded[2](x, 0, 0) = ((2 * m_disp_padded[2](x, 1, 0) - m_disp_padded[2](x, 2, 0)) + (2 * m_disp_padded[2](x, 0, 1) - m_disp_padded[2](x, 0, 2))) / 2;
      m_disp_padded[2](x, pad_pts_y - 1, 0) = ((2 * m_disp_padded[2](x, pad_pts_y - 2, 0) - m_disp_padded[2](x, pad_pts_y - 3, 0)) + (2 * m_disp_padded[2](x, pad_pts_y - 1, 1) - m_disp_padded[2](x, pad_pts_y - 1, 2))) / 2;
      m_disp_padded[2](x, 0, pad_pts_z - 1) = ((2 * m_disp_padded[2](x, 1, pad_pts_z - 1) - m_disp_padded[2](x, 2, pad_pts_z - 1)) + (2 * m_disp_padded[2](x, 0, pad_pts_z - 2) - m_disp_padded[2](x, 0, pad_pts_z - 3))) / 2;
      m_disp_padded[2](x, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[2](x, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[2](x, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[2](x, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[2](x, pad_pts_y - 1, pad_pts_z - 3))) / 2;
    }

    for (int y = 1; y < pad_pts_y - 1; y++)
    {
      m_disp_padded[0](0, y, 0) = ((2 * m_disp_padded[0](1, y, 0) - m_disp_padded[0](2, y, 0)) + (2 * m_disp_padded[0](0, y, 1) - m_disp_padded[0](0, y, 2))) / 2;
      m_disp_padded[0](pad_pts_x - 1, y, 0) = ((2 * m_disp_padded[0](pad_pts_x - 2, y, 0) - m_disp_padded[0](pad_pts_x - 3, y, 0)) + (2 * m_disp_padded[0](pad_pts_x - 1, y, 1) - m_disp_padded[0](pad_pts_x - 1, y, 2))) / 2;
      m_disp_padded[0](0, y, pad_pts_z - 1) = ((2 * m_disp_padded[0](1, y, pad_pts_z - 1) - m_disp_padded[0](2, y, pad_pts_z - 1)) + (2 * m_disp_padded[0](0, y, pad_pts_z - 2) - m_disp_padded[0](0, y, pad_pts_z - 3))) / 2;
      m_disp_padded[0](pad_pts_x - 1, y, pad_pts_z - 1) = ((2 * m_disp_padded[0](pad_pts_x - 2, y, pad_pts_z - 1) - m_disp_padded[0](pad_pts_x - 3, y, pad_pts_z - 1)) + (2 * m_disp_padded[0](pad_pts_x - 1, y, pad_pts_z - 2) - m_disp_padded[0](pad_pts_x - 1, y, pad_pts_z - 3))) / 2;
      m_disp_padded[1](0, y, 0) = ((2 * m_disp_padded[1](1, y, 0) - m_disp_padded[1](2, y, 0)) + (2 * m_disp_padded[1](0, y, 1) - m_disp_padded[1](0, y, 2))) / 2;
      m_disp_padded[1](pad_pts_x - 1, y, 0) = ((2 * m_disp_padded[1](pad_pts_x - 2, y, 0) - m_disp_padded[1](pad_pts_x - 3, y, 0)) + (2 * m_disp_padded[1](pad_pts_x - 1, y, 1) - m_disp_padded[1](pad_pts_x - 1, y, 2))) / 2;
      m_disp_padded[1](0, y, pad_pts_z - 1) = ((2 * m_disp_padded[1](1, y, pad_pts_z - 1) - m_disp_padded[1](2, y, pad_pts_z - 1)) + (2 * m_disp_padded[1](0, y, pad_pts_z - 2) - m_disp_padded[1](0, y, pad_pts_z - 3))) / 2;
      m_disp_padded[1](pad_pts_x - 1, y, pad_pts_z - 1) = ((2 * m_disp_padded[1](pad_pts_x - 2, y, pad_pts_z - 1) - m_disp_padded[1](pad_pts_x - 3, y, pad_pts_z - 1)) + (2 * m_disp_padded[1](pad_pts_x - 1, y, pad_pts_z - 2) - m_disp_padded[1](pad_pts_x - 1, y, pad_pts_z - 3))) / 2;
      m_disp_padded[2](0, y, 0) = ((2 * m_disp_padded[2](1, y, 0) - m_disp_padded[2](2, y, 0)) + (2 * m_disp_padded[2](0, y, 1) - m_disp_padded[2](0, y, 2))) / 2;
      m_disp_padded[2](pad_pts_x - 1, y, 0) = ((2 * m_disp_padded[2](pad_pts_x - 2, y, 0) - m_disp_padded[2](pad_pts_x - 3, y, 0)) + (2 * m_disp_padded[2](pad_pts_x - 1, y, 1) - m_disp_padded[2](pad_pts_x - 1, y, 2))) / 2;
      m_disp_padded[2](0, y, pad_pts_z - 1) = ((2 * m_disp_padded[2](1, y, pad_pts_z - 1) - m_disp_padded[2](2, y, pad_pts_z - 1)) + (2 * m_disp_padded[2](0, y, pad_pts_z - 2) - m_disp_padded[2](0, y, pad_pts_z - 3))) / 2;
      m_disp_padded[2](pad_pts_x - 1, y, pad_pts_z - 1) = ((2 * m_disp_padded[2](pad_pts_x - 2, y, pad_pts_z - 1) - m_disp_padded[2](pad_pts_x - 3, y, pad_pts_z - 1)) + (2 * m_disp_padded[2](pad_pts_x - 1, y, pad_pts_z - 2) - m_disp_padded[2](pad_pts_x - 1, y, pad_pts_z - 3))) / 2;
    }

    for (int z = 1; z < pad_pts_z - 1; z++)
    {
      m_disp_padded[0](0, 0, z) = ((2 * m_disp_padded[0](1, 0, z) - m_disp_padded[0](2, 0, z)) + (2 * m_disp_padded[0](0, 1, z) - m_disp_padded[0](0, 2, z))) / 2;
      m_disp_padded[0](pad_pts_x - 1, 0, z) = ((2 * m_disp_padded[0](pad_pts_x - 2, 0, z) - m_disp_padded[0](pad_pts_x - 3, 0, z)) + (2 * m_disp_padded[0](pad_pts_x - 1, 1, z) - m_disp_padded[0](pad_pts_x - 1, 2, z))) / 2;
      m_disp_padded[0](0, pad_pts_y - 1, z) = ((2 * m_disp_padded[0](1, pad_pts_y - 1, z) - m_disp_padded[0](2, pad_pts_y - 1, z)) + (2 * m_disp_padded[0](0, pad_pts_y - 2, z) - m_disp_padded[0](0, pad_pts_y - 3, z))) / 2;
      m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, z) = ((2 * m_disp_padded[0](pad_pts_x - 2, pad_pts_y - 1, z) - m_disp_padded[0](pad_pts_x - 3, pad_pts_y - 1, z)) + (2 * m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 2, z) - m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 3, z))) / 2;
      m_disp_padded[1](0, 0, z) = ((2 * m_disp_padded[1](1, 0, z) - m_disp_padded[1](2, 0, z)) + (2 * m_disp_padded[1](0, 1, z) - m_disp_padded[1](0, 2, z))) / 2;
      m_disp_padded[1](pad_pts_x - 1, 0, z) = ((2 * m_disp_padded[1](pad_pts_x - 2, 0, z) - m_disp_padded[1](pad_pts_x - 3, 0, z)) + (2 * m_disp_padded[1](pad_pts_x - 1, 1, z) - m_disp_padded[1](pad_pts_x - 1, 2, z))) / 2;
      m_disp_padded[1](0, pad_pts_y - 1, z) = ((2 * m_disp_padded[1](1, pad_pts_y - 1, z) - m_disp_padded[1](2, pad_pts_y - 1, z)) + (2 * m_disp_padded[1](0, pad_pts_y - 2, z) - m_disp_padded[1](0, pad_pts_y - 3, z))) / 2;
      m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, z) = ((2 * m_disp_padded[1](pad_pts_x - 2, pad_pts_y - 1, z) - m_disp_padded[1](pad_pts_x - 3, pad_pts_y - 1, z)) + (2 * m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 2, z) - m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 3, z))) / 2;
      m_disp_padded[2](0, 0, z) = ((2 * m_disp_padded[2](1, 0, z) - m_disp_padded[2](2, 0, z)) + (2 * m_disp_padded[2](0, 1, z) - m_disp_padded[2](0, 2, z))) / 2;
      m_disp_padded[2](pad_pts_x - 1, 0, z) = ((2 * m_disp_padded[2](pad_pts_x - 2, 0, z) - m_disp_padded[2](pad_pts_x - 3, 0, z)) + (2 * m_disp_padded[2](pad_pts_x - 1, 1, z) - m_disp_padded[2](pad_pts_x - 1, 2, z))) / 2;
      m_disp_padded[2](0, pad_pts_y - 1, z) = ((2 * m_disp_padded[2](1, pad_pts_y - 1, z) - m_disp_padded[2](2, pad_pts_y - 1, z)) + (2 * m_disp_padded[2](0, pad_pts_y - 2, z) - m_disp_padded[2](0, pad_pts_y - 3, z))) / 2;
      m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, z) = ((2 * m_disp_padded[2](pad_pts_x - 2, pad_pts_y - 1, z) - m_disp_padded[2](pad_pts_x - 3, pad_pts_y - 1, z)) + (2 * m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 2, z) - m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 3, z))) / 2;
    }

    m_disp_padded[0](0, 0, 0) = ((2 * m_disp_padded[0](1, 0, 0) - m_disp_padded[0](2, 0, 0)) + (2 * m_disp_padded[0](0, 1, 0) - m_disp_padded[0](0, 2, 0)) + (2 * m_disp_padded[0](0, 0, 1) - m_disp_padded[0](0, 0, 2))) / 3;
    m_disp_padded[0](pad_pts_x - 1, 0, 0) = ((2 * m_disp_padded[0](pad_pts_x - 2, 0, 0) - m_disp_padded[0](pad_pts_x - 3, 0, 0)) + (2 * m_disp_padded[0](pad_pts_x - 1, 1, 0) - m_disp_padded[0](pad_pts_x - 1, 2, 0)) + (2 * m_disp_padded[0](pad_pts_x - 1, 0, 1) - m_disp_padded[0](pad_pts_x - 1, 0, 2))) / 3;
    m_disp_padded[0](0, pad_pts_y - 1, 0) = ((2 * m_disp_padded[0](1, pad_pts_y - 1, 0) - m_disp_padded[0](2, pad_pts_y - 1, 0)) + (2 * m_disp_padded[0](0, pad_pts_y - 2, 0) - m_disp_padded[0](0, pad_pts_y - 3, 0)) + (2 * m_disp_padded[0](0, pad_pts_y - 1, 1) - m_disp_padded[0](0, pad_pts_y - 1, 2))) / 3;
    m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, 0) = ((2 * m_disp_padded[0](pad_pts_x - 2, pad_pts_y - 1, 0) - m_disp_padded[0](pad_pts_x - 3, pad_pts_y - 1, 0)) + (2 * m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 2, 0) - m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 3, 0)) + (2 * m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, 1) - m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, 2))) / 3;
    m_disp_padded[1](0, 0, 0) = ((2 * m_disp_padded[1](1, 0, 0) - m_disp_padded[1](2, 0, 0)) + (2 * m_disp_padded[1](0, 1, 0) - m_disp_padded[1](0, 2, 0)) + (2 * m_disp_padded[1](0, 0, 1) - m_disp_padded[1](0, 0, 2))) / 3;
    m_disp_padded[1](pad_pts_x - 1, 0, 0) = ((2 * m_disp_padded[1](pad_pts_x - 2, 0, 0) - m_disp_padded[1](pad_pts_x - 3, 0, 0)) + (2 * m_disp_padded[1](pad_pts_x - 1, 1, 0) - m_disp_padded[1](pad_pts_x - 1, 2, 0)) + (2 * m_disp_padded[1](pad_pts_x - 1, 0, 1) - m_disp_padded[1](pad_pts_x - 1, 0, 2))) / 3;
    m_disp_padded[1](0, pad_pts_y - 1, 0) = ((2 * m_disp_padded[1](1, pad_pts_y - 1, 0) - m_disp_padded[1](2, pad_pts_y - 1, 0)) + (2 * m_disp_padded[1](0, pad_pts_y - 2, 0) - m_disp_padded[1](0, pad_pts_y - 3, 0)) + (2 * m_disp_padded[1](0, pad_pts_y - 1, 1) - m_disp_padded[1](0, pad_pts_y - 1, 2))) / 3;
    m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, 0) = ((2 * m_disp_padded[1](pad_pts_x - 2, pad_pts_y - 1, 0) - m_disp_padded[1](pad_pts_x - 3, pad_pts_y - 1, 0)) + (2 * m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 2, 0) - m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 3, 0)) + (2 * m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, 1) - m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, 2))) / 3;
    m_disp_padded[2](0, 0, 0) = ((2 * m_disp_padded[2](1, 0, 0) - m_disp_padded[2](2, 0, 0)) + (2 * m_disp_padded[2](0, 1, 0) - m_disp_padded[2](0, 2, 0)) + (2 * m_disp_padded[2](0, 0, 1) - m_disp_padded[2](0, 0, 2))) / 3;
    m_disp_padded[2](pad_pts_x - 1, 0, 0) = ((2 * m_disp_padded[2](pad_pts_x - 2, 0, 0) - m_disp_padded[2](pad_pts_x - 3, 0, 0)) + (2 * m_disp_padded[2](pad_pts_x - 1, 1, 0) - m_disp_padded[2](pad_pts_x - 1, 2, 0)) + (2 * m_disp_padded[2](pad_pts_x - 1, 0, 1) - m_disp_padded[2](pad_pts_x - 1, 0, 2))) / 3;
    m_disp_padded[2](0, pad_pts_y - 1, 0) = ((2 * m_disp_padded[2](1, pad_pts_y - 1, 0) - m_disp_padded[2](2, pad_pts_y - 1, 0)) + (2 * m_disp_padded[2](0, pad_pts_y - 2, 0) - m_disp_padded[2](0, pad_pts_y - 3, 0)) + (2 * m_disp_padded[2](0, pad_pts_y - 1, 1) - m_disp_padded[2](0, pad_pts_y - 1, 2))) / 3;
    m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, 0) = ((2 * m_disp_padded[2](pad_pts_x - 2, pad_pts_y - 1, 0) - m_disp_padded[2](pad_pts_x - 3, pad_pts_y - 1, 0)) + (2 * m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 2, 0) - m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 3, 0)) + (2 * m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, 1) - m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, 2))) / 3;

    m_disp_padded[0](0, 0, pad_pts_z - 1) = ((2 * m_disp_padded[0](1, 0, pad_pts_z - 1) - m_disp_padded[0](2, 0, pad_pts_z - 1)) + (2 * m_disp_padded[0](0, 1, pad_pts_z - 1) - m_disp_padded[0](0, 2, pad_pts_z - 1)) + (2 * m_disp_padded[0](0, 0, pad_pts_z - 2) - m_disp_padded[0](0, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[0](pad_pts_x - 1, 0, pad_pts_z - 1) = ((2 * m_disp_padded[0](pad_pts_x - 2, 0, pad_pts_z - 1) - m_disp_padded[0](pad_pts_x - 3, 0, pad_pts_z - 1)) + (2 * m_disp_padded[0](pad_pts_x - 1, 1, pad_pts_z - 1) - m_disp_padded[0](pad_pts_x - 1, 2, pad_pts_z - 1)) + (2 * m_disp_padded[0](pad_pts_x - 1, 0, pad_pts_z - 2) - m_disp_padded[0](pad_pts_x - 1, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[0](0, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[0](1, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[0](2, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[0](0, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[0](0, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[0](0, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[0](0, pad_pts_y - 1, pad_pts_z - 3))) / 3;
    m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[0](pad_pts_x - 2, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[0](pad_pts_x - 3, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[0](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 3))) / 3;
    m_disp_padded[1](0, 0, pad_pts_z - 1) = ((2 * m_disp_padded[1](1, 0, pad_pts_z - 1) - m_disp_padded[1](2, 0, pad_pts_z - 1)) + (2 * m_disp_padded[1](0, 1, pad_pts_z - 1) - m_disp_padded[1](0, 2, pad_pts_z - 1)) + (2 * m_disp_padded[1](0, 0, pad_pts_z - 2) - m_disp_padded[1](0, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[1](pad_pts_x - 1, 0, pad_pts_z - 1) = ((2 * m_disp_padded[1](pad_pts_x - 2, 0, pad_pts_z - 1) - m_disp_padded[1](pad_pts_x - 3, 0, pad_pts_z - 1)) + (2 * m_disp_padded[1](pad_pts_x - 1, 1, pad_pts_z - 1) - m_disp_padded[1](pad_pts_x - 1, 2, pad_pts_z - 1)) + (2 * m_disp_padded[1](pad_pts_x - 1, 0, pad_pts_z - 2) - m_disp_padded[1](pad_pts_x - 1, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[1](0, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[1](1, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[1](2, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[1](0, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[1](0, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[1](0, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[1](0, pad_pts_y - 1, pad_pts_z - 3))) / 3;
    m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[1](pad_pts_x - 2, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[1](pad_pts_x - 3, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[1](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 3))) / 3;
    m_disp_padded[2](0, 0, pad_pts_z - 1) = ((2 * m_disp_padded[2](1, 0, pad_pts_z - 1) - m_disp_padded[2](2, 0, pad_pts_z - 1)) + (2 * m_disp_padded[2](0, 1, pad_pts_z - 1) - m_disp_padded[2](0, 2, pad_pts_z - 1)) + (2 * m_disp_padded[2](0, 0, pad_pts_z - 2) - m_disp_padded[2](0, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[2](pad_pts_x - 1, 0, pad_pts_z - 1) = ((2 * m_disp_padded[2](pad_pts_x - 2, 0, pad_pts_z - 1) - m_disp_padded[2](pad_pts_x - 3, 0, pad_pts_z - 1)) + (2 * m_disp_padded[2](pad_pts_x - 1, 1, pad_pts_z - 1) - m_disp_padded[2](pad_pts_x - 1, 2, pad_pts_z - 1)) + (2 * m_disp_padded[2](pad_pts_x - 1, 0, pad_pts_z - 2) - m_disp_padded[2](pad_pts_x - 1, 0, pad_pts_z - 3))) / 3;
    m_disp_padded[2](0, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[2](1, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[2](2, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[2](0, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[2](0, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[2](0, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[2](0, pad_pts_y - 1, pad_pts_z - 3))) / 3;
    m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 1) = ((2 * m_disp_padded[2](pad_pts_x - 2, pad_pts_y - 1, pad_pts_z - 1) - m_disp_padded[2](pad_pts_x - 3, pad_pts_y - 1, pad_pts_z - 1)) + (2 * m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 2, pad_pts_z - 1) - m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 3, pad_pts_z - 1)) + (2 * m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 2) - m_disp_padded[2](pad_pts_x - 1, pad_pts_y - 1, pad_pts_z - 3))) / 3;
  }
}
