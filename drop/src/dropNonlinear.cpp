#include "dropNonlinear.h"
#include "dropFFD.h"
#include "miaImageProcessing.h"
#include "RandomField.h"
#include "Optimization.h"
#include <tbb/parallel_for.h>
#include <iostream>
#include <cmath>

namespace drop
{
  using namespace mrfopt;

  class CostFunctionNonlinear
  {
  public:

    void evaluate(const FFD &ffd, const std::vector<mia::Image> &field, const std::vector<Eigen::Vector3d> &labels, const Eigen::Vector3d &label_scaling, Eigen::MatrixXd &costs)
    {
      Eigen::Matrix4d transformInImageSpace = m_source.worldToImageTransform() * m_transform * m_target.imageToWorldTransform();

      auto sampler_linear = [&](const Eigen::Vector3d& p, float f)
      {
        return m_source.linear(p, f);
      };
      auto sampler_nearest = [&](const Eigen::Vector3d& p, float f)
      {
        return m_source.nearest(p, f);
      };

      if (m_interpolation == mia::LINEAR)
      {
        switch (m_similarity)
        {
        case MAD:
          evaluate_mad(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_linear);
          break;
        case CC:
          evaluate_cc(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_linear);
          break;
        case ECC:
          evaluate_ecc(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_linear);
          break;
        }
      }
      else
      {
        switch (m_similarity)
        {
        case MAD:
          evaluate_mad(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_nearest);
          break;
        case CC:
          evaluate_cc(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_nearest);
          break;
        case ECC:
          evaluate_ecc(ffd, field, labels, label_scaling, costs, transformInImageSpace, sampler_nearest);
          break;
        }
      }
    }

    void target(const mia::Image &image, const mia::Image &mask)
    {
      m_target = image;
      m_target_mask = mask;
      m_target_min = min(m_target);
      m_target_max = max(m_target);
    }

    void source(const mia::Image &image, const mia::Image &mask)
    {
      m_source = image;
      m_source_mask = mask;
      m_source_min = min(m_source);
      m_source_max = max(m_source);
      m_source_oob = m_source_min - 1.0f;
    }

    void transform(const Eigen::Matrix4d &transform)
    {
      m_transform = transform;
    }

    void similarity(const Similarity &type)
    {
      m_similarity = type;
    }

    void interpolation(const mia::Interpolation &type)
    {
      m_interpolation = type;
    }

  private:

    template <typename Sampler>
    void evaluate_mad(const FFD &ffd, const std::vector<mia::Image> &field, const std::vector<Eigen::Vector3d> &labels, const Eigen::Vector3d &label_scaling, Eigen::MatrixXd &costs, const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      costs.setZero();

      std::vector<Eigen::Vector2d> weights_linear_x = ffd.weightsLinearX();
      std::vector<Eigen::Vector2d> weights_linear_y = ffd.weightsLinearY();
      std::vector<Eigen::Vector2d> weights_linear_z = ffd.weightsLinearZ();

      double scale_x = static_cast<double>(ffd.sizeX() - 1) / static_cast<double>(m_target.sizeX());
      double scale_y = static_cast<double>(ffd.sizeY() - 1) / static_cast<double>(m_target.sizeY());
      double scale_z = static_cast<double>(ffd.sizeZ() - 1) / static_cast<double>(m_target.sizeZ());

      tbb::parallel_for(0, static_cast<int>(labels.size()), [&](int label)
      {
        Eigen::Vector3d disp = labels[label].cwiseProduct(label_scaling);

        std::vector<double> sum_of_weights(ffd.size());
        for (auto &i : sum_of_weights) i = 0;

        for (int z = 0; z < m_target.sizeZ(); z++)
        {
          for (int y = 0; y < m_target.sizeY(); y++)
          {
            for (int x = 0; x < m_target.sizeX(); x++)
            {
              Eigen::Vector4d p = transformInImageSpace * Eigen::Vector4d(x + (field[0](x, y, z) + disp[0]) / m_target.spacing()[0], y + (field[1](x, y, z) + disp[1]) / m_target.spacing()[1], z + (field[2](x, y, z) + disp[2]) / m_target.spacing()[2], 1.0);
              Eigen::Vector3d pixel = p.segment(0, 3);

              if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
              {
                float trg_value = m_target(x, y, z);
                float src_value = sampler(pixel, m_source_oob);
                if (src_value != m_source_oob)
                {
                  auto value = std::abs(trg_value - src_value);

                  auto ind_x = static_cast<int>((x + 0.5) * scale_x);
                  auto ind_y = static_cast<int>((y + 0.5) * scale_y);
                  auto ind_z = static_cast<int>((z + 0.5) * scale_z);

                  // project value on FFD control points
                  for (int dz = 0; dz < 2; dz++)
                  {
                    for (int dy = 0; dy < 2; dy++)
                    {
                      for (int dx = 0; dx < 2; dx++)
                      {
                        auto ptx = ind_x + dx;
                        auto pty = ind_y + dy;
                        auto ptz = ind_z + dz;
                        auto pti = ptx + pty * ffd.sizeX() + ptz * ffd.sizeX() * ffd.sizeY();

                        double weight = weights_linear_x[x][dx] * weights_linear_y[y][dy] * weights_linear_z[z][dz];
                        costs(pti, label) += value * weight;
                        sum_of_weights[pti] += weight;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        for (size_t p = 0; p < ffd.size(); p++)
        {
          if (sum_of_weights[p] > 0.0) costs(p, label) /= sum_of_weights[p];
        }
      });
    }

    template <typename Sampler>
    void evaluate_cc(const FFD &ffd, const std::vector<mia::Image> &field, const std::vector<Eigen::Vector3d> &labels, const Eigen::Vector3d &label_scaling, Eigen::MatrixXd &costs, const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      costs.setZero();

      std::vector<Eigen::Vector2d> weights_linear_x = ffd.weightsLinearX();
      std::vector<Eigen::Vector2d> weights_linear_y = ffd.weightsLinearY();
      std::vector<Eigen::Vector2d> weights_linear_z = ffd.weightsLinearZ();

      double scale_x = static_cast<double>(ffd.sizeX() - 1) / static_cast<double>(m_target.sizeX());
      double scale_y = static_cast<double>(ffd.sizeY() - 1) / static_cast<double>(m_target.sizeY());
      double scale_z = static_cast<double>(ffd.sizeZ() - 1) / static_cast<double>(m_target.sizeZ());

      tbb::parallel_for(0, static_cast<int>(labels.size()), [&](int label)
      {
        Eigen::Vector3d disp = labels[label].cwiseProduct(label_scaling);

        Eigen::MatrixXd cc(ffd.size(), 5);
        cc.setZero();
        std::vector<double> counts(ffd.size());
        for (auto &i : counts) i = 0;

        for (int z = 0; z < m_target.sizeZ(); z++)
        {
          for (int y = 0; y < m_target.sizeY(); y++)
          {
            for (int x = 0; x < m_target.sizeX(); x++)
            {
              Eigen::Vector4d p = transformInImageSpace * Eigen::Vector4d(x + (field[0](x, y, z) + disp[0]) / m_target.spacing()[0], y + (field[1](x, y, z) + disp[1]) / m_target.spacing()[1], z + (field[2](x, y, z) + disp[2]) / m_target.spacing()[2], 1.0);
              Eigen::Vector3d pixel = p.segment(0, 3);

              if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
              {
                float trg_value = m_target(x, y, z);
                float src_value = sampler(pixel, m_source_oob);
                if (src_value != m_source_oob)
                {
                  auto ind_x = static_cast<int>((x + 0.5) * scale_x);
                  auto ind_y = static_cast<int>((y + 0.5) * scale_y);
                  auto ind_z = static_cast<int>((z + 0.5) * scale_z);

                  // project value on FFD control points
                  for (int dz = 0; dz < 2; dz++)
                  {
                    for (int dy = 0; dy < 2; dy++)
                    {
                      for (int dx = 0; dx < 2; dx++)
                      {
                        auto ptx = ind_x + dx;
                        auto pty = ind_y + dy;
                        auto ptz = ind_z + dz;
                        auto pti = ptx + pty * ffd.sizeX() + ptz * ffd.sizeX() * ffd.sizeY();

                        cc(pti, 0) += trg_value * trg_value;
                        cc(pti, 1) += trg_value;
                        cc(pti, 2) += src_value * src_value;
                        cc(pti, 3) += src_value;
                        cc(pti, 4) += trg_value * src_value;
                        counts[pti]++;

                      }
                    }
                  }
                }
              }
            }
          }
        }

        for (size_t p = 0; p < ffd.size(); p++)
        {
          auto trg_std = sqrt(static_cast<double>(counts[p])* cc(p, 0) - cc(p, 1) * cc(p, 1));
          auto src_std = sqrt(static_cast<double>(counts[p])* cc(p, 2) - cc(p, 3) * cc(p, 3));
          auto covar = static_cast<double>(counts[p])* cc(p, 4) - cc(p, 1) * cc(p, 3);
          if (trg_std > 0 && src_std > 0)
          {
            auto cc = covar / trg_std / src_std;
            costs(p, label) = -cc;
          }
        }
      });
    }

    template <typename Sampler>
    void evaluate_ecc(const FFD &ffd, const std::vector<mia::Image> &field, const std::vector<Eigen::Vector3d> &labels, const Eigen::Vector3d &label_scaling, Eigen::MatrixXd &costs, const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      costs.setZero();

      std::vector<Eigen::Vector2d> weights_linear_x = ffd.weightsLinearX();
      std::vector<Eigen::Vector2d> weights_linear_y = ffd.weightsLinearY();
      std::vector<Eigen::Vector2d> weights_linear_z = ffd.weightsLinearZ();

      double scale_x = static_cast<double>(ffd.sizeX() - 1) / static_cast<double>(m_target.sizeX());
      double scale_y = static_cast<double>(ffd.sizeY() - 1) / static_cast<double>(m_target.sizeY());
      double scale_z = static_cast<double>(ffd.sizeZ() - 1) / static_cast<double>(m_target.sizeZ());

      tbb::parallel_for(0, static_cast<int>(labels.size()), [&](int label)
      {
        Eigen::Vector3d disp = labels[label].cwiseProduct(label_scaling);

        std::vector<mia::JointHistogram> hist_joint(ffd.size());
        for (auto &h : hist_joint) h = mia::JointHistogram(16, 16, m_target_min, m_target_max, m_source_min, m_source_max);

        std::vector<mia::Histogram> hist_trg(ffd.size());
        for (auto &h : hist_trg) h = mia::Histogram(16, m_target_min, m_target_max);

        std::vector<mia::Histogram> hist_src(ffd.size());
        for (auto &h : hist_src) h = mia::Histogram(16, m_source_min, m_source_max);

        for (int z = 0; z < m_target.sizeZ(); z++)
        {
          for (int y = 0; y < m_target.sizeY(); y++)
          {
            for (int x = 0; x < m_target.sizeX(); x++)
            {
              Eigen::Vector4d p = transformInImageSpace * Eigen::Vector4d(x + (field[0](x, y, z) + disp[0]) / m_target.spacing()[0], y + (field[1](x, y, z) + disp[1]) / m_target.spacing()[1], z + (field[2](x, y, z) + disp[2]) / m_target.spacing()[2], 1.0);
              Eigen::Vector3d pixel = p.segment(0, 3);

              if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
              {
                float trg_value = m_target(x, y, z);
                float src_value = sampler(pixel, m_source_oob);
                if (src_value != m_source_oob)
                {
                  auto ind_x = static_cast<int>((x + 0.5) * scale_x);
                  auto ind_y = static_cast<int>((y + 0.5) * scale_y);
                  auto ind_z = static_cast<int>((z + 0.5) * scale_z);

                  // project value on FFD control points
                  for (int dz = 0; dz < 2; dz++)
                  {
                    for (int dy = 0; dy < 2; dy++)
                    {
                      for (int dx = 0; dx < 2; dx++)
                      {
                        auto ptx = ind_x + dx;
                        auto pty = ind_y + dy;
                        auto ptz = ind_z + dz;
                        auto pti = ptx + pty * ffd.sizeX() + ptz * ffd.sizeX() * ffd.sizeY();

                        hist_joint[pti].add(trg_value, src_value);
                        hist_trg[pti].add(trg_value);
                        hist_src[pti].add(src_value);

                      }
                    }
                  }
                }
              }
            }
          }
        }

        for (size_t p = 0; p < ffd.size(); p++)
        {
          double entropy_joint = mia::entropy(hist_joint[p].counts(), hist_joint[p].totalCount());
          double entropy_trg = mia::entropy(hist_trg[p].counts(), hist_trg[p].totalCount());
          double entropy_src = mia::entropy(hist_src[p].counts(), hist_src[p].totalCount());
          double ecc = 2 - 2 * (entropy_joint + std::numeric_limits<double>::epsilon()) / (entropy_src + entropy_trg + std::numeric_limits<double>::epsilon());
          costs(p, label) = -ecc;
        }
      });
    }

    mia::Image m_target;
    mia::Image m_source;
    mia::Image m_target_mask;
    mia::Image m_source_mask;
    float m_target_min;
    float m_target_max;
    float m_source_min;
    float m_source_max;
    float m_source_oob;
    Eigen::Matrix4d m_transform;
    Similarity m_similarity;
    mia::Interpolation m_interpolation;
  };

  class FFDEnergyFunction : public EnergyFunction<double, int>
  {
  public:
    FFDEnergyFunction(const Eigen::MatrixXd &unary_costs, const Eigen::Vector3d &ffd_spacing, const std::vector<mia::Image> &displacements, const std::vector<Eigen::Vector3d> &labels, const Eigen::Vector3d &label_scaling, double lambda)
      : m_unary_costs(unary_costs)
      , m_ffd_spacing(ffd_spacing)
      , m_labels(labels)
      , m_label_scaling(label_scaling)
      , m_lambda(lambda)
    {
      // need copy for direct data access with 1D indices (displacements are sub-images and don't allow for linear memory access)
      m_disp = std::vector<mia::Image>(displacements.size());
      for (size_t i = 0; i < m_disp.size(); i++)
      {
        m_disp[i] = displacements[i].clone();
      }
    }

    double unary_potential(int nodeIndex, int label) const override
    {
      return m_unary_costs(nodeIndex, label);
    }

    double clique_potential(const Clique& clique, const std::vector<int>& labels) const override
    {
      auto clique_size = static_cast<int>(clique.size());
      if (clique_size == 2)
      {
        return pairwise_potential(clique.nodes[0], clique.nodes[1], labels[0], labels[1]);
      }
      else if (clique_size == 3)
      {
        return triplet_potential(clique.nodes[0], clique.nodes[1], clique.nodes[2], labels[0], labels[1], labels[2]);
      }
      else
      {
        return 0;
      }
    }

    double pairwise_potential(int nodeA, int nodeB, int labelA, int labelB) const
    {
      if (m_lambda == 0) return 0;

      Eigen::Vector3d currentA(m_disp[0].data()[nodeA], m_disp[1].data()[nodeA], m_disp[2].data()[nodeA]);
      Eigen::Vector3d currentB(m_disp[0].data()[nodeB], m_disp[1].data()[nodeB], m_disp[2].data()[nodeB]);

      Eigen::Vector3d dispA = m_labels[labelA].cwiseProduct(m_label_scaling);
      Eigen::Vector3d dispB = m_labels[labelB].cwiseProduct(m_label_scaling);

      Eigen::Vector3d diff = (currentA + dispA - currentB - dispB).cwiseQuotient(m_ffd_spacing);

      return m_lambda * diff.squaredNorm();
    }

    double triplet_potential(int nodeA, int nodeB, int nodeC, int labelA, int labelB, int labelC) const
    {
      if (m_lambda == 0) return 0;

      Eigen::Vector3d currentA(m_disp[0].data()[nodeA], m_disp[1].data()[nodeA], m_disp[2].data()[nodeA]);
      Eigen::Vector3d currentB(m_disp[0].data()[nodeB], m_disp[1].data()[nodeB], m_disp[2].data()[nodeB]);
      Eigen::Vector3d currentC(m_disp[0].data()[nodeC], m_disp[1].data()[nodeC], m_disp[2].data()[nodeC]);

      Eigen::Vector3d dispA = m_labels[labelA].cwiseProduct(m_label_scaling);
      Eigen::Vector3d dispB = m_labels[labelB].cwiseProduct(m_label_scaling);
      Eigen::Vector3d dispC = m_labels[labelC].cwiseProduct(m_label_scaling);

      Eigen::Vector3d diff = ((currentA + dispA) - 2 * (currentB + dispB) + (currentC + dispC)).cwiseQuotient(m_ffd_spacing.cwiseProduct(m_ffd_spacing));

      return m_lambda * diff.squaredNorm();
    }

  private:
    Eigen::MatrixXd m_unary_costs;
    Eigen::Vector3d m_ffd_spacing;
    std::vector<mia::Image> m_disp;
    std::vector<Eigen::Vector3d> m_labels;
    Eigen::Vector3d m_label_scaling;
    double m_lambda;
  };

  std::vector<Eigen::Vector3d> generate_labels(int steps, bool mode2d)
  {
		if (!mode2d)
		{
			int directions = 14;

			std::vector<Eigen::Vector3d> labels(steps * directions + 1);
			labels[0].setZero();

			double stepsize = 1.0 / static_cast<double>(steps);
			for (int i = 1, j = 1; i < labels.size(); i += directions, j++)
			{
				double length = j * stepsize;
				double diag = length * length / sqrt(3 * length * length);

				labels[i + 0] = Eigen::Vector3d(length, 0, 0);
				labels[i + 1] = Eigen::Vector3d(0, length, 0);
				labels[i + 2] = Eigen::Vector3d(0, 0, length);
				labels[i + 3] = Eigen::Vector3d(-length, 0, 0);
				labels[i + 4] = Eigen::Vector3d(0, -length, 0);
				labels[i + 5] = Eigen::Vector3d(0, 0, -length);

				labels[i + 6] = Eigen::Vector3d(diag, diag, diag);
				labels[i + 7] = Eigen::Vector3d(diag, diag, -diag);
				labels[i + 8] = Eigen::Vector3d(diag, -diag, diag);
				labels[i + 9] = Eigen::Vector3d(diag, -diag, -diag);
				labels[i + 10] = Eigen::Vector3d(-diag, diag, diag);
				labels[i + 11] = Eigen::Vector3d(-diag, diag, -diag);
				labels[i + 12] = Eigen::Vector3d(-diag, -diag, diag);
				labels[i + 13] = Eigen::Vector3d(-diag, -diag, -diag);
			}

			return labels;
		}
		else
		{
			int directions = 8;

			std::vector<Eigen::Vector3d> labels(steps * directions + 1);
			labels[0].setZero();

			double stepsize = 1.0 / static_cast<double>(steps);
			for (int i = 1, j = 1; i < labels.size(); i += directions, j++)
			{
				double length = j * stepsize;
				double diag = length * length / sqrt(3 * length * length);

				labels[i + 0] = Eigen::Vector3d(length, 0, 0);
				labels[i + 1] = Eigen::Vector3d(0, length, 0);
				labels[i + 2] = Eigen::Vector3d(-length, 0, 0);
				labels[i + 3] = Eigen::Vector3d(0, -length, 0);

				labels[i + 4] = Eigen::Vector3d(diag, diag, 0);
				labels[i + 5] = Eigen::Vector3d(diag, -diag, 0);
				labels[i + 6] = Eigen::Vector3d(-diag, diag, 0);
				labels[i + 7] = Eigen::Vector3d(-diag, -diag, 0);
			}

			return labels;
		}
  }

  void fix_boundaries(mia::Image &field, float oob)
  {
    int sizeX = field.sizeX();
    int sizeY = field.sizeY();
    int sizeZ = field.sizeZ();
    if (sizeZ > 1)
    {
      int sz1 = 0;
      while (sz1 < (sizeZ - 1) && field(0, 0, sz1) == oob) sz1++;

      int sz2 = sizeZ - 1;
      while (sz2 > 0 && field(0, 0, sz2) == oob) sz2--;

      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          for (int z = sz1 - 1; z >= 0; z--)
          {
            if (field(x, y, z) == oob)
            {
              field(x, y, z) = field(x, y, z + 1);
            }
          }
          for (int z = sz2 + 1; z < sizeZ; z++)
          {
            if (field(x, y, z) == oob)
            {
              field(x, y, z) = field(x, y, z - 1);
            }
          }
        }
      }
    }
    else
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = 0; x < sizeX; x++)
        {
          if (field(x, y, 0) == oob)
          {
            field(x, y, 0) = 0.0f;
          }
          if (field(x, y, sizeZ - 1) == oob)
          {
            field(x, y, sizeZ - 1) = 0.0f;
          }
        }
      }
    }

    int sy1 = 0;
    while (sy1 < (sizeY - 1) && field(0, sy1, 0) == oob) sy1++;

    int sy2 = sizeY - 1;
    while (sy2 > 0 && field(0, sy2, 0) == oob) sy2--;

    for (int z = 0; z < sizeZ; z++)
    {
      for (int x = 0; x < sizeX; x++)
      {
        for (int y = sy1 - 1; y >= 0; y--)
        {
          if (field(x, y, z) == oob)
          {
            field(x, y, z) = field(x, y + 1, z);
          }
        }
        for (int y = sy2 + 1; y < sizeY; y++)
        {
          if (field(x, y, z) == oob)
          {
            field(x, y, z) = field(x, y - 1, z);
          }
        }
      }
    }

    int sx1 = 0;
    while (sx1 < (sizeX - 1) && field(sx1, 0, 0) == oob) sx1++;

    int sx2 = sizeX - 1;
    while (sx2 > 0 && field(sx2, 0, 0) == oob) sx2--;

    for (int z = 0; z < sizeZ; z++)
    {
      for (int y = 0; y < sizeY; y++)
      {
        for (int x = sx1 - 1; x >= 0; x--)
        {
          if (field(x, y, z) == oob)
          {
            field(x, y, z) = field(x + 1, y, z);
          }
        }
        for (int x = sx2 + 1; x < sizeX; x++)
        {
          if (field(x, y, z) == oob)
          {
            field(x, y, z) = field(x - 1, y, z);
          }
        }
      }
    }
  }

	void nonlinear(const mia::Image &source, const mia::Image &target, const mia::Image &sourceMask, const mia::Image &targetMask, std::vector<mia::Image> &field, const Eigen::Matrix4d &transform, const MRFType &type, const Similarity &similarity, double ffdSpacing, const RegularizationType &regularization, double regularizationWeight, std::vector<Eigen::Vector3d> &levelFactors, int iterationsPerLevel, double sampling, const mia::Interpolation &interpolation, bool pinBoundary, bool mode2d)
  {
    CostFunctionNonlinear cost_function;
    cost_function.transform(transform);
    cost_function.source(source, sourceMask);
    cost_function.similarity(similarity);
    cost_function.interpolation(interpolation);

    auto target_oob = min(target) - 1.0f;
		int min_dim = !mode2d ? 32 : 1;

    bool first_order = type == FIRST_ORDER;
    bool reset_field = regularization == FLUID;

    int num_label_sampling_steps = 3;
    double label_scaling_factor = 0.66;

    std::vector<Eigen::Vector3d> labels = generate_labels(num_label_sampling_steps, mode2d);
    size_t num_labels = labels.size();

    std::vector<mia::Image> initial_field = field;

    std::cout << "MRF Type:       " << (first_order ? "FIRST_ORDER" : "SECOND_ORDER") << std::endl;
    std::cout << "Similarity:     " << (similarity == MAD ? "MAD" : similarity == CC ? "CC" : "ECC") << std::endl;
    std::cout << "Regularization: " << (reset_field ? "FLUID" : "ELASTIC") << std::endl;
    std::cout << "Interpolation:  " << (interpolation == mia::LINEAR ? "LINEAR" : "NEAREST") << std::endl;

    // registration over multiple FFD and image resolutions
    int num_levels = static_cast<int>(levelFactors.size());
    for (int level = 0; level < num_levels; level++)
    {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Level " << (level + 1) << " of " << num_levels << std::endl;
      mia::Image trg = mia::resample(target, target.spacing().cwiseProduct(levelFactors[level]), interpolation, target_oob, min_dim);

      std::cout << "Image " << trg.spacing()[0] << "x" << trg.spacing()[1] << "x" << trg.spacing()[2] << " mm (" << trg.sizeX() << "x" << trg.sizeY() << "x" << trg.sizeZ() << ")" << std::endl;

      mia::Image msk = trg.clone();
      resample(targetMask, msk, mia::NEAREST);

      mia::Image binary = trg.clone();
      threshold(trg, binary, target_oob, target_oob);
      invert_binary(binary, binary);
      mul(msk, binary, msk);

      if (sampling < 1.0)
      {
        random_binary(binary, static_cast<float>(sampling));
        mul(msk, binary, msk);
      }

      cost_function.target(trg, msk);

      // configure FFD
      double ffd_spacing = ffdSpacing / pow(2.0, level);
      FFD ffd(ffd_spacing, trg.sizeX() * trg.spacing()[0], trg.sizeY() * trg.spacing()[1], trg.sizeZ() * trg.spacing()[2]);
      ffd.compute_weights(trg.sizeX(), trg.sizeY(), trg.sizeZ());

      std::cout << "FFD " << ffd_spacing << " mm (" << ffd.sizeX() << "x" << ffd.sizeY() << "x" << ffd.sizeZ() << ")" << std::endl;
      std::cout << "----------------------------------------" << std::endl;
      // initialize field and field update
      std::vector<mia::Image> fld(field.size());
      std::vector<mia::Image> upd(field.size());
      for (size_t i = 0; i < field.size(); i++)
      {
        fld[i] = trg.clone();
        upd[i] = trg.clone();
        resample(initial_field[i], fld[i], mia::LINEAR, std::numeric_limits<float>::lowest());
        fix_boundaries(fld[i], std::numeric_limits<float>::lowest());
      }

      // create MRF with either pairwise or triple potentials
      mrfopt::Graph graph = mrfopt::grid_graph(ffd.sizeX(), ffd.sizeY(), ffd.sizeZ(), first_order, !first_order);

      std::vector<int> labeling(graph.num_nodes());

      Eigen::MatrixXd unary_costs(ffd.size(), num_labels);

      // multiple iterations on one level
      double label_factor = 1;
      double energy_previous = std::numeric_limits<double>::max();
      for (size_t iter = 0; iter < iterationsPerLevel; iter++)
      {
        std::cout << "# " << (iter+1) << " of " << iterationsPerLevel << ": ";

        if (reset_field)
        {
          // reset displacement field yields fluid-type regularization of updates only
          ffd.reset();
        }
        else
        {
          // fast projection of current field to FFD (for approximate regularization costs)
          ffd.project_linear(fld);
        }

        // setup label scaling
        Eigen::Vector3d label_scaling = ffd.spacing() * 0.4 * label_factor;

        // compute unaries
        cost_function.evaluate(ffd, fld, labels, label_scaling, unary_costs);

        // setup MRF energy function
        FFDEnergyFunction energy_function(unary_costs, ffd.spacing(), ffd.displacements(), labels, label_scaling, regularizationWeight);

        // solve MRF
        for (auto &l : labeling) l = 0;

        double energy_before = mrfopt::compute_energy(graph, energy_function, labeling);

        if (energy_before < energy_previous)
        {
          mrfopt::multi_label(graph, energy_function, labeling, static_cast<int>(num_labels), 3, 10, mrfopt::HOCR, false);
          double energy_after = mrfopt::compute_energy(graph, energy_function, labeling);

          // apply labeling
          for (int z = 0; z < ffd.sizeZ(); z++)
          {
            for (int y = 0; y < ffd.sizeY(); y++)
            {
              for (int x = 0; x < ffd.sizeX(); x++)
              {
                int index = x + y * ffd.sizeX() + z * ffd.sizeX() * ffd.sizeY();
                Eigen::Vector3d disp = labels[labeling[index]].cwiseProduct(label_scaling);
                for (size_t i = 0; i < field.size(); i++)
                {
                  // no linear memory access possible, as FFD displacements are sub-images of padded versions
                  ffd.displacements()[i](x, y, z) = static_cast<float>(disp[i]);
                }
              }
            }
          }

          // update displacment field
          for (auto &u : upd) zeros(u);
					ffd.evaluate(upd, pinBoundary);
          compose(fld[0], fld[1], fld[2], upd[0], upd[1], upd[2]);
          for (size_t i = 0; i < field.size(); i++)
          {
            copy(upd[i], fld[i]);
          }

          std::cout << std::to_string(energy_before / ffd.size()) << " --> " << std::to_string(energy_after / ffd.size()) << " | label factor: " << label_factor << std::endl;

          energy_previous = energy_before;
        }
        else
        {
          std::cout << std::to_string(energy_before / ffd.size()) << " --> " << std::to_string(energy_before / ffd.size()) << " | label factor: " << label_factor << std::endl;

          // update label factor
          label_factor *= label_scaling_factor;

          energy_previous = std::numeric_limits<double>::max();
        }
      }

      // set field for next level
      initial_field = fld;
    }

    // generate final field
    std::cout << "resampling field to full resolution...";
    for (size_t i = 0; i < field.size(); i++)
    {
      resample(initial_field[i], field[i], mia::LINEAR, std::numeric_limits<float>::lowest());
      fix_boundaries(field[i], std::numeric_limits<float>::lowest());
    }
    std::cout << "done." << std::endl;
  }
}
