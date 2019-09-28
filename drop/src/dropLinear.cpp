#include "dropLinear.h"
#include "dropSimplex.h"
#include "miaImageProcessing.h"
#include <iostream>
#include <cmath>

namespace drop
{
  Eigen::Matrix4d parameters_to_matrix(const std::vector<double> &parameters)
  {
    Eigen::Matrix4d translation;
    translation.setIdentity();

    translation(0, 3) = parameters[TRANS_1];
    translation(1, 3) = parameters[TRANS_2];
    translation(2, 3) = parameters[TRANS_3];

    Eigen::Matrix4d rotZ1;
    rotZ1.setIdentity();

    rotZ1(0, 0) = cos(parameters[ANGLE_1]);
    rotZ1(0, 1) = -sin(parameters[ANGLE_1]);
    rotZ1(1, 0) = sin(parameters[ANGLE_1]);
    rotZ1(1, 1) = cos(parameters[ANGLE_1]);

    Eigen::Matrix4d rotX;
    rotX.setIdentity();

    rotX(1, 1) = cos(parameters[ANGLE_2]);
    rotX(1, 2) = -sin(parameters[ANGLE_2]);
    rotX(2, 1) = sin(parameters[ANGLE_2]);
    rotX(2, 2) = cos(parameters[ANGLE_2]);

    Eigen::Matrix4d rotZ2;
    rotZ2.setIdentity();

    rotZ2(0, 0) = cos(parameters[ANGLE_3]);
    rotZ2(0, 1) = -sin(parameters[ANGLE_3]);
    rotZ2(1, 0) = sin(parameters[ANGLE_3]);
    rotZ2(1, 1) = cos(parameters[ANGLE_3]);

    if (parameters.size() == 6)
    {
      return rotZ1 * rotX * rotZ2 * translation;
    }
    else if (parameters.size() == 7)
    {
      Eigen::Matrix4d scaling;
      scaling.setIdentity();

      scaling(0, 0) = exp(parameters[SCALE_1]);
      scaling(1, 1) = scaling(0, 0);
      scaling(2, 2) = scaling(0, 0);

      return scaling * rotZ1 * rotX * rotZ2 * translation;
    }
    else
    {
      Eigen::Matrix4d scaling;
      scaling.setIdentity();

      scaling(0, 0) = exp(parameters[SCALE_1]);
      scaling(1, 1) = exp(parameters[SCALE_2]);
      scaling(2, 2) = exp(parameters[SCALE_3]);

      Eigen::Matrix4d shear1;
      shear1.setIdentity();

      shear1(1, 1) = cos(parameters[SHEAR_1]);
      shear1(1, 2) = -sin(parameters[SHEAR_1]);
      shear1(2, 1) = sin(parameters[SHEAR_1]);
      shear1(2, 2) = cos(parameters[SHEAR_1]);

      Eigen::Matrix4d shear2;
      shear2.setIdentity();

      shear2(0, 0) = cos(parameters[SHEAR_2]);
      shear2(0, 2) = sin(parameters[SHEAR_2]);
      shear2(2, 0) = -sin(parameters[SHEAR_2]);
      shear2(2, 2) = cos(parameters[SHEAR_2]);

      Eigen::Matrix4d shear3;
      shear3.setIdentity();

      shear3(0, 0) = cos(parameters[SHEAR_3]);
      shear3(0, 1) = -sin(parameters[SHEAR_3]);
      shear3(1, 0) = sin(parameters[SHEAR_3]);
      shear3(1, 1) = cos(parameters[SHEAR_3]);

      Eigen::Matrix4d shear;
      shear = shear1 * shear2 * shear3;

      Eigen::Matrix4d shear_inv;
      shear_inv = shear.transpose();

      return shear * scaling * shear_inv * rotZ1 * rotX * rotZ2 * translation;
    }
  }

  class CostFunctionLinear : public CostFunctionSimplex
  {
  public:

    double evaluate(const std::vector<double> &parameters) override
    {
      Eigen::Matrix4d transform = m_transform * m_centerwarp * parameters_to_matrix(parameters) * m_centerwarp_inv;
      Eigen::Matrix4d transformInImageSpace = m_source.worldToImageTransform() * transform * m_target.imageToWorldTransform();

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
          return evaluate_mad(transformInImageSpace, sampler_linear);
        case CC:
          return evaluate_cc(transformInImageSpace, sampler_linear);
        case ECC:
          return evaluate_ecc(transformInImageSpace, sampler_linear);
        default:
          return 0.0f;
        }
      }
      else
      {
        switch (m_similarity)
        {
        case MAD:
          return evaluate_mad(transformInImageSpace, sampler_nearest);
        case CC:
          return evaluate_cc(transformInImageSpace, sampler_nearest);
        case ECC:
          return evaluate_ecc(transformInImageSpace, sampler_nearest);
        default:
          return 0.0f;
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

    void transforms(const Eigen::Matrix4d &transform, const Eigen::Matrix4d &centerWarp)
    {
      m_transform = transform;
      m_centerwarp = centerWarp;
      m_centerwarp_inv = centerWarp.inverse();
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
    double evaluate_mad(const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      Eigen::Vector3d deltaX = Eigen::Vector3d(transformInImageSpace(0, 0), transformInImageSpace(1, 0), transformInImageSpace(2, 0));
      Eigen::Vector3d deltaY = Eigen::Vector3d(transformInImageSpace(0, 1), transformInImageSpace(1, 1), transformInImageSpace(2, 1));

      double cost = 0;
      size_t counter = 0;
      for (int z = 0; z < m_target.sizeZ(); z++)
      {
        Eigen::Vector4d p_h = transformInImageSpace * Eigen::Vector4d(0, 0, z, 1);
        Eigen::Vector3d p = p_h.segment(0, 3);
        for (int y = 0; y < m_target.sizeY(); y++)
        {
          Eigen::Vector3d pixel = p;
          for (int x = 0; x < m_target.sizeX(); x++)
          {
            if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
            {
              float trg_value = m_target(x, y, z);
              float src_value = sampler(pixel, m_source_oob);
              if (src_value != m_source_oob)
              {
                cost += std::abs(trg_value - src_value);
                counter++;
              }
            }
            pixel += deltaX;
          }
          p += deltaY;
        }
      }

      if (counter > 0) cost /= static_cast<double>(counter);

      return cost;
    }

    template <typename Sampler>
    double evaluate_cc(const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      Eigen::Vector3d deltaX = Eigen::Vector3d(transformInImageSpace(0, 0), transformInImageSpace(1, 0), transformInImageSpace(2, 0));
      Eigen::Vector3d deltaY = Eigen::Vector3d(transformInImageSpace(0, 1), transformInImageSpace(1, 1), transformInImageSpace(2, 1));

      std::vector<double> cc(5);
      for (auto &v : cc) v = 0;
      size_t counter = 0;
      for (int z = 0; z < m_target.sizeZ(); z++)
      {
        Eigen::Vector4d p_h = transformInImageSpace * Eigen::Vector4d(0, 0, z, 1);
        Eigen::Vector3d p = p_h.segment(0, 3);
        for (int y = 0; y < m_target.sizeY(); y++)
        {
          Eigen::Vector3d pixel = p;
          for (int x = 0; x < m_target.sizeX(); x++)
          {
            if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
            {
              float trg_value = m_target(x, y, z);
              float src_value = sampler(pixel, m_source_oob);
              if (src_value != m_source_oob)
              {
                cc[0] += trg_value * trg_value;
                cc[1] += trg_value;
                cc[2] += src_value * src_value;
                cc[3] += src_value;
                cc[4] += trg_value * src_value;
                counter++;
              }
            }
            pixel += deltaX;
          }
          p += deltaY;
        }
      }

      auto trg_std = sqrt(static_cast<double>(counter)* cc[0] - cc[1] * cc[1]);
      auto src_std = sqrt(static_cast<double>(counter)* cc[2] - cc[3] * cc[3]);
      auto covar = static_cast<double>(counter)* cc[4] - cc[1] * cc[3];
      if (trg_std > 0 && src_std > 0)
      {
        auto cc = covar / trg_std / src_std;
        return -cc;
      }
      else
      {
        return 0;
      }
    }

    template <typename Sampler>
    double evaluate_ecc(const Eigen::Matrix4d &transformInImageSpace, Sampler sampler)
    {
      Eigen::Vector3d deltaX = Eigen::Vector3d(transformInImageSpace(0, 0), transformInImageSpace(1, 0), transformInImageSpace(2, 0));
      Eigen::Vector3d deltaY = Eigen::Vector3d(transformInImageSpace(0, 1), transformInImageSpace(1, 1), transformInImageSpace(2, 1));

      mia::JointHistogram hist_joint(32, 32, m_target_min, m_target_max, m_source_min, m_source_max);
      mia::Histogram hist_trg(32, m_target_min, m_target_max);
      mia::Histogram hist_src(32, m_source_min, m_source_max);

      for (int z = 0; z < m_target.sizeZ(); z++)
      {
        Eigen::Vector4d p_h = transformInImageSpace * Eigen::Vector4d(0, 0, z, 1);
        Eigen::Vector3d p = p_h.segment(0, 3);
        for (int y = 0; y < m_target.sizeY(); y++)
        {
          Eigen::Vector3d pixel = p;
          for (int x = 0; x < m_target.sizeX(); x++)
          {
            if (m_target_mask(x, y, z) != 0.0f && m_source_mask.nearest(pixel) != 0.0f)
            {
              float trg_value = m_target(x, y, z);
              float src_value = sampler(pixel, m_source_oob);
              if (src_value != m_source_oob)
              {
                hist_joint.add(trg_value, src_value);
                hist_trg.add(trg_value);
                hist_src.add(src_value);
              }
            }
            pixel += deltaX;
          }
          p += deltaY;
        }
      }

      double entropy_joint = mia::entropy(hist_joint.counts(), hist_joint.totalCount());
      double entropy_trg = mia::entropy(hist_trg.counts(), hist_trg.totalCount());
      double entropy_src = mia::entropy(hist_src.counts(), hist_src.totalCount());
      double ecc = 2 - 2 * (entropy_joint + std::numeric_limits<double>::epsilon()) / (entropy_src + entropy_trg + std::numeric_limits<double>::epsilon());

      return -ecc;
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
    Eigen::Matrix4d m_centerwarp;
    Eigen::Matrix4d m_centerwarp_inv;
    Similarity m_similarity;
    mia::Interpolation m_interpolation;
  };

  Eigen::Matrix4d linear(const mia::Image &source, const mia::Image &target, const mia::Image &sourceMask, const mia::Image &targetMask, const Eigen::Matrix4d &transform, const TransformType &type, const Similarity &similarity, std::vector<Eigen::Vector3d> &levelFactors, int iterationsPerLevel, double sampling, const mia::Interpolation &interpolation, bool mode2d)
  {
    int num_params = type == RIGID ? 6 : type == SIMILARITY ? 7 : 12;

    std::vector<double> step_sizes(12);
 //   step_sizes[TRANS_1] = 20.0;
 //   step_sizes[TRANS_2] = 20.0;
	//step_sizes[TRANS_3] = !mode2d ? 20.0 : 0.0;
  step_sizes[TRANS_1] = target.sizeX() * target.spacing()[0] * 0.1;
  step_sizes[TRANS_2] = target.sizeY() * target.spacing()[1] * 0.1;
  step_sizes[TRANS_3] = !mode2d ? target.sizeZ() * target.spacing()[2] * 0.1 : 0.0;
  step_sizes[ANGLE_1] = DROP_PI / 8.0;
	step_sizes[ANGLE_2] = !mode2d ? DROP_PI / 8.0 : 0.0;
	step_sizes[ANGLE_3] = !mode2d ? DROP_PI / 8.0 : 0.0;
    step_sizes[SCALE_1] = 0.1;
    step_sizes[SCALE_2] = 0.1;
	step_sizes[SCALE_3] = !mode2d ? 0.1 : 0.0;
	step_sizes[SHEAR_1] = !mode2d ? DROP_PI / 8.0 : 0.0;
	step_sizes[SHEAR_2] = !mode2d ? DROP_PI / 8.0 : 0.0;
	step_sizes[SHEAR_3] = DROP_PI / 8.0;

    Eigen::Matrix4d centerWarp = Eigen::Matrix4d::Identity();
    centerWarp.col(3) = transform.inverse() * source.imageCenterToWorld().homogeneous();

    std::shared_ptr<CostFunctionLinear> cost_function(new CostFunctionLinear());
    cost_function->transforms(transform, centerWarp);
    cost_function->source(source, sourceMask);
    cost_function->similarity(similarity);
    cost_function->interpolation(interpolation);

    SimplexOptimizer optimizer(num_params, step_sizes);
    optimizer.max_iterations(iterationsPerLevel);
    optimizer.cost_function(cost_function);

    auto target_oob = min(target) - 1.0f;
	int min_dim = !mode2d ? 32 : 1;

    std::cout << "Transform:     " << (type == RIGID ? "RIGID" : type == SIMILARITY ? "SIMILARITY" : "AFFINE") << std::endl;
    std::cout << "Similarity:    " << (similarity == MAD ? "MAD" : similarity == CC ? "CC" : "ECC") << std::endl;
    std::cout << "Interpolation: " << (interpolation == mia::LINEAR ? "LINEAR" : "NEAREST") << std::endl;

    int num_levels = static_cast<int>(levelFactors.size());
    for (int level = 0; level < num_levels; level++)
    {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Level " << (level + 1) << " of " << num_levels << std::endl;
	  mia::Image trg = mia::resample(target, target.spacing().cwiseProduct(levelFactors[level]), interpolation, target_oob, min_dim);

      std::cout << "Image " << trg.spacing()[0] << "x" << trg.spacing()[1] << "x" << trg.spacing()[2] << " mm (" << trg.sizeX() << "x" << trg.sizeY() << "x" << trg.sizeZ() << ")" << std::endl;

      mia::Image binary = trg.clone();
      threshold(trg, binary, target_oob, target_oob);
      invert_binary(binary, binary);

      mia::Image msk = trg.clone();
      resample(targetMask, msk, mia::NEAREST);
      mul(msk, binary, msk);

      if (sampling < 1.0)
      {
        random_binary(binary, static_cast<float>(sampling));
        mul(msk, binary, msk);
      }

      cost_function->target(trg, msk);

      double energy_before = cost_function->evaluate(optimizer.parameters());
      optimizer.run();
      double energy_after = cost_function->evaluate(optimizer.parameters());
      std::cout << "Energy: " << std::to_string(energy_before) << " --> " << std::to_string(energy_after) << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }

    return transform * centerWarp * parameters_to_matrix(optimizer.parameters()) * centerWarp.inverse();
  }
}
