#pragma once

#include <Eigen/Dense>
#include "dropCommon.h"
#include "miaImage.h"

namespace drop
{
  const double DROP_PI = 3.141592653589793;

  enum TransformType
  {
    RIGID,
    SIMILARITY,
    AFFINE,
  };

  enum TransformParameters
  {
    TRANS_1,
    TRANS_2,
    TRANS_3,
    ANGLE_1,
    ANGLE_2,
    ANGLE_3,
    SCALE_1,
    SCALE_2,
    SCALE_3,
    SHEAR_1,
    SHEAR_2,
    SHEAR_3,
  };

  Eigen::Matrix4d parameters_to_matrix(const std::vector<double> &parameters);

  Eigen::Matrix4d linear(const mia::Image &source, const mia::Image &target, const mia::Image &sourceMask, const mia::Image &targetMask, const Eigen::Matrix4d &transform, const TransformType &type, const Similarity &similarity, std::vector<Eigen::Vector3d> &levelFactors, int iterationsPerLevel, double sampling, const mia::Interpolation &interpolation, bool mode2d = false);
}
