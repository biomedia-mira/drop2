#pragma once

#include <Eigen/Dense>
#include "dropCommon.h"
#include "miaImage.h"

namespace drop
{
  enum MRFType
  {
    FIRST_ORDER,
    SECOND_ORDER,
  };

  enum RegularizationType
  {
    FLUID,
    ELASTIC,
  };

  void nonlinear(const mia::Image &source, const mia::Image &target, const mia::Image &sourceMask, const mia::Image &targetMask, std::vector<mia::Image> &field, const Eigen::Matrix4d &transform, const MRFType &type, const Similarity &similarity, double ffdSpacing, const RegularizationType &regularization, double regularizationWeight, std::vector<Eigen::Vector3d> &levelFactors, int iterationsPerLevel, double sampling, const mia::Interpolation &interpolation, bool pinBoundary = false, bool mode2d = false);
}
