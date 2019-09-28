#pragma once

#include <Eigen/Dense>
#include "miaImage.h"

namespace drop
{
  enum Similarity
  {
    MAD, // mean absolute differences
    CC,  // correlation coefficient
    ECC, // entropy correlation coefficient
  };

  Eigen::Matrix4d align_center_of_mass(const mia::Image &source, const mia::Image &target);
}
