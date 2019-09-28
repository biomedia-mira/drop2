#include "dropCommon.h"
#include "miaImageProcessing.h"

namespace drop
{
  Eigen::Matrix4d align_center_of_mass(const mia::Image &source, const mia::Image &target)
  {
    Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
	mia::Image src = mia::resample(source, Eigen::Vector3d(4.0, 4.0, 4.0), mia::LINEAR, 0.0f, 32);
	mia::Image trg = mia::resample(target, Eigen::Vector3d(4.0, 4.0, 4.0), mia::LINEAR, 0.0f, 32);
    transform.col(3) = (mia::center_of_mass(src) - mia::center_of_mass(trg)).homogeneous();
    return transform;
  }
}
