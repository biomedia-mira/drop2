/*
* mia - a lightweight C++ image processing library.
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

#include "itkio.h"
#include "miaImage.h"
#include "miaImageProcessing.h"

#include <iostream>
#include <vector>
#include <chrono>

using namespace mia;

int main(int argc, char* argv[])
{

  Image image1(64, 64, 64);
  checkerboard(image1, 4, 4, 4);

  namespace ch = std::chrono;
  auto start = ch::high_resolution_clock::now();

  Image image2 = resample(image1, Eigen::Vector3d(2.0, 2.0, 2.0), Interpolation::LINEAR);
  Image image3 = subimage(image1, image1.sizeX() / 4, image1.sizeY() / 4, image1.sizeZ() / 4, image1.sizeX() / 2, image1.sizeY() / 2, image1.sizeZ() / 2);

  gauss(image3, image3, 1.0, 1.0, 1.0);

  itkio::save(image1, "test1.nii.gz");
  itkio::save(image2, "test2.nii.gz");
  itkio::save(image3, "test3.nii.gz");

  auto stop = ch::high_resolution_clock::now();
  std::cout << "duration: " << ch::duration_cast< ch::milliseconds >(stop-start).count() << " ms" << std::endl;
}
