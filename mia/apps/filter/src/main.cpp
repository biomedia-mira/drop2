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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

using namespace mia;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template <typename T>
void strings_to_values(const std::vector<std::string>& string_seq, std::vector<T>& values)
{
  for (std::vector<std::string>::const_iterator it = string_seq.begin(); it != string_seq.end(); ++it)
  {
    std::stringstream ss(*it);
    std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), back_inserter(values));
  }
}

void params_to_vector(const std::vector<double>& params, Eigen::Vector3d& vector)
{
  if (params.size() > 0)
  {
    vector[0] = (params.size() > 0) ? params[0] : 0;
    vector[1] = (params.size() > 1) ? params[1] : params[0];
    vector[2] = (params.size() > 2) ? params[2] : params[0];
  }
  else
  {
    vector.setZero();
  }
}

int main(int argc, char* argv[])
{
  std::string filename_in;
  std::string filename_out;
  std::vector<std::string> skernel;

  try
  {
    // Declare the supported options.
    po::options_description options("options");
    options.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>(&filename_in), "filename of input image")
    ("output,o", po::value<std::string>(&filename_out), "filename of output image")
    ("gauss", po::value<std::vector<std::string>>(&skernel)->multitoken(), "kernel size (in mm)")
    ;

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help") || vm.size() == 0)
    {
      std::cout << options << std::endl;
      return 0;
    }
  }
  catch (std::exception& e)
  {
    std::cout << e.what() << std::endl;
    return 1;
  }

  namespace ch = std::chrono;
  auto start = ch::high_resolution_clock::now();
  std::cout << "filtering image...";

  auto input = itkio::load(filename_in);

  std::vector<double> kernel;
  strings_to_values(skernel, kernel);

  Image output = input.clone();
  output.dataType(mia::FLOAT);
  zeros(output);

  gauss(input, output, kernel[0], kernel[1], kernel[2]);

  itkio::save(output, filename_out);

  auto stop = ch::high_resolution_clock::now();
  std::cout << "done. took " << ch::duration_cast< ch::milliseconds >(stop - start).count() << " ms" << std::endl;
}
