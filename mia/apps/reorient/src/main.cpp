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

int main(int argc, char* argv[])
{
  std::string filename_in;
  std::string filename_out;

  try
  {
    // Declare the supported options.
    po::options_description options("options");
    options.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>(&filename_in), "filename of input image")
    ("output,o", po::value<std::string>(&filename_out), "filename of output image")
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
  std::cout << "reorienting image...";

  auto input = itkio::load(filename_in);

  auto output = mia::reorient(input);

  itkio::save(output, filename_out);

  auto stop = ch::high_resolution_clock::now();
  std::cout << "done. took " << ch::duration_cast< ch::milliseconds >(stop - start).count() << " ms" << std::endl;
}
