#include "itkio.h"
#include "miaImage.h"
#include "miaImageProcessing.h"
#include "dropLinear.h"
#include "dropFFD.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

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

std::vector<Eigen::Vector3d> parse_level_factors(const std::vector<std::string>& string_seq)
{
  std::vector<Eigen::Vector3d> level_factors;
  std::vector<double> level_values;
  strings_to_values(string_seq, level_values);
  for (int l = 0; l < level_values.size() / 3; l++)
  {
    Eigen::Vector3d factors(level_values[l * 3 + 0], level_values[l * 3 + 1], level_values[l * 3 + 2]);
    level_factors.push_back(factors);
  }

  return level_factors;
}

int main(int argc, char* argv[])
{
  std::string filename_input;
  std::string filename_output;

  int interpolation;
  float fill_value;

  float max_translation;
  float max_rotation;
  float max_scaling;
  float max_shearing;
  float max_disp;
  double spacing;
  bool pin_boundary;
  bool mode_2d;
  bool no_matrix;
  bool no_field;

  try
  {
    // Declare the supported options.
    po::options_description options("options");
    options.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>(&filename_input), "filename of input image")
    ("output,o", po::value<std::string>(&filename_output), "filename of output image")
    ("mode2d", po::bool_switch(&mode_2d)->default_value(false), "enable 2D mode")        
    ("max_trans", po::value<float>(&max_translation)->default_value(25.0f), "maximum translation")
    ("max_rot", po::value<float>(&max_rotation)->default_value(10.0f), "maximum rotation")
    ("max_scale", po::value<float>(&max_scaling)->default_value(0.2f), "maximum scaling")
    ("max_shear", po::value<float>(&max_shearing)->default_value(10.0f), "maximum shearing")    
    ("max_disp", po::value<float>(&max_disp)->default_value(10.0f), "maximum displacement")
    ("spacing", po::value<double>(&spacing)->default_value(10), "FFD spacing")
    ("pin", po::bool_switch(&pin_boundary)->default_value(false), "pin boundaries")
    ("interp", po::value<int>(&interpolation)->default_value(1), "image interpolation (0=NEAREST, 1=LINEAR)")
    ("fill", po::value<float>(&fill_value)->default_value(0.0f), "fill value for out of bounds")
    ("no_matrix", po::bool_switch(&no_matrix)->default_value(false), "no matrix output")
    ("no_field", po::bool_switch(&no_field)->default_value(false), "no field output")
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

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<float> sampler(-1.0f, 1.0f);

  // load image
  mia::Image image = itkio::load(filename_input);

  // initialize parameters
  std::vector<double> parameters(12);
  parameters[drop::TRANS_1] = sampler(mt) * max_translation;
  parameters[drop::TRANS_2] = sampler(mt) * max_translation;
  parameters[drop::TRANS_3] = !mode_2d ? sampler(mt) * max_translation : 0.0f;
  parameters[drop::ANGLE_1] = sampler(mt) * max_rotation / 180.0f * drop::DROP_PI;
  parameters[drop::ANGLE_2] = !mode_2d ? sampler(mt) * max_rotation / 180.0f * drop::DROP_PI : 0.0f;
  parameters[drop::ANGLE_3] = !mode_2d ? sampler(mt) * max_rotation / 180.0f * drop::DROP_PI : 0.0f;
  parameters[drop::SCALE_1] = sampler(mt) * max_scaling;
  parameters[drop::SCALE_2] = sampler(mt) * max_scaling;
  parameters[drop::SCALE_3] = !mode_2d ? sampler(mt) * max_scaling : 0.0f;
  parameters[drop::SHEAR_1] = !mode_2d ? sampler(mt) * max_shearing / 180.0f * drop::DROP_PI : 0.0f;
  parameters[drop::SHEAR_2] = !mode_2d ? sampler(mt) * max_shearing / 180.0f * drop::DROP_PI : 0.0f;
  parameters[drop::SHEAR_3] = sampler(mt) * max_shearing / 180.0f * drop::DROP_PI;

  Eigen::Matrix4d centerWarp = Eigen::Matrix4d::Identity();
  centerWarp.col(3) = image.imageCenterToWorld().homogeneous();

  Eigen::Matrix4d transform = centerWarp * drop::parameters_to_matrix(parameters) * centerWarp.inverse();

  // initialize field
  std::vector<mia::Image> field(3);
  field[0] = image.clone();
  zeros(field[0]);
  field[1] = image.clone();
  zeros(field[1]);
  field[2] = image.clone();
  zeros(field[2]);
  for (auto &f : field)
  {
    f.dataType(mia::FLOAT);
  }

  drop::FFD ffd(spacing, image.sizeX() * image.spacing()[0], image.sizeY() * image.spacing()[1], image.sizeZ() * image.spacing()[2]);
  ffd.compute_weights(image.sizeX(), image.sizeY(), image.sizeZ());

  // assign random displacement
  for (int z = 0; z < ffd.sizeZ(); z++)
  {
    for (int y = 0; y < ffd.sizeY(); y++)
    {
      for (int x = 0; x < ffd.sizeX(); x++)
      {
        ffd.displacements()[0](x, y, z) = sampler(mt) * max_disp;
        ffd.displacements()[1](x, y, z) = sampler(mt) * max_disp;
        if (!mode_2d)
        {
          ffd.displacements()[2](x, y, z) = sampler(mt) * max_disp;
        }          
      }
    }
  }

  // update displacment field
  ffd.evaluate(field, pin_boundary);

  // image warping and saving results
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "WARPING image and saving results..." << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  warp(image.clone(), image, transform, field[0], field[1], field[2], (mia::Interpolation)interpolation, fill_value);

  fs::path output_path(filename_output);
  std::string basename = fs::basename(output_path);
  if (fs::extension(basename) != "") basename = fs::basename(basename);

  std::string out_folder = output_path.parent_path().string();
  if (out_folder == "") out_folder = ".";
  else if (!fs::exists(out_folder)) fs::create_directories(out_folder);

  itkio::save(image, filename_output);

  if (!no_matrix)
  {
    std::stringstream filename_transform;
    filename_transform << out_folder << "/" << basename << "_transform.txt";

    std::ofstream ofs(filename_transform.str());
    for (int r = 0; r < 4; r++)
    {
      ofs << transform(r, 0);
      for (int c = 1; c < 4; c++)
      {
        ofs << "\t" << transform(r, c);
      }
      ofs << std::endl;
    }
    ofs.close();
  }

  if (!no_field)
  {
    std::stringstream filename_field_x;
    filename_field_x << out_folder << "/" << basename << "_field_x.nii.gz";
    itkio::save(field[0], filename_field_x.str());

    std::stringstream filename_field_y;
    filename_field_y << out_folder << "/" << basename << "_field_y.nii.gz";
    itkio::save(field[1], filename_field_y.str());

    std::stringstream filename_field_z;
    filename_field_z << out_folder << "/" << basename << "_field_z.nii.gz";
    itkio::save(field[2], filename_field_z.str());
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Finished. Timing (s): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() / 1000.0 << std::endl;
}
