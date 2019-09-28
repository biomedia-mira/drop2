#include "itkio.h"
#include "miaImage.h"
#include "miaImageProcessing.h"
#include "dropFFD.h"
#include "dropLinear.h"
#include "dropNonlinear.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

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
  std::string filename_source;
  std::string filename_target;
  std::string filename_output;

  std::string filename_source_mask;
  std::string filename_target_mask;

  std::string filename_transform;
  std::string filename_field_x;
  std::string filename_field_y;
  std::string filename_field_z;

  std::vector<std::string> llevels;
  std::vector<std::string> nlevels;

  int o_interpolation;
  float o_fill_value;

  bool run_com_alignment;

  bool run_linear;
  int l_type;
  int l_similarity;
  int l_iterations;
  int l_interpolation;
  double l_sampling;

  bool run_nonlinear;
  int n_type;
  int n_reg;
  double n_ffd;
  int n_similarity;
  int n_iterations;
  int n_interpolation;
  double n_sampling;
  double n_lambda;
  bool n_pin_boundary;

  bool mode_2d;

  try
  {
    // Declare the supported options.
    po::options_description options("options");
    options.add_options()
    ("help,h", "produce help message")
    ("source,s", po::value<std::string>(&filename_source), "filename of source image")
    ("target,t", po::value<std::string>(&filename_target), "filename of target image")
    ("output,o", po::value<std::string>(&filename_output), "filename of output image")
    ("smask", po::value<std::string>(&filename_source_mask)->default_value(""), "filename of source mask")
    ("tmask", po::value<std::string>(&filename_target_mask)->default_value(""), "filename of target mask")
    ("transform", po::value<std::string>(&filename_transform)->default_value(""), "filename of transformation")
    ("fx", po::value<std::string>(&filename_field_x)->default_value(""), "filename of displacement field for X-component")
    ("fy", po::value<std::string>(&filename_field_y)->default_value(""), "filename of displacement field for Y-component")
    ("fz", po::value<std::string>(&filename_field_z)->default_value(""), "filename of displacement field for Z-component")
    ("ointerp", po::value<int>(&o_interpolation)->default_value(1), "OUTPUT: image interpolation (0=NEAREST, 1=LINEAR)")
	  ("ofill", po::value<float>(&o_fill_value)->default_value(0.0f), "OUTPUT: fill value for out of bounds")
    ("mode2d", po::bool_switch(&mode_2d)->default_value(false), "enable 2D mode")    
    ("com,c", po::bool_switch(&run_com_alignment)->default_value(false), "run CENTER OF MASS alignment")
    ("linear,l", po::bool_switch(&run_linear)->default_value(false), "run LINEAR registration")
    ("ltype", po::value<int>(&l_type)->default_value(0), "LINEAR: transformation type (0=RIGID, 1=SIMILARITY, 2=AFFINE)")
    ("lsim", po::value<int>(&l_similarity)->default_value(0), "LINEAR: similarity measure (0=MAD, 1=CC, 2=ECC)")
    ("llevels", po::value<std::vector<std::string>>(&llevels)->multitoken(), "LINEAR: image level factors")
    ("liters", po::value<int>(&l_iterations)->default_value(100), "LINEAR: number of iterations per level")
    ("linterp", po::value<int>(&l_interpolation)->default_value(1), "LINEAR: image interpolation (0=NEAREST, 1=LINEAR)")
    ("lsampling", po::value<double>(&l_sampling)->default_value(0.05), "LINEAR: image sampling probability [0,1]")
    ("nonlinear,n", po::bool_switch(&run_nonlinear)->default_value(false), "run NONLINEAR registration")
    ("ntype", po::value<int>(&n_type)->default_value(0), "NONLINEAR: MRF type (0=FIRST_ORDER, 1=SECOND_ORDER)")
    ("nreg", po::value<int>(&n_reg)->default_value(0), "NONLINEAR: regularization type (0=FLUID, 1=ELASTIC)")
    ("nffd", po::value<double>(&n_ffd)->default_value(80), "NONLINEAR: FFD spacing on first level")
    ("nsim", po::value<int>(&n_similarity)->default_value(0), "NONLINEAR: similarity measure (0=MAD, 1=CC, 2=ECC)")
    ("nlevels", po::value<std::vector<std::string>>(&nlevels)->multitoken(), "NONLINEAR: image level factors")
    ("niters", po::value<int>(&n_iterations)->default_value(10), "NONLINEAR: number of iterations per level")
    ("ninterp", po::value<int>(&n_interpolation)->default_value(1), "NONLINEAR: image interpolation (0=NEAREST, 1=LINEAR)")
    ("nsampling", po::value<double>(&n_sampling)->default_value(1), "NONLINEAR: image sampling probability [0,1]")
    ("nlambda", po::value<double>(&n_lambda)->default_value(0), "NONLINEAR: regularization weight")
    ("npin", po::bool_switch(&n_pin_boundary)->default_value(false), "NONLINEAR: pin boundaries")
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

  // load images
  mia::Image source = itkio::load(filename_source);
  mia::Image target = itkio::load(filename_target);

  // load image masks
  mia::Image source_mask;
  if (filename_source_mask != "")
  {
    source_mask = itkio::load(filename_source_mask);
  }
  else
  {
    source_mask = source.clone();
    ones(source_mask);
  }

  mia::Image target_mask;
  if (filename_target_mask != "")
  {
    target_mask = itkio::load(filename_target_mask);
  }
  else
  {
    target_mask = target.clone();
    ones(target_mask);
  }

  // load transformation
  Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
  if (filename_transform != "")
  {
    std::ifstream ifs(filename_transform);
    for (int r = 0; r < 4; r++)
    {
      std::string line;
      getline(ifs, line);
      std::vector<std::string> str;
      boost::split(str, line, boost::is_any_of("\t"), boost::token_compress_on);
      if (str.size() == 4)
      {
        for (int c = 0; c < 4; c++)
        {
          transform(r, c) = atof(str[c].c_str());
        }
      }
    }
    ifs.close();
  }

  // load displacement field
  std::vector<mia::Image> field(3);
  if (filename_field_x != "")
  {
    field[0] = itkio::load(filename_field_x);
  }
  else
  {
    field[0] = target.clone();
    zeros(field[0]);
  }
  if (filename_field_y != "")
  {
    field[1] = itkio::load(filename_field_y);
  }
  else
  {
    field[1] = target.clone();
    zeros(field[1]);
  }
  if (filename_field_z != "")
  {
    field[2] = itkio::load(filename_field_z);
  }
  else
  {
    field[2] = target.clone();
    zeros(field[2]);
  }
  for (auto &f : field)
  {
    f.dataType(mia::FLOAT);
  }

  // center of mass alignment
  if (run_com_alignment)
  {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Running CENTER OF MASS alignment..." << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    transform = drop::align_center_of_mass(source, target);

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Finished. Timing (s): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() / 1000.0 << std::endl;
    std::cout << std::endl;
  }

  // linear registration
  if (run_linear)
  {
    std::vector<Eigen::Vector3d> l_level_factors = parse_level_factors(llevels);

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Running LINEAR registration..." << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    transform = drop::linear(source, target, source_mask, target_mask, transform, (drop::TransformType)l_type, (drop::Similarity)l_similarity, l_level_factors, l_iterations, l_sampling, (mia::Interpolation)l_interpolation, mode_2d);

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Finished. Timing (s): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() / 1000.0 << std::endl;
    std::cout << std::endl;
  }

  // nonlinear registration
  if (run_nonlinear)
  {
    std::vector<Eigen::Vector3d> n_level_factors = parse_level_factors(nlevels);

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Running NONLINEAR registration..." << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    drop::nonlinear(source, target, source_mask, target_mask, field, transform, (drop::MRFType)n_type, (drop::Similarity)n_similarity, n_ffd, (drop::RegularizationType)n_reg, n_lambda, n_level_factors, n_iterations, n_sampling, (mia::Interpolation)n_interpolation, n_pin_boundary, mode_2d);

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Finished. Timing (s): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop - start).count() / 1000.0 << std::endl;
    std::cout << std::endl;
  }

  // image warping and saving results
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "WARPING image and saving results..." << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  warp(source, target, transform, field[0], field[1], field[2], (mia::Interpolation)o_interpolation, o_fill_value);

  fs::path output_path(filename_output);
  std::string basename = fs::basename(output_path);
  if (fs::extension(basename) != "") basename = fs::basename(basename);

  std::string out_folder = output_path.parent_path().string();
  if (out_folder == "") out_folder = ".";
  else if (!fs::exists(out_folder)) fs::create_directories(out_folder);

  target.dataType(source.dataType());
  itkio::save(target, filename_output);

  if (run_com_alignment || run_linear)
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

  if (run_nonlinear)
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
