// This an MRF optimization example for image denoising.
//
// The example illustrates the use of fusion moves for solving integer multi-labeling problems.
// The first-order energy model is taken from: http://vision.middlebury.edu/MRF/results/denoise/penguin/
// Additionally, a second-order term is implemented which penalizes second derivatives
//
// A test image is provided in the data folder.
// The suggested values for the pairwise potential are: truncation = 200, lambda = 25
// The suggested values for the triplet potential are: truncation = 500, lambda = 10
//

#include "RandomField.h"
#include "Optimization.h"
#include "miaImage.h"
#include "itkio.h"

#include <vector>
#include <iostream>
#include <chrono>

using namespace mrfopt;
using namespace mia;

class DenoisingEnergyFunction : public EnergyFunction<int, int>
{
  public:
    DenoisingEnergyFunction(Image image, Image mask, int truncationPairwise, int lambdaPairwise, int truncationTriplet, int lambdaTriplet)
      : m_image(image)
      , m_mask(mask)
      , m_truncation_pairwise(truncationPairwise)
      , m_lambda_pairwise(lambdaPairwise)
      , m_truncation_triplet(truncationTriplet)
      , m_lambda_triplet(lambdaTriplet)
    {}

    int unary_potential(int nodeIndex, int label) const override
    {
      if (m_mask.data()[nodeIndex] > 0)
        return static_cast<int>(pow(static_cast<double>(m_image.data()[nodeIndex]) - static_cast<double>(label), 2.0));
      else
        return 0;
    }

    int clique_potential(const Clique& clique, const std::vector<int>& labels) const override
    {
      auto clique_size = static_cast<int>(clique.size());
      if (clique_size == 2)
      {
        return pairwise_potential(clique.nodes[0], clique.nodes[1], labels[0], labels[1]);
      }
      else if (clique_size == 3)
      {
        return triplet_potential(clique.nodes[0], clique.nodes[1], clique.nodes[2], labels[0], labels[1], labels[2]);
      }
      else
      {
        return 0;
      }
    }

    int pairwise_potential(int nodeA, int nodeB, int labelA, int labelB) const
    {
      if (m_lambda_pairwise == 0) return 0;

      auto potential = static_cast<int>(pow(static_cast<double>(labelA - labelB), 2.0));
      auto energy = std::min(potential, m_truncation_pairwise);
      return energy * m_lambda_pairwise;
    }

    int triplet_potential(int nodeA, int nodeB, int nodeC, int labelA, int labelB, int labelC) const
    {
      if (m_lambda_triplet == 0) return 0;

      auto potential = static_cast<int>(pow(static_cast<double>(labelA - 2*labelB + labelC), 2.0));
      auto energy = std::min(potential, m_truncation_triplet);
      return energy * m_lambda_triplet;
    }

  private:
    Image m_image;
    Image m_mask;
    int m_truncation_pairwise;
    int m_lambda_pairwise;
    int m_truncation_triplet;
    int m_lambda_triplet;
};

int main(int argc, char* argv[])
{
  std::string file_image = "../../../../../source/mira/mrfopt/examples/denoising/data/penguin.png";
  std::string file_mask = "../../../../../source/mira/mrfopt/examples/denoising/data/mask.png";
  int truncation_pairwise = 200;
  int lambda_pairwise = 25;
  int truncation_triplet = 500;
  int lambda_triplet = 0; //10;

  int numLabels = 256;
  int numSweeps = 2;
  int numImprovements = 0;
  Reduction reductionMode = HOCR;
  bool randomPermutation = false;

  Image image = itkio::load(file_image);
  Image mask = itkio::load(file_mask);

  Graph graph = grid_graph_second_order(image.sizeX(), image.sizeY(), image.sizeZ());

  DenoisingEnergyFunction function(image, mask, truncation_pairwise, lambda_pairwise, truncation_triplet, lambda_triplet);

  std::vector<int> labeling(graph.num_nodes());
  std::fill(labeling.begin(), labeling.end(), 0);

  auto start = std::chrono::high_resolution_clock::now();

  auto energy_initial = compute_energy(graph, function, labeling);
  std::cout << "Initial energy:\t" << energy_initial << std::endl;

  multi_label(graph, function, labeling, numLabels, numSweeps, numImprovements, reductionMode, randomPermutation);

  auto energy_final = compute_energy(graph, function, labeling);
  std::cout << "Final energy:\t" << energy_final << std::endl;

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Timing (ms): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop-start).count() << std::endl;


  Image result = image.clone();
  for (int i = 0; i < image.size(); i++)
  {
    result.data()[i] = labeling[i];
  }

  itkio::save(result, "result.png");
}
