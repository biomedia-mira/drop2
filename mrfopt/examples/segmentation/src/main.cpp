// This an MRF optimization example for image segmentation.
//
// The example illustrates the use of a single fusion move for solving binary labeling problems.
// The energy model is based on pre-computed unary likelihoods and contrast-sensitive pairwise penalties.
//
// A test image is provided in the data folder.
// The suggested values for the pairwise potential are: sigma = 50, lambda = 10
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

class SegmentationEnergyFunction : public EnergyFunction<double, int>
{
  public:
    SegmentationEnergyFunction(Image image, Image probabilities, double sigma, double lambda)
      : m_image(image)
      , m_probs(probabilities)
      , m_sigma(sigma)
      , m_lambda(lambda)
    {}

    double unary_potential(int nodeIndex, int label) const override
    {
      auto likelihood = m_probs.data()[nodeIndex];
      if (likelihood > 0.8 || likelihood < 0.2)
        return (label == 0) ? -log(1.0-likelihood + std::numeric_limits<double>::epsilon()) : -log(likelihood + std::numeric_limits<double>::epsilon());
      else
        return 0;
    }

    double clique_potential(const Clique& clique, const std::vector<int>& labels) const override
    {
      auto clique_size = static_cast<int>(clique.size());
      if (clique_size == 2)
      {
        return pairwise_potential(clique.nodes[0], clique.nodes[1], labels[0], labels[1]);
      }
      else
      {
        return 0;
      }
    }

    double pairwise_potential(int nodeA, int nodeB, int labelA, int labelB) const
    {
      if (labelA == labelB) return 0;

      double dist = (nodeB - nodeA > 1) ? m_image.spacing()[1] : m_image.spacing()[0];

      double valueA = m_image.data()[nodeA];
      double valueB = m_image.data()[nodeB];

      const double diffSq = (valueA - valueB) * (valueA - valueB);

      return m_lambda * exp(- diffSq / (2 * m_sigma * m_sigma)) / dist;
    }

  private:
    Image m_image;
    Image m_probs;
    double m_sigma;
    double m_lambda;
};

int main(int argc, char* argv[])
{
  std::string file_image = "../../../../../source/mira/mrfopt/examples/segmentation/data/image.nii.gz";
  std::string file_probmap = "../../../../../source/mira/mrfopt/examples/segmentation/data/likelihood.nii.gz";
  double sigma = 50;
  double lambda = 10;

  int numImprovements = 0;
  Reduction reductionMode = HOCR;

  Image image = itkio::load(file_image);
  Image probs = itkio::load(file_probmap);

  Graph graph = grid_graph_first_order(image.sizeX(), image.sizeY(), image.sizeZ());

  SegmentationEnergyFunction function(image, probs, sigma, lambda);

  std::vector<int> labeling(graph.num_nodes());
  std::fill(labeling.begin(), labeling.end(), 0);

  std::vector<int> proposal(graph.num_nodes());
  std::fill(proposal.begin(), proposal.end(), 1);

  auto start = std::chrono::high_resolution_clock::now();

  auto energy_initial = compute_energy(graph, function, labeling);
  std::cout << "Initial energy:\t" << energy_initial << std::endl;

  fusion_move(graph, function, labeling, proposal, numImprovements, reductionMode);

  auto energy_final = compute_energy(graph, function, labeling);
  std::cout << "Final energy:\t" << energy_final << std::endl;

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Timing (ms): " << std::chrono::duration_cast< std::chrono::milliseconds >(stop-start).count() << std::endl;

  Image result = image.clone();
  for (int i = 0; i < image.size(); i++)
  {
    result.data()[i] = labeling[i];
  }

  itkio::save(result, "result.nii.gz");
}
