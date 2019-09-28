#include "QPBO.h"
#include "ELC.h"
#include <time.h>
#include <random>
#include <vector>

using namespace ELCReduce;

namespace mrfopt
{
  namespace impl
  {
    void random_permutation(std::vector<int>& labels);

    template<typename energy_type>
    void reduce_and_convert(PBF<energy_type>& pbf, QPBO<energy_type>& qpbo, Reduction mode)
    {
      switch(mode)
      {
      case ELC_HOCR:
      {
        PBF<energy_type> qpbf;
        pbf.reduceHigher(); // Use the ELC technique to reduce higher-order terms without auxiliary variables
        pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce the remaining higher-order terms using HOCR adding auxiliary variables
        qpbf.convert(qpbo, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;
      }
      case ELC_APPROX:
      {
        pbf.reduceHigherApprox(); // Use the approximate ELC technique to reduce higher-order terms without auxiliary variables
        pbf.convert(qpbo, pbf.maxID()+1);
        pbf.clear();
        break;
      }
      case HOCR:
      {
        PBF<energy_type> qpbf;
        pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce to Quadratic pseudo-Boolean function using HOCR.
        qpbf.convert(qpbo, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;
      }
      }
    }

    template<typename energy_type, typename label_type>
    void fusion_move(QPBO<energy_type>& qpbo, const Graph& graph, const EnergyFunction<energy_type, label_type>& function, std::vector<label_type>& labeling, const std::vector<label_type>& proposal, int numImprovements, Reduction reductionMode)
    {
      PBF<energy_type> pbf;

      auto num_nodes = graph.num_nodes();
      for (int node = 0; node < num_nodes; node++)
      {
        energy_type e0 = function.unary_potential(node, labeling[node]);
        energy_type e1 = function.unary_potential(node, proposal[node]);
        pbf.AddUnaryTerm(node, e0, e1);
      }

      auto cliques = graph.cliques();
      for( auto& clique : cliques)
      {
        auto clique_size = static_cast<int>(clique.size());
        std::vector<label_type> clique_labeling(clique_size);
        std::vector<label_type> clique_proposal(clique_size);
        for (int i = 0; i < clique_size; i++)
        {
          clique_labeling[i] = labeling[clique.nodes[i]];
          clique_proposal[i] = proposal[clique.nodes[i]];
        }

        auto combinations = static_cast<int>(pow(2.0, static_cast<double>(clique_size)));
        std::vector<energy_type> energies(combinations);
        std::vector<label_type> clique_labels(clique_size);
        for (int c = 0; c < combinations; c++)
        {
          for (int i = 0; i < clique_size; i++)
          {
            bool b = (c >> (clique_size - i - 1)) % 2 == 1;
            clique_labels[i] = b ? clique_proposal[i] : clique_labeling[i];
          }
          energies[c] = function.clique_potential(clique, clique_labels);
        }

        if (clique_size == 2)
        {
          pbf.AddPairwiseTerm(clique.nodes[0], clique.nodes[1], energies[0], energies[1], energies[2], energies[3]);
        }
        else if (clique_size > 2)
        {
          pbf.AddHigherTerm(static_cast<int>(clique_size), &clique.nodes[0], &energies[0]);
        }
      }

      qpbo.Reset();

      reduce_and_convert(pbf, qpbo, reductionMode);

      // MergeParallelEdges causes problems when using ELC_HOCR and ELC_APPROX reduction modes in some cases.
      // Should not be needed if AddPairwiseTerm(i,j) is not called twice or more for the nodes i and j.
      //qpbo.MergeParallelEdges();
      qpbo.Solve();
      qpbo.ComputeWeakPersistencies();

      srand ( static_cast<unsigned int>(time(NULL)) );
      for (int i = 0; i < numImprovements; i++)
      {
        qpbo.Improve();
      }

      for (int node = 0; node < num_nodes; node++)
      {
        if (qpbo.GetLabel(node) == 1)
          labeling[node] = proposal[node];
      }
    }
  }

  template<typename energy_type>
  void multi_label_first_order_chain(const Graph& graph, const EnergyFunction<energy_type, int>& function, std::vector<int>& labeling, int numLabels)
  {
    int num_nodes = graph.num_nodes();

    std::vector<std::vector<int>> labeling_table(num_nodes);
    std::vector<std::vector<energy_type>> cumsum_table(num_nodes);

    std::vector<int> clique_labels(2);

    for (int i = 0; i < num_nodes; i++)
    {
      labeling_table[i] = std::vector<int>(numLabels);
      cumsum_table[i] = std::vector<energy_type>(numLabels);
    }

    for (int label = 0; label < numLabels; label++)
    {
      cumsum_table[0][label] = function.unary_potential(0, label);
    }

    for (int index = 1; index < num_nodes; index++)
    {
      for (int label = 0; label < numLabels; label++)
      {
        int bestPred = 0;
        clique_labels[1] = label;
        energy_type sumBestPred = std::numeric_limits<energy_type>::max();
        for (int labelPred = 0; labelPred < numLabels; labelPred++)
        {
          clique_labels[0] = labelPred;
          energy_type newCumSum = cumsum_table[index - 1][labelPred] + function.clique_potential(Clique(index - 1, index), clique_labels);

          if (newCumSum < sumBestPred)
          {
            bestPred = labelPred;
            sumBestPred = newCumSum;
          }
        }

        labeling_table[index][label] = bestPred;
        cumsum_table[index][label] = sumBestPred + function.unary_potential(index, label);
      }
    }

    int bestLabelLast = 0;
    energy_type minSum = cumsum_table[num_nodes - 1][bestLabelLast];
    for (int label = 1; label < numLabels; label++)
    {
      energy_type sum = cumsum_table[num_nodes - 1][label];
      if (sum < minSum)
      {
        minSum = sum;
        bestLabelLast = label;
      }
    }

    labeling[num_nodes - 1] = bestLabelLast;

    for (int index = num_nodes - 2; index >= 0; index--)
    {
      labeling[index] = labeling_table[index + 1][labeling[index + 1]];
    }
  }

  template<typename energy_type>
  void multi_label(const Graph& graph, const EnergyFunction<energy_type, int>& function, std::vector<int>& labeling, int numLabels, int numSweeps, int numImprovements, Reduction reductionMode, bool randomPermutation)
  {
    QPBO<energy_type> qpbo(graph.num_nodes(), graph.num_nodes() * 4);

    std::vector<int> labels(numLabels);
    for (int i = 0; i < numLabels; i++)
    {
      labels[i] = i;
    }

    std::vector<int> proposal(graph.num_nodes());
    for (int i = 0; i < numSweeps; i++)
    {
      if (randomPermutation) impl::random_permutation(labels);
      for (int labelIndex = 0; labelIndex < numLabels; labelIndex++)
      {
        std::fill(proposal.begin(), proposal.end(), labels[labelIndex]);
        impl::fusion_move(qpbo, graph, function, labeling, proposal, numImprovements, reductionMode);
      }
    }
  }

  template<typename energy_type, typename label_type>
  void fusion_move(const Graph& graph, const EnergyFunction<energy_type, label_type>& function, std::vector<label_type>& labeling, const std::vector<label_type>& proposal, int numImprovements, Reduction reductionMode)
  {
    QPBO<energy_type> qpbo(graph.num_nodes(), graph.num_nodes() * 4);
    impl::fusion_move(qpbo, graph, function, labeling, proposal, numImprovements, reductionMode);
  }
}
