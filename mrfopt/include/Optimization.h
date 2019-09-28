#pragma once

#include "RandomField.h"

namespace mrfopt
{
  enum Reduction
  {
    ELC_HOCR,
    ELC_APPROX,
    HOCR,
  };

  /**
  * \brief Implements globally optimal dynamic programming for solving pairwise chain graphs.
  * \param graph The pairwise chain graph.
  * \param function The energy function.
  * \param labeling The MRF labeling.
  * \param numLabels The number of integer labels.
  **/
  template<typename energy_type>
  void multi_label_first_order_chain(const Graph& graph, const EnergyFunction<energy_type, int>& function, std::vector<int>& labeling, int numLabels);

  /**
   * \brief Implements fusion move optimization for integer-based multi-labeling problems.
   * \param graph The MRF graph.
   * \param function The energy function.
   * \param labeling The MRF labeling.
   * \param numLabels The number of integer labels.
   * \param numSweeps The number of sweeps over the label set.
   * \param numImprovements The number of QPBO-I trials.
   * \param reductionMode The higher-order clique reduction mode.
   * \param randomPermutation Enables random label set permutations for each sweep.
   **/
  template<typename energy_type>
  void multi_label(const Graph& graph, const EnergyFunction<energy_type, int>& function, std::vector<int>& labeling, int numLabels, int numSweeps, int numImprovements, Reduction reductionMode, bool randomPermutation);

  /**
   * \brief Implements a single fusion move optimization over a given labeling and a proposal labeling.
   * \param graph The MRF graph.
   * \param function The energy function.
   * \param labeling The MRF labeling.
   * \param proposal The MRF proposal labeling.
   * \param numImprovements The number of QPBO-I trials.
   * \param reductionMode The higher-order clique reduction mode.
   **/
  template<typename energy_type, typename label_type>
  void fusion_move(const Graph& graph, const EnergyFunction<energy_type, label_type>& function, std::vector<label_type>& labeling, const std::vector<label_type>& proposal, int numImprovements, Reduction reductionMode);
}

#include "Optimization.hpp"
