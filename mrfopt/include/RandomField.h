#pragma once

#include <vector>
#include <algorithm>
#include <stddef.h>

namespace mrfopt
{
  struct Clique
  {
    Clique(int nodeA, int nodeB)
    {
      nodes.resize(2);
      nodes[0] = nodeA;
      nodes[1] = nodeB;
      std::sort(nodes.begin(), nodes.end());
    }

    Clique(int nodeA, int nodeB, int nodeC)
    {
      nodes.resize(3);
      nodes[0] = nodeA;
      nodes[1] = nodeB;
      nodes[2] = nodeC;
      std::sort(nodes.begin(), nodes.end());
    }

    Clique(const std::vector<int>& cliqueNodes)
    {
      nodes = cliqueNodes;
      std::sort(nodes.begin(), nodes.end());
    }

    size_t size() const
    {
      return nodes.size();
    }

    std::vector<int> nodes;
  };

  class Graph
  {
  public:

    /**
     * \brief Adds a clique to the graph.
     * \param clique The clique to be added.
     **/
    void add_clique(const Clique& clique)
    {
      m_cliques.push_back(clique);
    }

    /**
     * \brief Sets the number of graph nodes.
     * \param nodes The number of nodes.
     **/
    void num_nodes(int numNodes)
    {
      m_num_nodes = numNodes;
    }

    /**
     * \brief Returns the number of graph nodes.
     * \return The number of nodes.
     **/
    int num_nodes() const
    {
      return m_num_nodes;
    }

    /**
     * \brief Returns the number of cliques.
     * \return The number of cliques.
     **/
    size_t num_cliques() const
    {
      return m_cliques.size();
    }

    /**
     * \brief Returns the cliques.
     * \return The vector of cliques.
     **/
    const std::vector<Clique>& cliques() const
    {
      return m_cliques;
    }

  protected:
    int m_num_nodes;
    std::vector<Clique> m_cliques;
  };

  template<typename energy_type, typename label_type>
  class EnergyFunction
  {
  public:
    /**
     * \brief Computes the unary potential for a given node and label.
     * \param nodeIndex The node index.
     * \param label The label.
     * \return The potential.
     **/
    virtual energy_type unary_potential(int nodeIndex, label_type label) const = 0;

    /**
     * \brief Computes the clique potential for a given clique and labeling.
     * \param clique The clique.
     * \param labels The clique labeling.
     * \return The potential.
     **/
    virtual energy_type clique_potential(const Clique& clique, const std::vector<label_type>& labels) const = 0;
  };

  /**
   * \brief Computes the energy of a graph labeling.
   * \param graph The MRF graph.
   * \param function The energy function.
   * \param labeling The MRF labeling.
   * \return The energy.
   **/
  template <typename energy_type, typename label_type>
  energy_type compute_energy(const Graph& graph, const EnergyFunction<energy_type, label_type>& function, const std::vector<label_type>& labeling);

  /**
   * \brief Generates a first-order grid graph with pairwise cliques.
   * \param nodesX The number of nodes on the x-axis.
   * \param nodesY The number of nodes on the y-axis.
   * \param nodesZ The number of nodes on the z-axis.
   * \return The MRF graph.
   **/
  Graph grid_graph_first_order(int nodesX, int nodesY, int nodesZ);

  /**
   * \brief Generates a second-order grid graph with pairwise and triple cliques.
   * \param nodesX The number of nodes on the x-axis.
   * \param nodesY The number of nodes on the y-axis.
   * \param nodesZ The number of nodes on the z-axis.
   * \return The MRF graph.
   **/
  Graph grid_graph_second_order(int nodesX, int nodesY, int nodesZ);

  /**
  * \brief Generates a grid graph with selected cliques.
  * \param nodesX The number of nodes on the x-axis.
  * \param nodesY The number of nodes on the y-axis.
  * \param nodesZ The number of nodes on the z-axis.
  * \param pairwise Enables pairwise cliques.
  * \param triple Enables triple cliques.
  * \return The MRF graph.
  **/
  Graph grid_graph(int nodesX, int nodesY, int nodesZ, bool pairwise, bool triple);
}

#include "RandomField.hpp"
