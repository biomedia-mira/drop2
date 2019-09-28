namespace mrfopt
{
  template <typename energy_type, typename label_type>
  energy_type compute_energy(const Graph& graph, const EnergyFunction<energy_type, label_type>& function, const std::vector<label_type>& labeling)
  {
    energy_type energy = 0;

    auto num_nodes = graph.num_nodes();
    for (int node = 0; node < num_nodes; node++)
    {
      energy += function.unary_potential(node, labeling[node]);
    }

    auto cliques = graph.cliques();
    for( auto& clique : cliques)
    {
      auto clique_size = static_cast<int>(clique.size());
      std::vector<label_type> clique_labels(clique_size);
      for (int i = 0; i < clique_size; i++)
      {
        clique_labels[i] = labeling[clique.nodes[i]];
      }

      energy += function.clique_potential(clique, clique_labels);
    }

    return energy;
  }
}
