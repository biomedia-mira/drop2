#include "RandomField.h"

namespace mrfopt
{
  Graph grid_graph_first_order(int nodesX, int nodesY, int nodesZ)
  {
    Graph graph;
    graph.num_nodes(nodesX * nodesY * nodesZ);

    for (int z = 0; z < nodesZ; z++)
    {
      for (int y = 0; y < nodesY; y++)
      {
        for (int x = 0; x < nodesX; x++)
        {
          int index = x + y * nodesX + z * nodesX * nodesY;
          if (x < nodesX - 1)
          {
            graph.add_clique(Clique(index, index + 1));
          }
          if (y < nodesY - 1)
          {
            graph.add_clique(Clique(index, index + nodesX));
          }
          if (z < nodesZ - 1)
          {
            graph.add_clique(Clique(index, index + nodesX * nodesY));
          }
        }
      }
    }

    return graph;
  }

  Graph grid_graph_second_order(int nodesX, int nodesY, int nodesZ)
  {
    Graph graph;
    graph.num_nodes(nodesX * nodesY * nodesZ);

    for (int z = 0; z < nodesZ; z++)
    {
      for (int y = 0; y < nodesY; y++)
      {
        for (int x = 0; x < nodesX; x++)
        {
          int index = x + y * nodesX + z * nodesX * nodesY;
          if (x < nodesX - 1)
          {
            graph.add_clique(Clique(index, index + 1));
          }
          if (y < nodesY - 1)
          {
            graph.add_clique(Clique(index, index + nodesX));
          }
          if (z < nodesZ - 1)
          {
            graph.add_clique(Clique(index, index + nodesX * nodesY));
          }
          if (x < nodesX - 2)
          {
            graph.add_clique(Clique(index, index + 1, index + 2));
          }
          if (y < nodesY - 2)
          {
            graph.add_clique(Clique(index, index + nodesX, index + nodesX * 2));
          }
          if (z < nodesZ - 2)
          {
            graph.add_clique(Clique(index, index + nodesX * nodesY, index + nodesX * nodesY * 2));
          }
        }
      }
    }

    return graph;
  }

  Graph grid_graph(int nodesX, int nodesY, int nodesZ, bool pairwise, bool triple)
  {
    Graph graph;
    graph.num_nodes(nodesX * nodesY * nodesZ);

    for (int z = 0; z < nodesZ; z++)
    {
      for (int y = 0; y < nodesY; y++)
      {
        for (int x = 0; x < nodesX; x++)
        {
          int index = x + y * nodesX + z * nodesX * nodesY;
          if (pairwise)
          {
            if (x < nodesX - 1)
            {
              graph.add_clique(Clique(index, index + 1));
            }
            if (y < nodesY - 1)
            {
              graph.add_clique(Clique(index, index + nodesX));
            }
            if (z < nodesZ - 1)
            {
              graph.add_clique(Clique(index, index + nodesX * nodesY));
            }
          }
          if (triple)
          {
            if (x < nodesX - 2)
            {
              graph.add_clique(Clique(index, index + 1, index + 2));
            }
            if (y < nodesY - 2)
            {
              graph.add_clique(Clique(index, index + nodesX, index + nodesX * 2));
            }
            if (z < nodesZ - 2)
            {
              graph.add_clique(Clique(index, index + nodesX * nodesY, index + nodesX * nodesY * 2));
            }
          }
        }
      }
    }

    return graph;
  }
}
