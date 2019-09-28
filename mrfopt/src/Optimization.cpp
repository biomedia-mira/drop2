#include "Optimization.h"

namespace mrfopt
{
  namespace impl
  {
    void random_permutation(std::vector<int>& labels)
    {
      std::random_device rd;
      std::mt19937 mt(rd());

      int n = static_cast<int>(labels.size());
      for (int i = 0; i < n; i++)
      {
        std::uniform_int_distribution<int> uid(i, n - 1);
        int r = uid(mt);
        std::swap(labels[i], labels[r]);
      }
    }
  }
}
