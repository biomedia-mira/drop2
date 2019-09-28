#include "dropSimplex.h"
#include <cmath>

namespace drop
{
  SimplexOptimizer::SimplexOptimizer(int numParameters, const std::vector<double>& step_sizes)
  {
    m_parameters = std::vector<double>(numParameters);
    for (auto &p : m_parameters) p = 0.0;

    m_step_sizes = step_sizes;

    m_num_parameters = numParameters;
    m_num_points = m_num_parameters + 1;

    m_simplex_try = std::vector<double>(m_num_parameters);
    m_simplex_sum = std::vector<double>(m_num_parameters);
    m_simplex_rst = std::vector<double>(m_num_points);
    m_simplex = std::vector<std::vector<double>>(m_num_points);
    for (int j = 0; j < m_num_points; j++)
    {
      m_simplex[j] = std::vector<double>(m_num_parameters);
    }

    m_max_evals = 50000;
    m_max_iters = 10000;
    m_function_tolerance = 1e-20;
    m_parameter_tolerance = 1e-15;
    m_abort = false;
    m_termination = Termination::OPT_UNKNOWN_REASON;
  }

  SimplexOptimizer::Termination SimplexOptimizer::run()
  {
    m_num_evals = 0;
    m_num_iters = 0;

    for (int i = 0; i < m_num_parameters; i++)
    {
      m_simplex[0][i] = m_parameters[i];
    }

    m_simplex_rst[0] = costfct->evaluate(m_simplex[0]);

    for (int j = 1; j < m_num_points; j++)
    {
      for (int i = 0; i < m_num_parameters; i++)
      {
        m_simplex[j][i] = ((j - 1) == i) ? (m_simplex[0][i] + m_step_sizes[i]) : (m_simplex[0][i]);
      }
      m_simplex_rst[j] = costfct->evaluate(m_simplex[j]);
    }

    sum();

    int jhi, jnhi, jlo = 0;
    double simplexRstSave, simplexRstTry;

    while (!m_abort)
    {
      jlo = 1;

      //determine highest(worst = jhi), next-highest (jnhi) and lowest (best = jlo) point:
      if (m_simplex_rst[0] > m_simplex_rst[1])
      {
        jhi = 1;
        jnhi = 2;
      }
      else
      {
        jhi = 2;
        jnhi = 1;
      }

      for (int j = 0; j < m_num_points; j++)
      {
        if (m_simplex_rst[j] <= m_simplex_rst[jlo])
        {
          jlo = j;
        }
        if (m_simplex_rst[j] > m_simplex_rst[jhi])
        {
          jnhi = jhi;
          jhi = j;
        }
        else if (m_simplex_rst[j] > m_simplex_rst[jnhi] && j != jhi)
        {
          jnhi = j;
        }
      }

      // this is better than using TINY, since depending on the function TINY
      // might end up being actually rather big...
      m_ftol = 0;
      double absSum = std::abs(m_simplex_rst[jhi]) + std::abs(m_simplex_rst[jlo]);
      if (absSum != 0)
      {
        m_ftol = 2 * std::abs(m_simplex_rst[jhi] - m_simplex_rst[jlo]) / absSum;
      }

      if (m_ftol < m_function_tolerance)
      {
        m_simplex_rst[0] = m_simplex_rst[jlo];

        for (int i = 0; i < m_num_parameters; i++)
        {
          m_simplex[0][i] = m_simplex[jlo][i];
        }

        m_termination = Termination::OPT_SMALL_ERROR;

        break;
      }

      m_ptol = 0;
      for (int j = 0; j < m_num_points - 1; j++)
      {
        double diff;
        double norm = 0;
        for (int i = 0; i < m_num_parameters; i++)
        {
          diff = m_simplex[j][i] - m_simplex[j + 1][i];
          norm += diff * diff;
        }
        m_ptol = std::fmax(std::sqrt(norm), m_ptol);
      }

      if (m_ptol < m_parameter_tolerance)
      {
        m_simplex_rst[0] = m_simplex_rst[jlo];

        for (int i = 0; i < m_num_parameters; i++)
        {
          m_simplex[0][i] = m_simplex[jlo][i];
        }

        m_termination = Termination::OPT_SMALL_PARAMETER_DELTA;

        break;
      }

      if (m_num_evals >= m_max_evals)
      {
        m_termination = Termination::OPT_MAX_EVALS_REACHED;

        break;
      }

      m_num_evals += 2;

      if (m_num_iters >= m_max_iters)
      {
        m_termination = Termination::OPT_MAX_ITERS_REACHED;

        break;
      }

      //new iteration
      simplexRstTry = amotry(jhi, -1);
      if (simplexRstTry <= m_simplex_rst[jlo])
      {
        // better result than best point --> extrapolation by factor 2
        simplexRstTry = amotry(jhi, 2);
      }
      else if (simplexRstTry >= m_simplex_rst[jnhi])
      {
        // worse than second highest		--> 1-dim contraction
        simplexRstSave = m_simplex_rst[jhi];
        simplexRstTry = amotry(jhi, 0.5);

        if (simplexRstTry >= simplexRstSave)
        {
          //cannot get rid of point		--> contract around best point
          for (int j = 0; j < m_num_points; j++)
          {
            if (j != jlo)
            {
              for (int i = 0; i < m_num_parameters; i++)
              {
                m_simplex[j][i] = m_simplex_sum[i] = 0.5 * (m_simplex[j][i] + m_simplex[jlo][i]);
              }
              m_simplex_rst[j] = costfct->evaluate(m_simplex_sum);
            }
          }
          m_num_evals += m_num_parameters;
          sum();
        }
      }
      else
      {
        m_num_evals--;
      }

      m_num_iters++;
    }

    for (int i = 0; i < m_num_parameters; i++)
    {
      m_parameters[i] = m_simplex[jlo][i];
    }

    return m_termination;
  }

  double SimplexOptimizer::amotry(int jhi, double factor)
  {
    double fac1 = (1.0 - factor) / static_cast<double>(m_num_parameters);
    double fac2 = fac1 - factor;
    for (int i = 0; i < m_num_parameters; i++)
    {
      m_simplex_try[i] = m_simplex_sum[i] * fac1 - m_simplex[jhi][i] * fac2;
    }

    double simplexResultTry = costfct->evaluate(m_simplex_try);

    if (simplexResultTry < m_simplex_rst[jhi])
    {
      m_simplex_rst[jhi] = simplexResultTry;
      for (int i = 0; i < m_num_parameters; i++)
      {
        m_simplex_sum[i] += m_simplex_try[i] - m_simplex[jhi][i];
        m_simplex[jhi][i] = m_simplex_try[i];
      }
    }

    return simplexResultTry;
  }

  void SimplexOptimizer::sum()
  {
    for (int i = 0; i < m_num_parameters; i++)
    {
      m_simplex_sum[i] = 0;
      for (int j = 0; j < m_num_points; j++)
      {
        m_simplex_sum[i] += m_simplex[j][i];
      }
    }
  }
}
