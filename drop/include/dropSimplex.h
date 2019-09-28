#pragma once

#include <vector>
#include <memory>

namespace drop
{
  class CostFunctionSimplex
  {
  public:
    virtual double evaluate(const std::vector<double> &parameters) = 0;
  };

  class SimplexOptimizer
  {
  public:

    enum Termination
    {
      OPT_SMALL_ERROR,
      OPT_SMALL_PARAMETER_DELTA,
      OPT_MAX_EVALS_REACHED,
      OPT_MAX_ITERS_REACHED,
      OPT_UNKNOWN_REASON,
    };

    SimplexOptimizer(int numParameters, const std::vector<double>& step_sizes);

    /**
    * \brief Sets the cost function.
    * \param function The cost function.
    **/
    void cost_function(std::shared_ptr<CostFunctionSimplex> function)
    {
      costfct = function;
    }

    /**
    * \brief Gets the parameter vector.
    * \return A vector of parameters.
    **/
    const std::vector<double>& parameters() const
    {
      return m_parameters;
    }

    /**
    * \brief Sets the maximum number of iterations.
    * \param value The maximum number of iterations.
    **/
    void max_iterations(int value)
    {
      m_max_iters = value;
    }

    /**
    * \brief Sets the maximum number of evaluations.
    * \param value The maximum number of evaluations.
    **/
    void max_evaluations(int value)
    {
      m_max_evals = value;
    }

    /**
    * \brief Sets the function tolerance stopping criterion.
    * \param value The function tolerance.
    **/
    void function_tolerance(double value)
    {
      m_function_tolerance = value;
    }

    /**
    * \brief Sets the parameter tolerance stopping criterion.
    * \param value The parameter tolerance.
    **/
    void parameter_tolerance(double value)
    {
      m_parameter_tolerance = value;
    }

    /**
    * \brief Gets the number of iterations.
    * \return The number of iterations.
    **/
    int num_iterations() const
    {
      return m_num_iters;
    }

    /**
    * \brief Gets the number of evaluations.
    * \return The number of evaluations.
    **/
    int num_evaluations() const
    {
      return m_num_evals;
    }

    /**
    * \brief Runs the optimization.
    **/
    Termination run();

  private:

    double amotry(int jhi, double factor);

    void sum();

    std::shared_ptr<CostFunctionSimplex> costfct;
    std::vector<double> m_parameters;
    std::vector<double> m_step_sizes;
    std::vector<double> m_simplex_try;
    std::vector<double> m_simplex_sum;
    std::vector<double> m_simplex_rst;
    std::vector<std::vector<double>> m_simplex;

    Termination m_termination;
    int m_num_parameters;
    int m_num_points;
    int m_max_evals;
    int m_max_iters;
    int m_num_iters;
    int m_num_evals;
    double m_ftol;
    double m_ptol;
    double m_function_tolerance;
    double m_parameter_tolerance;
    bool m_abort;
  };
}
