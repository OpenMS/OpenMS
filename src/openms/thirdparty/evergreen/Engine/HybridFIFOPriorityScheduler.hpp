#ifndef _HYBRIDFIFOPRIORITYSCHEDULER_HPP
#define _HYBRIDFIFOPRIORITYSCHEDULER_HPP

#include "FIFOScheduler.hpp"
#include "PriorityScheduler.hpp"

// Combines FIFOScheduler and PriorityScheduler: first run
// FIFOScheduler (this is useful for solving ).
template <typename VARIABLE_KEY>
class HybridFIFOPriorityScheduler : public Scheduler<VARIABLE_KEY> {
protected:
  const InferenceGraph<VARIABLE_KEY> & _graph;
  FIFOScheduler<VARIABLE_KEY>* fs;
  PriorityScheduler<VARIABLE_KEY>* ps;

public:
  HybridFIFOPriorityScheduler(double dampening_lambda, double convergence_threshold, unsigned long maximum_iterations, const InferenceGraph<VARIABLE_KEY> & graph):
    Scheduler<VARIABLE_KEY>(dampening_lambda, convergence_threshold, maximum_iterations),
    _graph(graph),
    fs(NULL),
    ps(NULL)
  { }

  ~HybridFIFOPriorityScheduler() {
    if (fs != NULL)
      delete fs;
    if (ps != NULL)
      delete ps;
  }
  
  bool process_next_edges() {
    if (!fs.has_converged())
      fs->process_next_edges();
    else
      ps->process_next_edges();
  }

  unsigned long run_until_convergence() {
    // Use inf as the convergence threshold for FIFOScheduler,
    // guaranteening that every reachable edge is only visited once:
    fs = new FIFOScheduler<VARIABLE_KEY> (this->_dampening_lambda, std::numeric_limits<double>::infinity(), this->_maximum_iterations, _graph);
    // First, send messages across every reachable edge via the
    // FIFOScheduler:
    unsigned long iterations_used = fs->run_until_convergence();

    // iterations_used should be <= this->_maximum_iterations, so subtracting is safe:
    PriorityScheduler<VARIABLE_KEY> ps(this->_dampening_lambda, this->_convergence_threshold, this->_maximum_iterations - iterations_used, _graph);
    return iterations_used + ps.run_until_convergence();
  }

  bool has_converged() const {
    return ps->has_converged();
  }
};

#endif
