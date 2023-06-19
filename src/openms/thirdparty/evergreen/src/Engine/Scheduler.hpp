#ifndef _SCHEDULER_HPP
#define _SCHEDULER_HPP

#include "Edge.hpp"
#include "MessagePasser.hpp"
#include "InferenceGraph.hpp"

template <typename VARIABLE_KEY>
class Scheduler {
protected:
  // dampening_lambda = 0.0 uses only _current_message; with 1.0, it
  // only uses _old_message.
  double _dampening_lambda;
  double _convergence_threshold;
  unsigned long _maximum_iterations;

public:
  Scheduler(double dampening_lambda_param, double convergence_threshold_param, unsigned long maximum_iterations_param):
    _dampening_lambda(dampening_lambda_param),
    _convergence_threshold(convergence_threshold_param),
    _maximum_iterations(maximum_iterations_param)
  {
    assert(_dampening_lambda < 0.5 && "Dampening should be performed with lambda < 0.5 (higher lambda values will weight older messages over new messages, and may lead to oscillations [unproven])");
  }

  virtual ~Scheduler() {}

  double dampening_lambda() const {
    return _dampening_lambda;
  }

  double convergence_threshold() const {
    return _convergence_threshold;
  }

  void set_dampening_lambda(double lambda) {
    _dampening_lambda = lambda;
  }

  void set_convergence_threshold(double epsilon) {
    _convergence_threshold = epsilon;
  }

  void set_maximum_iterations(unsigned long n) {
    _maximum_iterations = n;
  }

  // Returns the number of iterations spent:
  virtual unsigned long process_next_edges() = 0;

  virtual bool has_converged() const = 0;
  // Add messages to the queue that are elligible to pass from the start:
  virtual void add_ab_initio_edges(InferenceGraph<VARIABLE_KEY> &) = 0;

  // Returns the number of iterations spent:
  virtual unsigned long run_until_convergence() {
    unsigned long iteration;
    for (iteration = 0; ! has_converged() && iteration < this->_maximum_iterations; ) {
      iteration += process_next_edges();
    }

    if (iteration >= this->_maximum_iterations)
      std::cerr << "Warning: Did not meet desired convergence threshold (stopping anyway after exceeding " << this->_maximum_iterations << " iterations)." << std::endl;
    return iteration;
  }
};

#endif
