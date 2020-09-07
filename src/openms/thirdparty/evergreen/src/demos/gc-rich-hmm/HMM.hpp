#ifndef _ISOTOPESOLVER_HPP
#define _ISOTOPESOLVER_HPP

#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/Clock.hpp"

class HMM {
private:
  PMF _prior;
  PMF _transition;
  PMF _emission;
  //  std::vector<unsigned long> _hidden_variables;
  //  std::vector<unsigned long> _observed_variables;
  std::vector<std::string> _hidden_variables;
  std::vector<std::string> _observed_variables;

  const std::string & _evidence;
  InferenceGraph<std::string> *_ig;

  Scheduler<std::string> & _sched;

  // Note: could be done faster via a table of 256 chars --> indices
  // in {G,A,T,C}.
  PMF create_nucleotide_evidence_pmf(char gatc) {
    Vector<double> evidence({0.0, 0.0, 0.0, 0.0});
    
    switch (gatc) {
    case 'G':
      evidence[0] = 1.0;
      break;
    case 'A':
      evidence[1] = 1.0;
      break;
    case 'T':
      evidence[2] = 1.0;
      break;
    case 'C':
      evidence[3] = 1.0;
      break;
    default:
      assert(false && "Not a valid nucleotide 'G' 'A' 'T' or 'C'");
      break;
    }

    return PMF({0L}, Tensor<double>({4ul}, evidence));    
  }

  void construct_graph(double p) {

    //    std::vector<MessagePasser<unsigned long>* > mps;
    std::vector<MessagePasser<std::string>* > mps;
    
    const unsigned long n = _evidence.size();

    HUGINMessagePasser<std::string>*current_node = new HUGINMessagePasser<std::string>(LabeledPMF<std::string>({_hidden_variables[0]}, _prior), p);
    
    for(unsigned long i=0; i<n; ++i) {
      // Create observed DNA evidence:
      HUGINMessagePasser<std::string>*hmp_data = new HUGINMessagePasser<std::string>(LabeledPMF<std::string>({_observed_variables[i]}, create_nucleotide_evidence_pmf(_evidence[i])), p);
      mps.push_back(hmp_data);

      // Create emission between hypotheses and observed DNA evidence:
      HUGINMessagePasser<std::string>*hmp_emission = new HUGINMessagePasser<std::string>(LabeledPMF<std::string>({_hidden_variables[i], _observed_variables[i]}, _emission), p);

      mps.push_back(hmp_emission);
      hmp_emission->bind_to(hmp_data, new std::vector<std::string>{_observed_variables[i]});
      
      // Note: the above two HUGINMessagePasser types could be
      // compressed into one, which basically inlines hmp_emission
      // conditional on data=_evidence[i].

      current_node->bind_to(hmp_emission, new std::vector<std::string>{_hidden_variables[i]});
      mps.push_back(current_node);

      // Create transition to next nucleotide (if not at the final node):
      if (i+1 < n) {
	HUGINMessagePasser<std::string>*hmp_transition = new HUGINMessagePasser<std::string>(LabeledPMF<std::string>({_hidden_variables[i], _hidden_variables[i+1]}, _transition), p);
	current_node->bind_to(hmp_transition, new std::vector<std::string>{_hidden_variables[i]});
	mps.push_back(hmp_transition);

	current_node = new HUGINMessagePasser<std::string>(p);
	hmp_transition->bind_to(current_node, new std::vector<std::string>{_hidden_variables[i+1]});
      }
    }

    _ig = new InferenceGraph<std::string>(std::move(mps));
  }
  
public:

  HMM(const PMF & prior, const PMF & transition, const PMF & emission, const std::string & evidence, double p, Scheduler<std::string> & sched):
    _prior(prior),
    _transition(transition),
    _emission(emission),
    _evidence(evidence),
    _sched(sched)
  {
    for(unsigned long i=0; i<_evidence.size(); ++i) {
      _hidden_variables.push_back("H" + to_string(i));
      _observed_variables.push_back("D" + to_string(i));
    }
    
    // create inference graph
    construct_graph(p);
  }

  ~HMM() {
    delete _ig;
  }

  std::vector<LabeledPMF<std::string> > solve() {
    std::cout << "solving..." << std::endl;

    // apply belief propagation to inference graph
    _sched.add_ab_initio_edges(*_ig);
    BeliefPropagationInferenceEngine<std::string> bpie(_sched, *_ig);

    Clock c;
    std::vector<std::vector<std::string> > hidden_variable_singletons;
    for(unsigned long i=0; i<_evidence.size(); ++i)
      hidden_variable_singletons.push_back({_hidden_variables[i]});
    
    auto result = bpie.estimate_posteriors(hidden_variable_singletons);
    c.ptock();

    return result;
  }

};

#endif
