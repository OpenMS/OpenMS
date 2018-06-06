#ifndef _MASSPEPTIDESOLVER_HPP
#define _MASSPEPTIDESOLVER_HPP

#include "../../Evergreen/evergreen.hpp"

#include "../../Utility/inference_utilities.hpp"
#include "../../Utility/Clock.hpp"

#include "Peptide.hpp"

#include "../../Utility/graph_to_dot.hpp"

class MassPeptideSolver {
private:
  Scheduler<std::string> & _sched;
  InferenceGraph<std::string> *_ig_ptr;

  static constexpr double DITHERING_SIGMA = 0.1;
  // The value beyond which Gaussian tails are no longer considered:
  static constexpr double GAUSSIAN_TAIL_EPSILON = 0.005;

public:
  
  MassPeptideSolver(double mass_goal, const double & p, const unsigned int max_num_copies, const double mass_discretization, Scheduler<std::string> & sched):
    _sched(sched)
  {
    ///////////////////////////
    ///// Construct Graph /////
    ///////////////////////////
    
    BetheInferenceGraphBuilder<std::string> igb;
    
    std::vector<std::string> amino_acid_strings(Peptide::amino_acids.size());
    for (unsigned int i=0; i<Peptide::amino_acids.size(); ++i)
      amino_acid_strings[i] += Peptide::amino_acids[i];
    
    // Vectors used later on for graph construction.
    std::vector<std::vector<std::string> > aa_mass_singletons;

    //// Add Table Dependencies ////
    // Make uniform distribution for each amino acid count
    for (const std::string & aa : amino_acid_strings) {
      aa_mass_singletons.push_back({"mass_" + aa});
     
      igb.insert_dependency( TableDependency<std::string>(make_nonneg_uniform(aa, max_num_copies), p) );
    }
    
    //// Add Constant Multiplication Dependencies ////
    for (unsigned long i=0; i<amino_acid_strings.size(); ++i)
      igb.insert_dependency( ConstantMultiplierDependency<std::string>({amino_acid_strings[i]}, {aa_mass_singletons[i]}, {Peptide::masses[i]*mass_discretization}, false, true, DITHERING_SIGMA) );
    
    // Make additive dep. for total mass.
    LabeledPMF<std::string> total_mass = LabeledPMF<std::string>( {"total_mass"}, scaled_pmf_dither(PMF({1L},Tensor<double>({1ul},{1.0})), {mass_goal*mass_discretization}, DITHERING_SIGMA) );

    igb.insert_dependency( TableDependency<std::string>(total_mass, p) );
    igb.insert_dependency( AdditiveDependency<std::string>(aa_mass_singletons, {"total_mass"}, p) );
    
    // create inference graph
    _ig_ptr = new InferenceGraph<std::string>(igb.to_graph());

    write_graph_to_dot_file(*_ig_ptr, "mass_peptide_graph.dot");
  }

  ~MassPeptideSolver() {
    delete _ig_ptr;
  }

  void solve_and_print() {
    ///////////////////////
    ///// Solve Graph /////
    ///////////////////////
    
    //ig.print(std::cout);
    std::cout << "solving..." << std::endl;
    
    // apply message scheduler to inference graph
    _sched.add_ab_initio_edges(*_ig_ptr);
    
    // apply belief propagation to inference graph
    BeliefPropagationInferenceEngine<std::string> bpie(_sched, *_ig_ptr);

    std::vector<std::vector<std::string> > aa_singletons;
    for (char aa : Peptide::amino_acids)
      aa_singletons.push_back({std::string("")+aa});
    
    Clock c;
    c.tick();
    auto result = bpie.estimate_posteriors(aa_singletons);
    std::cout << "Time " << c.tock() << " in seconds" << std::endl;
    for (auto res : result)
      std::cout << res << std::endl;
  }

};

#endif
