#ifndef _HYDROPHOBICITYPEPTIDESOLVER_HPP
#define _HYDROPHOBICITYPEPTIDESOLVER_HPP

#include "../../Evergreen/evergreen.hpp"

#include "../../Utility/inference_utilities.hpp"
#include "../../Utility/Clock.hpp"

#include "Peptide.hpp"

#include "../../Utility/graph_to_dot.hpp"

class HydrophobicityPeptideSolver {
private:
  Scheduler<std::string> & _sched;
  InferenceGraph<std::string> *_ig_ptr;

  static constexpr double DITHERING_SIGMA = 0.1;
  // The value beyond which Gaussian tails are no longer considered:
  static constexpr double GAUSSIAN_TAIL_EPSILON = 0.005;

public:
  
  HydrophobicityPeptideSolver(double hydrophobicity_goal, const double & p, const unsigned int max_num_copies, const double hydrophobicity_discretization, Scheduler<std::string> & sched):
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
    std::vector<std::vector<std::string> > aa_hydrophobicity_singletons;

    //// Add Table Dependencies ////
    // Make uniform distribution for each amino acid count
    for (const std::string & aa : amino_acid_strings) {
      aa_hydrophobicity_singletons.push_back({"hydrophobicity_" + aa});
     
      igb.insert_dependency( TableDependency<std::string>(make_nonneg_uniform(aa, max_num_copies), p) );
    }
    
    //// Add Constant Multiplication Dependencies ////
    for (unsigned long i=0; i<amino_acid_strings.size(); ++i) {
      igb.insert_dependency( ConstantMultiplierDependency<std::string>({amino_acid_strings[i]}, {aa_hydrophobicity_singletons[i]}, {Peptide::hydrophobicities[i]*hydrophobicity_discretization}, false, true, DITHERING_SIGMA) );
    }
    
    // Make additive dep. for total hydrophobicity.
    LabeledPMF<std::string> total_hydrophobicity = LabeledPMF<std::string>( {"total_hydrophobicity"}, scaled_pmf_dither(PMF({1L},Tensor<double>({1ul},{1.0})), {hydrophobicity_goal*hydrophobicity_discretization}, DITHERING_SIGMA) );
    igb.insert_dependency( TableDependency<std::string>(total_hydrophobicity, p) );
    igb.insert_dependency( AdditiveDependency<std::string>(aa_hydrophobicity_singletons, {"total_hydrophobicity"}, p) );

    // create inference graph
    _ig_ptr = new InferenceGraph<std::string>(igb.to_graph());

    write_graph_to_dot_file(*_ig_ptr, "hydro_peptide_graph.dot");
  }

  ~HydrophobicityPeptideSolver() {
    delete _ig_ptr;
  }

  void solve_and_print() {
    ///////////////////////
    ///// Solve Graph /////
    ///////////////////////
    
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
