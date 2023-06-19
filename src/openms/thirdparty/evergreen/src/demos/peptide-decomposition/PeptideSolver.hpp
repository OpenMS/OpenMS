#ifndef _PEPTIESOLVER_HPP
#define _PEPTIESOLVER_HPP

#include "../../Evergreen/evergreen.hpp"

#include "../../Utility/inference_utilities.hpp"
#include "../../Utility/Clock.hpp"

#include "Peptide.hpp"

#include "../../Utility/graph_to_dot.hpp"

class PeptideSolver {
private:
  Scheduler<std::string> & _sched;
  InferenceGraph<std::string> *_ig_ptr;

  static constexpr double DITHERING_SIGMA = 0.1;

  // When goal mass or goal hydrophobicity is non-integral, distribute
  // mass equally over both adjacent bins.
  static constexpr double DITHERING_SIGMA_GOALS = 10000.0;

public:
  
  // Note: This could easily use total mass / amino acid mass (+ some
  // small amount for stability) to make a custom maximum number of
  // copies for each amino acid. Same could be done with hydrophobicity,
  // and the minimum of both maxes could be used:
  PeptideSolver(double mass_goal, double hydrophobicity_goal, const double & p, const unsigned int max_num_copies, const double mass_discretization, const double hydrophobicity_discretization, Scheduler<std::string> & sched):
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
    std::vector<std::vector<std::string> > aa_hydrophobicity_singletons;

    //// Add Table Dependencies ////
    // Make uniform distribution for each amino acid count
    for (const std::string & aa : amino_acid_strings) {
      // Note: max_num_copies could be inferred for each amino acid by
      // dividing goal mass (or hydrophobicity) by the mass of each
      // amino acid; but to start with, make it simple:
      aa_mass_singletons.push_back({"mass " + aa});
      aa_hydrophobicity_singletons.push_back({"hydrophobicity " + aa});
     
      igb.insert_dependency( TableDependency<std::string>(make_nonneg_uniform(aa, max_num_copies), p) );
    }
    
    //// Add Constant Multiplication Dependencies ////
    // For each amino acid, make constant mult. dep. for both mass and hydrophobicity.
    for (unsigned long i=0; i<amino_acid_strings.size(); ++i) {
      igb.insert_dependency( ConstantMultiplierDependency<std::string>({amino_acid_strings[i]}, {aa_mass_singletons[i]}, {Peptide::masses[i]*mass_discretization}, false, true, DITHERING_SIGMA) );
      igb.insert_dependency( ConstantMultiplierDependency<std::string>({amino_acid_strings[i]}, {aa_hydrophobicity_singletons[i]}, {Peptide::hydrophobicities[i]*hydrophobicity_discretization}, false, true, DITHERING_SIGMA) );
    }
    
    // Make additive dep. for total mass.
    LabeledPMF<std::string> total_mass = LabeledPMF<std::string>( {"total_mass"}, scaled_pmf_dither(PMF({1L},Tensor<double>({1ul},{1.0})), {mass_goal*mass_discretization}, DITHERING_SIGMA_GOALS) );

    igb.insert_dependency( TableDependency<std::string>(total_mass, p) );
    igb.insert_dependency( AdditiveDependency<std::string>(aa_mass_singletons, {"total_mass"}, p) );
    
    // Make additive dep. for total hydrophobicity.
    LabeledPMF<std::string> total_hydrophobicity = LabeledPMF<std::string>( {"total_hydrophobicity"}, scaled_pmf_dither(PMF({1L},Tensor<double>({1ul},{1.0})), {hydrophobicity_goal*hydrophobicity_discretization}, DITHERING_SIGMA_GOALS) );
    igb.insert_dependency( TableDependency<std::string>(total_hydrophobicity, p) );
    igb.insert_dependency( AdditiveDependency<std::string>(aa_hydrophobicity_singletons, {"total_hydrophobicity"}, p) );

    // create inference graph
    _ig_ptr = new InferenceGraph<std::string>(igb.to_graph());

    write_graph_to_dot_file(*_ig_ptr, "peptide_graph.dot");
  }

  ~PeptideSolver() {
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
