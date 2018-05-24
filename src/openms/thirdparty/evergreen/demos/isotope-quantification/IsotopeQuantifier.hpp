#ifndef _ISOTOPEQUANTIFIER_HPP
#define _ISOTOPEQUANTIFIER_HPP

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../../Evergreen/evergreen.hpp"

#include "Elements.hpp"
#include "../../Utility/inference_utilities.hpp"
#include "../../Utility/to_string.hpp"
#include "../../Utility/Clock.hpp"
#include "../../Utility/L1Regularization.hpp"

#include "../../Utility/graph_to_dot.hpp"
#include <fstream>

// To consider missing peaks, insert them into the spectra as values
// with small or zero intensity.

class IsotopeQuantifier {
private:
  const Elements _elements;
  const unsigned int _prior_maximum_copies_of_element;
  const unsigned int _maximum_number_unique_elements;
  const unsigned int _intensity_discretization;
  const double _sigma_observed_intensities;

  double _mass_discretization;

  static constexpr double DITHERING_SIGMA = 0.1;
  // The value beyond which Gaussian tails are no longer considered:
  static constexpr double GAUSSIAN_TAIL_EPSILON = 1e-32;

  std::map<double, std::vector<Isotope> > _theoretical_peaks_to_isotopes;
  bool _include_unobserved_peaks;

  // observed:
  std::map<double, double> _observed_peak_masses_to_intensities;
  
  std::set<std::string> _used_elements;
  std::set<std::string> _used_isotopes;
  
  Scheduler<std::string> & _scheduler;
  InferenceGraph<std::string>* _ig_ptr;

  static const std::string intensity_prefix;

  void load_peaks_from_file_and_discretize(const std::string & peak_file){
    std::ifstream fin(peak_file);
    assert(fin.is_open() && "Error: File not found");

    std::string garbage;
    fin >> garbage;
    assert(garbage == "mass_discretization");
    fin >> _mass_discretization;
    
    std::string line;
    double mass;
    double intensity;
    
    while ( fin >> mass >> intensity ) {
      if (_observed_peak_masses_to_intensities.find(mass) == _observed_peak_masses_to_intensities.end())
	_observed_peak_masses_to_intensities[mass] = 0.0;

      _observed_peak_masses_to_intensities[mass] += intensity;
    }
    fin.close();

    _observed_peak_masses_to_intensities = mass_discretized_peaks(_observed_peak_masses_to_intensities, _mass_discretization, _include_unobserved_peaks);
  }
  
  void map_observed_peaks_to_isotopes_with_similar_mass() {
    for (const std::pair<std::string, std::vector<Isotope> > & ele: _elements) {
      for (const Isotope & iso: ele.second) {
	const double discretized_mass = round(iso.mass * _mass_discretization) / _mass_discretization;

	auto iter = _observed_peak_masses_to_intensities.find(discretized_mass);
	// If there is an observed peak at this discretized_mass:
	if (iter != _observed_peak_masses_to_intensities.end()) {
	  // theoretical mass for isotope matches an observed mass
	  _theoretical_peaks_to_isotopes[discretized_mass].push_back(iso);
	  _used_isotopes.insert(iso.name + " " + to_string(iso.mass));
	  _used_elements.insert(ele.first);
	}
      }
    }
  }

  void add_regularization(InferenceGraphBuilder<std::string> & igb, double p) {
    LabeledPMF<std::string> sum_of_indicators = make_nonneg_uniform<std::string>("SumOfIndicators", _maximum_number_unique_elements);
    
    std::vector<std::string> indicators_for_used_elements(_used_elements.size());
    std::vector<std::string> used_elements_vector(_used_elements.begin(), _used_elements.end());
    
    for (unsigned long i=0; i<indicators_for_used_elements.size(); ++i)
      indicators_for_used_elements[i] = "Indicator[ " + used_elements_vector[i] + ">0 ]";
    
    L1Regularization<std::string>::apply(igb, used_elements_vector, indicators_for_used_elements, sum_of_indicators, p, _prior_maximum_copies_of_element);
  }

  void print_isotopes_matching_observed_peaks() {
    std::cout << "discretized data & matching isotopes" << std::endl;
    for (auto pr : _observed_peak_masses_to_intensities) {
      std::cout << pr.first << " " << pr.second << " ";
      
      auto iter = _theoretical_peaks_to_isotopes.find(pr.first);
      if (iter != _theoretical_peaks_to_isotopes.end()) {
	const std::vector<Isotope> & matching_isos = iter->second;
	for (const Isotope & iso : matching_isos) {
	  std::cout << iso << " ";
	}
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  void add_constant_multipliers(InferenceGraphBuilder<std::string> & igb) {
    // Make constant multiplier dependencies that say isotope
    // abundance is some constant times the element abundance.
    std::set<Isotope> isotopes_matching_any_observed;
    for (const std::pair<double, std::vector<Isotope> > & peak_and_isotopes: _theoretical_peaks_to_isotopes)
      for (const Isotope & iso: peak_and_isotopes.second)
	isotopes_matching_any_observed.insert(iso);

    for (const Isotope & iso : isotopes_matching_any_observed) {
      std::string isotope_id = intensity_prefix + iso.name + " " + to_string(iso.mass);
      // false, true --> when multiplying don't interpolate (since
      // we're starting with counts), but interpolate when dividing:
      igb.insert_dependency( ConstantMultiplierDependency<std::string>({iso.name}, {isotope_id}, {iso.abundance * _intensity_discretization}, false, true, DITHERING_SIGMA) );
    }
  }

  void add_gaussians_for_observed_peaks(InferenceGraphBuilder<std::string> & igb, double p) {
    // Make table dependency for intensity of each peak_i, where
    // intensity is a nonnegative gaussian distribution with
    // mean=observed intensity and standard
    // deviation=_sigma_observed_intensities.
    for (const std::pair<double, double> & peak: _observed_peak_masses_to_intensities ) {
      double observed_mass = peak.first;
      std::string peak_var = intensity_prefix + "peak" + to_string(observed_mass);
      
      double pre_discretized_observed_intensity = peak.second * _intensity_discretization;
      auto nonneg_gaussian_for_peak = make_nonneg_pseudo_gaussian(peak_var, pre_discretized_observed_intensity, _sigma_observed_intensities, GAUSSIAN_TAIL_EPSILON, long(pre_discretized_observed_intensity*10), 1e-5);
      igb.insert_dependency( TableDependency<std::string>(nonneg_gaussian_for_peak, p));
    }
  }

  void add_additive_dependencies(InferenceGraphBuilder<std::string> & igb, double p) {
    // Make additive dep. for intensity of peak_i (it should equal the
    // sum of the quantities of the element isotopes matching it).
    for (const std::pair<double, std::vector<Isotope> > & peak : _theoretical_peaks_to_isotopes) {
      double observed_mass = peak.first;
      std::string peak_var = intensity_prefix + "peak" + to_string(observed_mass);

      std::vector<std::vector<std::string> > isotopes_that_sum_to_this_peak;
      for(const Isotope & responsible_iso : peak.second) {
        assert(peak.second.size() != 0 && "Observed peak did not match any theoretical element isotope peaks");

        isotopes_that_sum_to_this_peak.push_back({ intensity_prefix + responsible_iso.name + " " + to_string(responsible_iso.mass) });
      }
      igb.insert_dependency( AdditiveDependency<std::string>(isotopes_that_sum_to_this_peak, {peak_var}, p) );
    }
  }

  void build_graph(const double p) {
    BetheInferenceGraphBuilder<std::string> igb;

    // Add uniform priors for each candidate element:
    for (const std::string el : _used_elements)
      igb.insert_dependency( TableDependency<std::string>(make_nonneg_uniform(el, _prior_maximum_copies_of_element), p) );
    
    // Add regularization if it is used:
    if (_maximum_number_unique_elements != 0)
      add_regularization(igb, p);

    add_constant_multipliers(igb);
      
    add_gaussians_for_observed_peaks(igb, p);

    add_additive_dependencies(igb, p);

    // Create inference graph from the graph builder:
    _ig_ptr = new InferenceGraph<std::string>(igb.to_graph());

    write_graph_to_dot_file(*_ig_ptr, "isotope_graph.dot");
  }

public:
  // Default value of _maximum_number_unique_elements=0 --> don't use regularization. 
  IsotopeQuantifier(const std::string & peak_file, const Elements & ele, Scheduler<std::string> & scheduler, const double p, unsigned long intensity_discretization, const double standard_deviation_observed_intensities, unsigned long prior_maximum_copies_of_element, bool include_unobserved_peaks, unsigned long maximum_number_unique_elements=0):
    _elements(ele),
    _prior_maximum_copies_of_element(prior_maximum_copies_of_element),
    _maximum_number_unique_elements(maximum_number_unique_elements),
    _intensity_discretization(intensity_discretization),
    _sigma_observed_intensities(standard_deviation_observed_intensities*intensity_discretization),
    _include_unobserved_peaks(include_unobserved_peaks),
    _scheduler(scheduler)
  {
    load_peaks_from_file_and_discretize(peak_file);
    map_observed_peaks_to_isotopes_with_similar_mass();

    build_graph(p);
    print_isotopes_matching_observed_peaks();
  }

  static std::map<double, double> theoretical_peaks_from_chemical_formula(const std::map<std::string, unsigned int> & formula, const Elements & element_collection) {
    std::map<double, double> result;
    for (const std::pair<std::string, unsigned int> & element: formula) {
      assert(element.second != 0 && "Error: Element count must be >0");
      
      for(const Isotope & iso: element_collection.get(element.first) ) {
	auto iter = result.find(iso.mass);
	// Just in case two values have identical masses:
	if (iter == result.end())
	  result[iso.mass] = 0.0;

	result[iso.mass] += iso.abundance*element.second;
      }
    }
    return result;
  }

  static std::map<double, double> mass_discretized_peaks(const std::map<double, double> & exact, double mass_discretization, bool include_unobserved_peaks) {
    // Get the maximum by using the fact that map is sorted ascending (add 1 because of 0 bin):
    std::vector<double> pre_result( (unsigned long)ceil(exact.rbegin()->first * mass_discretization) + 1, 0.0 );

    for (const std::pair<double, double> & mass_and_intensity : exact) {
      const double mass = mass_and_intensity.first;
      const double intensity = mass_and_intensity.second;

      const long discretized_mass = round(mass*mass_discretization);
      pre_result[discretized_mass] += intensity;
    }

    std::map<double, double> result;
    for (unsigned long i=0; i<pre_result.size(); ++i) {
      if (pre_result[i] > 0.0 || include_unobserved_peaks) {
	double mass = double(i) / mass_discretization;

	// add to result map
	if (result.find(mass) == result.end())
	  result[mass] = 0.0;
	result[mass] += pre_result[i];
      }
    }
    
    return result;
  }

  // mass_discretization = 100 means that accuracy is to 1/100 dalton
  // (pre rounding).
  static std::map<double, double> mass_discretized_theoretical_peaks_from_chemical_formula(const std::map<std::string, unsigned int> & formula, const Elements & element_collection, double mass_discretization, bool include_unobserved_peaks) {
    std::map<double, double> exact = theoretical_peaks_from_chemical_formula(formula, element_collection);
    return mass_discretized_peaks(exact, mass_discretization, include_unobserved_peaks);
  }

  void run_and_print_results() {
    // apply message scheduler to inference graph
    _scheduler.add_ab_initio_edges(*_ig_ptr);
    
    // apply belief propagation to inference graph
    BeliefPropagationInferenceEngine<std::string> bpie(_scheduler, *_ig_ptr);
    
    Clock c;
    c.tick();

    std::vector<std::vector<std::string> > element_singletons;
    for (const std::string & el : _used_elements)
      element_singletons.push_back( {el} );

    auto result = bpie.estimate_posteriors(element_singletons);

    std::cout << "Time " << c.tock() << " in seconds" << std::endl;
    for (auto res : result)
      std::cout << res << std::endl;
    
    std::cout << "Elements matching no observed peaks (treat as having 0 abundance with probability ~1):" << std::endl;
    for (const std::pair<std::string, std::vector<Isotope> > & ele: _elements ) {
      if ( _used_elements.find(ele.first) == _used_elements.end() ){
        std::cout << ele.first << " " << PMF({0L}, Tensor<double>({1ul},{1.0})) << std::endl;
      }
    }
  }
};

const std::string IsotopeQuantifier::intensity_prefix = "intensity ";

#endif
