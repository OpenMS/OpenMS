#ifndef _INFERENCE_UTILITIES_HPP
#define _INFERENCE_UTILITIES_HPP

#include "Clock.hpp"
#include <iostream>
#include <algorithm>

template<template <typename, typename...> class CONTAINER, typename T, typename ...OTHER_ARGS>
std::vector<std::vector<T> > make_singletons(const CONTAINER<T, OTHER_ARGS...> & var_container) {
  std::vector<std::vector<T> > result;
  for (const T & t : var_container)
    result.push_back({t});
  return result;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T> > & container) {
  std::vector<T> result;
  for (const std::vector<T> & row : container)
    for (const T & var : row)
      result.push_back(var);
  return result;
}

template <typename T>
std::vector<T> from_singletons(const std::vector<std::vector<T> > & singletons) {
  for (const std::vector<T> & t : singletons)
    assert( t.size() == 1 );

  return flatten(singletons);
}

template <typename T>
void estimate_and_print_posteriors(BruteForceInferenceEngine<T> & bf, const std::vector<std::vector<T> > & singletons){
  Clock c;
  auto result = bf.estimate_posteriors(singletons);
  c.ptock();
  for (auto res : result)
    std::cout << res << std::endl;
}

template <typename T>
void estimate_and_print_posteriors(BeliefPropagationInferenceEngine<T> & bpie, const std::vector<std::vector<T> > & singletons){
  Clock c;
  auto result = bpie.estimate_posteriors(singletons);
  c.ptock();
  for (auto res : result)
    std::cout << res << std::endl;
}

template <typename T>
LabeledPMF<T> make_nonneg_uniform(const T & var_name, unsigned long max_val) {
  Tensor<double> ten({max_val+1});
  ten.flat().fill(1.0);
  PMF pmf( {0L}, ten );
  return LabeledPMF<T>({var_name}, pmf);
}

inline double gaussian_density(double x, const double mu, const double sigma){
  double var = sigma*sigma;
  double dev = (x - mu);
  return ( exp(-(dev*dev) / (2*var)) / sigma ) / (sigma * sqrt(2*M_PI));
}

inline double inverse_standard_norm_cdf(const double p) {
  assert(p >= 0.5);
  return 5.5556*(1-pow((1-p)/p, 0.1186));
}

template <typename T>
TableDependency<T> table_dependency_by_gaussian(const T & label, const Vector<long> & support, double goal, const double p, const double sd){
  Tensor<double> pmf({support.size()});
  for(unsigned int i=0; i<support.size(); ++i)
    pmf[i] = gaussian_density(support[i], goal, sd);
  TableDependency<T> td(LabeledPMF<T>({label}, PMF({(long)floor(support[0])}, pmf)), p);
  return td;
}


template <typename T>
LabeledPMF<T> make_gaussian(const T & label, double mu, const double sigma, const double epsilon){
  const double max_z = inverse_standard_norm_cdf(1-epsilon);
  const double min_z = -max_z;

  // z-score = (x - mu) / sigma
  // --> x = z*sigma + mu

  // Find minimum and maximum integer values beyond which tails of
  // Gaussian are < epsilon.
  double min_double_support = mu + min_z*sigma;
  double max_double_support = mu + max_z*sigma;

  long min_support = (long) floor(min_double_support);
  long max_support = (long) ceil(max_double_support);

  Tensor<double> table({ (unsigned long)(max_support - min_support + 1) });
  for(unsigned long i=0; i<table.flat_size(); ++i)
    table[i] = gaussian_density(min_support + long(i), mu, sigma);

  return LabeledPMF<T>( {label}, PMF({min_support}, table) );
}

template <typename T>
LabeledPMF<T> make_nonneg_gaussian(const T & label, double mu, const double sigma, const double epsilon){
  const double max_z = inverse_standard_norm_cdf(1-epsilon);
  const double min_z = -max_z;

  // z-score = (x - mu) / sigma
  // --> x = z*sigma + mu

  // Find minimum and maximum integer values beyond which tails of
  // Gaussian are < epsilon.
  double min_double_support = mu + min_z*sigma;
  double max_double_support = mu + max_z*sigma;

  long min_support = (long) floor(min_double_support);
  long max_support = (long) ceil(max_double_support);

  min_support = std::max(0L, min_support);
  assert(max_support >= min_support);

  Tensor<double> table({ (unsigned long)(max_support - min_support + 1) });
  for(unsigned long i=0; i<table.flat_size(); ++i)
    table[i] = gaussian_density(min_support + long(i), mu, sigma);
  return LabeledPMF<T>({label}, PMF({min_support}, table));
}

// Like a Gaussian but with guaranteed minimum probability (to
// increase the density of tails):
template <typename T>
LabeledPMF<T> make_nonneg_pseudo_gaussian(const T & label, double mu, const double sigma, const double epsilon, long max_support, double pseudo_count){
  const double max_z = inverse_standard_norm_cdf(1-epsilon);

  // z-score = (x - mu) / sigma
  // --> x = z*sigma + mu

  // Find minimum and maximum integer values beyond which tails of
  // Gaussian are < epsilon.
  double max_double_support = mu + max_z*sigma;

  long min_support = 0;
  max_support = std::max( (long) ceil(max_double_support), max_support);

  assert(max_support >= min_support);

  Tensor<double> table({ (unsigned long)(max_support - min_support + 1) });
  for(unsigned long i=0; i<table.flat_size(); ++i)
    table[i] = std::max( gaussian_density(min_support + long(i), mu, sigma), pseudo_count);
  return LabeledPMF<T>({label}, PMF({min_support}, table));
}

template <typename T>
LabeledPMF<T> make_bernoulli(const T & var_name, double probability_x_equals_one) {
  assert(probability_x_equals_one >= 0.0 && probability_x_equals_one <= 1.0 && "make_bernoulli must receive a valid probability");
  return LabeledPMF<T>({var_name}, PMF({0L}, Tensor<double>({2ul}, {1-probability_x_equals_one, probability_x_equals_one})));
}

#endif
