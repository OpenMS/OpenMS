#include <string>
#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/inference_utilities.hpp"

// A simple demo of brute force inference

// Problem explained in
// https://en.wikipedia.org/wiki/Prisoner's_dilemma
int main() {
  const double p = 2.0;

  //////////////////////////////////
  ///// Construct Dependencies /////
  //////////////////////////////////
  // prior distribution of person1
  TableDependency<std::string> td1(LabeledPMF<std::string>({"person1"}, PMF({0L}, Tensor<double>({2ul}, {0.8, 0.2}))), p);
  
  // prior distribution of person2
  TableDependency<std::string> td2(LabeledPMF<std::string>({"person2"}, PMF({0L}, Tensor<double>({2ul}, {0.2, 0.8}))), p);
  
  // conditional dependency of person1 and pearson2

  TableDependency<std::string> td3(LabeledPMF<std::string>({"person1", "person2"}, PMF({0L,0L}, Tensor<double>({2ul,2ul}, {.87, .13, .74, .26}))), p);
  
  ///////////////////////
  ///// Solve Graph /////
  ///////////////////////
  
  BruteForceInferenceEngine<std::string> bf({td1, td2, td3},p);
  
  Clock c;
  std::vector<LabeledPMF<std::string> > result = bf.estimate_posteriors({{"person1"}, {"person2"}});
  std::cout << "BF Time: " << c.tock() << " in seconds" << std::endl;
  for (auto res : result)
    std::cout << res << std::endl;
  std::cout << std::endl;

  return 0;
}
