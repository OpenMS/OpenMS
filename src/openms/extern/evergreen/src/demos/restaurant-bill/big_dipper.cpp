#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/inference_utilities.hpp"
#include <fstream>

class BigDipperIceCream {
private:
  static std::set<double> load_prices(const std::string & menu_filename) {
    std::set<double> result;
    
    std::ifstream fin(menu_filename);
    std::string item_name;
    double price;
    while (fin >> item_name >> price)
      result.insert(price);
    
    return result;
  }
  
  static std::set<unsigned int> load_prices_in_quarters(const std::string & menu_filename) {
    std::set<double> prices = load_prices(menu_filename);
    std::set<unsigned int> result;
    for (double price : prices)
      // Prices are all divisible by 0.25-- thanks Big Dipper!
      // Regardless, round just to be safe (the value 0.99999 would
      // cast to integer 0).
      result.insert( (unsigned int)round(price / 0.25) );
    return result;
  }

  std::vector<unsigned int> _prices_in_quarters;
  
public:
  BigDipperIceCream(const std::string & menu_filename) {
    std::set<unsigned int> price_set = load_prices_in_quarters(menu_filename);
    _prices_in_quarters = std::vector<unsigned int>(price_set.begin(), price_set.end());
    std::cout << "K=" << _prices_in_quarters[_prices_in_quarters.size()-1] << std::endl;
  }

  PMF generate_pmf_of_preferences() {
    // Distribution will be in {0, 1, ... maximum price}. Use sorted
    // order of set to get maximum value and add 1:
    Tensor<double> probability_table( {*_prices_in_quarters.rbegin()+1ul} );

    for (unsigned int price : _prices_in_quarters) {
      // Choose a probability that the person buys this item (note: it
      // is not yet a true probability, since we do not know if it
      // sums to 1 with the other items, but that will be normalized
      // in the PMF constructor).
      double prob = rand() % 10000 / 9999.0 + 0.1;

      probability_table[price] = prob;
    }

    return PMF({0L}, probability_table);
  }
};

unsigned int randomly_sample_from_1d_pmf(const PMF & pmf) {
  double uniform = rand() % 10000 / 9999.0;

  double cumulative = 0.0;
  for (unsigned long i=0; i<pmf.table().flat_size(); ++i) {
    cumulative += pmf.table()[i];

    if (cumulative >= uniform)
      return i + pmf.first_support()[0];
  }

  // Should be impossible (sum of masses should = 1.0), but just in
  // case:
  return pmf.last_support()[0];
}

int main(int argc, char**argv) {
  if (argc != 3) {
    std::cerr << "Usage: bill_solver <N> <p>" << std::endl;
    exit(1);
  }

  const unsigned long N = atoi(argv[1]);
  const double p = atof(argv[2]);

  BetheInferenceGraphBuilder<std::string> igb;
  
  /*
    Prices from
    Big Dipper Ice Cream
    631S Higgins Ave.
    Missoula Montana
  */
  BigDipperIceCream bdic("big-dipper-prices.txt");

  unsigned long total_spent_in_quarters = 0;
  for (unsigned long i=0; i<N; ++i) {
    PMF pmf = bdic.generate_pmf_of_preferences();
    unsigned int person_spent = randomly_sample_from_1d_pmf(pmf);
    total_spent_in_quarters += person_spent;

    LabeledPMF<std::string> lpmf( {"X_" + to_string(i)}, pmf );
    igb.insert_dependency( TableDependency<std::string>(lpmf, p) );

    std::cout << lpmf << " " << person_spent << std::endl;
  }
  // We know that Y = total_spent_in_quarters with 100% probability:
  igb.insert_dependency( TableDependency<std::string>(LabeledPMF<std::string>({"Y"}, PMF({long(total_spent_in_quarters)}, Tensor<double>({1ul},{1.0}))), p) );

  // We know that Y = X_0 + X_1 + ... + X_{n-1}
  std::vector<std::vector<std::string> > input_singletons;
  for (unsigned long i=0; i<N; ++i)
    input_singletons.push_back( {"X_" + to_string(i)} );
  igb.insert_dependency( AdditiveDependency<std::string>(input_singletons, {"Y"}, p) );

  InferenceGraph<std::string> ig = igb.to_graph();

  FIFOScheduler<std::string> sched(0.0, 1e-8, N*8ul);
  sched.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<std::string> bpie(sched, ig);

  estimate_and_print_posteriors(bpie, {{"X_0"}});
  
  return 0;
}
