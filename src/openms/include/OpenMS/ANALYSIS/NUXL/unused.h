
//#define ANNOTATED_QUANTILES 1
#ifdef ANNOTATED_QUANTILES

#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
using namespace boost::accumulators;

typedef accumulator_set<double, stats<tag::p_square_quantile> > quantile_accu_t;
typedef accumulator_set<double, stats<tag::extended_p_square_quantile(quadratic)> > accumulator_t_quadratic;
          
struct SpectrumLevelScoreQuantiles
{
  SpectrumLevelScoreQuantiles():
    acc_(extended_p_square_probabilities = std::vector<double>{ 0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999, 1.00 })
  {   
  }
  void insert(double v) { acc_(v); }
  
  double quantileOfValue(double v)
  {
    double l = 0.0;
    double h = 1.0;
    double mid = l + (h - l) / 2.0; // we start with the median (0.5)
    double p_value = quantile(acc_, quantile_probability = mid); // value of quantile p

    size_t iter(0);
    while (fabs(p_value - v) > 0.01 && iter < 100)
    {
      mid = l + (h - l) / 2.0;
      p_value = quantile(acc_, quantile_probability = mid); // value of quantile p (e.g., 1234.56)
      if (p_value > v) // if the current quantile value (e.g., of the median) has a value larger than our value of interest
      {
        h = mid;  // then we need to search in the lower quantile range
      }
      else 
      {
        l = mid;
      }
      ++iter;
    }
    return l;
  }
  private:
    accumulator_t_quadratic acc_;
};

#endif


