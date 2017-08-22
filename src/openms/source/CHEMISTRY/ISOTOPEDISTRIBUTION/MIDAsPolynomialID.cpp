#include <boost/utility.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/MIDAsPolynomialID.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Polynomial.h>


#include <vector>
#include <utility>

#define DEBUG

using namespace std;
namespace OpenMS
{

  void addCounter(CounterSet& c, const double& abundance, const UInt& size, const UInt& N)
  {
    double expectation = size * abundance;
    double var = size * abundance *(1 - abundance);
    UInt U = expectation + (N * sqrt(1 + var));
    UInt B = expectation > (N * sqrt(1 + var)) ? ceil(expectation - (N * sqrt(1 + var))) : 0;
#ifdef DEBUG
    LOG_INFO << "Added counter with values " << B << " " << U <<endl;
#endif
    c.addCounter(B, U);
  }


  MIDAsPolynomialID::MIDAsPolynomialID(EmpiricalFormula& formula, double resolution):
    MIDAs(formula, resolution, 10),
    lighter_isotope(0)
  {
    for(EmpiricalFormula::const_iterator el = formula_.begin(); el != formula_.end(); ++el)
    {
      lighter_isotope += lightest_element(*(el->first)) * (el->second);
    }

    LOG_INFO << "Fine resolution: " << resolution_ << endl;

    mw_resolution = 1e-12;
    


  }

  inline double MIDAsPolynomialID::fact_ln(UInt x)
  {
    return boost::math::lgamma(x+1);
  }

  void MIDAsPolynomialID::run()
  {
    vector<Polynomial> el_dist;
    for(EmpiricalFormula::ConstIterator element = formula_.begin(); element != formula_.end(); ++element)
    {
      el_dist.push_back(generatePolynomial(*(element->first), element->second));
#ifdef DEBUG
      LOG_INFO << element->first->getName() <<" has " << el_dist.back().size() << " data points " << endl;
#endif
    }
    Polynomial& T = *(el_dist.begin());

    for(vector<Polynomial>::iterator pol = boost::next(el_dist.begin()); pol != el_dist.end(); ++pol)
    {
      multiplyPolynomials(T, *pol);
    }

#ifdef DEBUG
    LOG_INFO << "Distributian after multiplication has " << T.size() <<" elements" << endl;
#endif

    double probability = 0;
    for(Polynomial::const_iterator it = T.begin(); it != T.end(); ++it)
    {
      probability += it->getIntensity();
    }

    sort(T.begin(), T.end(), by_power);
    for(auto& pmember : T)
    {
      pmember.setMZ(pmember.getMZ() * mw_resolution);
    }
    distribution_.clear();
    move(T.begin(), T.end(), back_inserter(distribution_));

    merge(distribution_, 0.0001);

    LOG_INFO << "Lightest theoretical element " << lighter_isotope << endl;

    trimRight(0.0001);
    trimLeft(0.0001);

#ifdef DEBUG
    LOG_INFO << "Isotope Distribution of " << formula_.toString() << " successfully computed " << endl;
    LOG_INFO << "Isotope Distribution has " << T.size() << " data points " << endl;
#endif
  }

  MIDAsPolynomialID::Polynomial MIDAsPolynomialID::generatePolynomial(const Element& p, const SignedSize size)
  {
    vector<unsigned long> base_power;
    vector<double> log_prob;
    const IsotopeDistribution::ContainerType& isotope = p.getIsotopeDistribution().getContainer();
    CounterSet c(size);
    Polynomial pol;

    for(IsotopeDistribution::ConstIterator iso_it = isotope.begin(); iso_it != isotope.end(); ++iso_it)
    {
      if(iso_it->getIntensity() == 0)
      {
        continue;
      }
      addCounter(c, iso_it->getIntensity(), size, N);
      base_power.push_back(round(iso_it->getMZ() / mw_resolution));
      log_prob.push_back(log(iso_it->getIntensity()));
    }

    for(const CounterSet::ContainerType& counters = c.getCounters(); c.hasNext(); ++c)
    {
      Peak1D member;
      member.setMZ(0);
      member.setIntensity(fact_ln(size));
      UInt index = 0;
      for(CounterSet::ContainerType::const_iterator iso_count = counters.begin(); iso_count != counters.end(); ++iso_count, ++index)
      {
        member.setIntensity(member.getIntensity() + ((*iso_count) * log_prob[index]) - fact_ln((*iso_count)));
      }

      member.setIntensity(exp(member.getIntensity()));

      if(member.getIntensity() < min_prob)
      {
        continue;
      }
      // check if it is faster having another iteration
      index = 0;
      for(CounterSet::ContainerType::const_iterator iso_count = counters.begin(); iso_count != counters.end(); ++iso_count, ++index)
      {
        member.setMZ( member.getMZ() + (*iso_count)*base_power[index]);
      }

      pol.push_back(member);

    }

    return pol;
  }

  void MIDAsPolynomialID::multiplyPolynomials(Polynomial& f, Polynomial& g)
  {

    // sort polynomials to accelerate multiplication loop 
    sort(f.begin(), f.end(), desc_prob);
    sort(g.begin(), g.end(), desc_prob);

    double min_mass = min_element(g.begin(), g.end(), by_power)->getMZ()
                      + min_element(f.begin(), f.end(), by_power)->getMZ();
    double max_mass = max_element(g.begin(), g.end(), by_power)->getMZ()
                      + max_element(f.begin(), f.end(), by_power)->getMZ();
    double delta_mass = resolution_/mw_resolution;
    UInt size = round((max_mass-min_mass)/delta_mass);
    Polynomial fgid(size, Peak1D(0, 0));

    for(Polynomial::iterator g_it = g.begin(); g_it != g.end(); ++g_it)
    {
      for(Polynomial::iterator f_it = f.begin(); f_it != f.end(); ++f_it)
      {
        double prob = f_it->getIntensity() * g_it->getIntensity();
        if(prob > min_prob)
        {
          double mass = f_it->getMZ() + g_it->getMZ();
          UInt bin = round((mass - min_mass) / delta_mass);
          fgid[bin].setIntensity( fgid[bin].getIntensity() + prob);
          fgid[bin].setMZ( fgid[bin].getMZ() + mass * prob);
        }
        else
        {
          // Polynomials are sorted based on probability so we can safely break
          break;
        }
      }
    }

    fgid.erase(remove_if(fgid.begin(), fgid.end(), zero_prob), fgid.end());
    for(Polynomial::iterator f_it = fgid.begin(); f_it != fgid.end(); ++f_it)
    {
      f_it->setMZ(f_it->getMZ() / f_it->getIntensity());
    }
    f = fgid ;
  }

}
