
#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEPATTERNSOLVER_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEPATTERNSOLVER_H



#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/CONCEPT/Types.h>

#include <deque>
#include <algorithm>

namespace OpenMS
{
  class OPENMS_DLLAPI IsotopePatternGenerator : public IsotopeDistribution
  {
 public:
    IsotopePatternGenerator();
    IsotopePatternGenerator(double probability_cutoff);
    IsotopePatternGenerator(const IsotopeDistribution&);
    virtual void run(const EmpiricalFormula&) = 0;
    void merge(double);
 protected:
    double min_prob_;
    
  };

  class OPENMS_DLLAPI MIDAs : public IsotopePatternGenerator
  {
 public:
    
    typedef std::deque<Peak1D> Polynomial;
    MIDAs(double resolution, double probability_cutoff, UInt N);
    MIDAs();
    MIDAs(const IsotopeDistribution& isotope_distribution);
 protected:
    UInt N;
    double resolution_;
  };

  inline bool desc_prob(const Peak1D& p0, const Peak1D& p)
  {
    return p0.getIntensity() > p.getIntensity();
  }

  inline bool by_power(const Peak1D& p0, const Peak1D& p)
  {
    return p0.getMZ() < p.getMZ();
  }

  inline bool zero_prob(const Peak1D& m)
  {
    return m.getIntensity() == 0;
  }

  inline bool zero_power(const Peak1D& m)
  {
    return m.getMZ() == 0;
  }

  inline bool lightest(const IsotopeDistribution::MassAbundance& a, 
                       const IsotopeDistribution::MassAbundance& b)
  {
    return a.getMZ() < b.getMZ();
  }

  inline double lightest_element(const Element& el)
  {
    return min_element(el.getIsotopeDistribution().begin(), 
                       el.getIsotopeDistribution().end(), lightest)->getMZ();
  }


}

#endif
