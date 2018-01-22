

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  IsotopePatternGenerator::IsotopePatternGenerator(double probability_cutoff) : 
    IsotopeDistribution(),
    min_prob_(probability_cutoff)
  {

  }

  IsotopePatternGenerator::IsotopePatternGenerator() : 
    IsotopeDistribution(),
    min_prob_(1e-15)
  {

  }

  IsotopePatternGenerator::IsotopePatternGenerator(const IsotopeDistribution& rhs) :
    IsotopeDistribution(rhs)
  {
    
  }
  void IsotopePatternGenerator::merge(double resolution)
  {
    //raw must be ordered to work correctly ascending order on power field
    sortByMass();

    ContainerType raw = distribution_;
    double mass_range = (raw.back().getMZ() - raw.front().getMZ());
    UInt output_size = ceil(mass_range / resolution);
    distribution_.clear();
    ContainerType distribution(output_size, Peak1D(0, 0));
    double delta = mass_range / output_size;

    for (auto & p : raw)
    {
      UInt index = round((p.getMZ() - raw.front().getMZ()) / resolution);
      if (index >= distribution.size()) { continue; }
      double mass = raw.front().getMZ() + (index * delta);
      distribution[index].setMZ(mass);
      distribution[index].setIntensity(distribution[index].getIntensity() + p.getIntensity());
    }
    distribution_ = distribution;
    trimIntensities(min_prob_);
  }

  
  /* Start of the midas interface */
  MIDAs::MIDAs(double resolution, double min_prob, UInt N_):
    IsotopePatternGenerator(min_prob),
    N(N_),
    resolution_(resolution)
  {    
  }

  MIDAs::MIDAs() : IsotopePatternGenerator()
  {
  }

  MIDAs::MIDAs(const IsotopeDistribution& isotope_distribution) : 
    IsotopePatternGenerator(isotope_distribution)
  {
  }


}
