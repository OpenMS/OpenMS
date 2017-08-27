

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  IsotopePatternGenerator::IsotopePatternGenerator() : IsotopeDistribution()
  {}

  IsotopePatternGenerator::IsotopePatternGenerator(const IsotopeDistribution& rhs) :
    IsotopeDistribution(rhs)
  {
    
  }
  void IsotopePatternGenerator::merge(ContainerType& raw, double resolution)
  {
    //raw must be ordered to work correctly ascending order on power field

    UInt output_size = ceil((raw.back().getMZ() - raw.front().getMZ())/resolution);

    distribution_.clear();
    distribution_.resize(output_size, Peak1D(0, 0));

    for(auto& p : raw)
    {
      // Is this the case?
      
      UInt index = round((p.getMZ() - raw.front().getMZ())/resolution);
      if(index >= distribution_.size()){

        LOG_INFO << index <<endl;
        continue;
      }
      auto& mass = distribution_[index].getPosition() ;

      mass = mass == 0 ?
             raw.front().getMZ() * index :
             mass;
      distribution_[index].setIntensity(distribution_[index].getIntensity() + p.getIntensity());
    }
  }

  void IsotopePatternGenerator::dumpIDToFile(String file)
  {
    ofstream out(file.c_str());
    for(auto& sample : distribution_)
    {
      out << sample.getMZ() << " "<<sample.getIntensity() << endl;
    }
    
    out.close();
  }

  
  /* Start of the midas interface */
  MIDAs::MIDAs(double resolution, UInt N_):
    IsotopePatternGenerator(),
    min_prob(1e-16),
    resolution_(resolution),
    N(N_)
  {
    
  }

  MIDAs::MIDAs() : IsotopePatternGenerator()
  {
  }

  //MIDAs::MIDAs(const IsotopeDistribution& isotope_distribution) : 
  //  IsotopePatternGenerator(isotope_distribution)
  //{
  //}


}
