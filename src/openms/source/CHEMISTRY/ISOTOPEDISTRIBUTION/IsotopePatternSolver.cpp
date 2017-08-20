
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternSolver.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  void IsotopePatternSolver::merge(ContainerType& raw, double resolution)
  {
    //raw must be ordered to work correctly ascending order on power field

    UInt output_size = ceil((raw.back().getMZ() - raw.front().getMZ())/resolution);
    LOG_INFO << "output size " << output_size << endl;
    LOG_INFO << "raw size " << raw.size() <<endl;

    distribution_.clear();
    distribution_.resize(output_size, Peak1D(0, 0));

    for(auto& p : raw)
    {
      // Is this the case?
      
      UInt index = round((p.getMZ() - raw.front().getMZ())/resolution);
      if(index >= distribution_.size()){

        LOG_INFO << index <<endl;
        
      }
      auto& mass = distribution_[index].getPosition() ;

      mass = mass == 0 ?
             raw.front().getMZ() * index :
             mass;
      distribution_[index].setIntensity(distribution_[index].getIntensity() + p.getIntensity());
    }
  }
  

  void IsotopePatternSolver::dumpIDToFile(String file)
  {
    ofstream out(file.c_str());
    for(auto& sample : distribution_)
    {
      out << sample.getMZ() << " "<<sample.getIntensity() << endl;
    }
    
    out.close();
  }
  
  /* Start of the midas interface */
  MIDAs::MIDAs(EmpiricalFormula& formula, double resolution, UInt N_):
    IsotopePatternSolver(),
    min_prob(1e-16),
    formula_(formula),
    resolution_(resolution),
    N(N_)
  {
    
  }

  MIDAs::MIDAs() : IsotopePatternSolver()
  {
  }

  //MIDAs::MIDAs(const IsotopeDistribution& isotope_distribution) : 
  //  IsotopePatternSolver(isotope_distribution)
  //{
  //}


}
