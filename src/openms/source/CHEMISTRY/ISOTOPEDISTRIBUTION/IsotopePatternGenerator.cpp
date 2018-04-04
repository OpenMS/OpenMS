

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  IsotopePatternGenerator::IsotopePatternGenerator(double probability_cutoff) :
    min_prob_(probability_cutoff)
  {
  }

  IsotopePatternGenerator::IsotopePatternGenerator() :
    min_prob_(1e-15)
  {
  }
  
  IsotopePatternGenerator::~IsotopePatternGenerator()
  {
  }

}
