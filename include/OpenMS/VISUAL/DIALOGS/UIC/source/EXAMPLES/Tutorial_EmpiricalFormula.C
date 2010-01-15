#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  EmpiricalFormula methanol("CH3OH"), water("H2O");

  EmpiricalFormula sum = methanol + water;

  cout   << sum << " " 
        << sum.getNumberOf("Carbon") << " " 
        << sum.getAverageWeight() << endl;
  
  IsotopeDistribution iso_dist = sum.getIsotopeDistribution(3);

  for (IsotopeDistribution::ConstIterator it = iso_dist.begin(); it != iso_dist.end(); ++it)
  {
    cout << it->first << " " << it->second << endl;
  }

  return 0;
} //end of main
