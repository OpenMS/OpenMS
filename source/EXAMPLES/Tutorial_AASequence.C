#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  AASequence seq("DFPIANGER");

  AASequence prefix(seq.getPrefix(4));
  AASequence suffix(seq.getSuffix(5));

  cout   << seq << " "
        << prefix << " "
        << suffix << " " 
        << seq.getAverageWeight() << endl;
          
  return 0;
} //end of main
