#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  MSSpectrum<> spectrum;
  Peak1D peak;

  for (Real mz=1500.0; mz>=500; mz-=100.0)
  { 
    peak.setMZ(mz);
    spectrum.push_back(peak);
  }

  spectrum.sortByPosition();

  MSSpectrum<>::Iterator it;
  for(it=spectrum.MZBegin(800.0); it!=spectrum.MZEnd(1000.0); ++it)
  {
    cout << it->getMZ() << endl;
  }

  return 0;
} //end of main
