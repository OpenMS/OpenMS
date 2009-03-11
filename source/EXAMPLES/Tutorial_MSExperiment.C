#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  PeakMap exp;

  for (Size i=0; i<4; ++i)
  {
    PeakSpectrum spectrum;
    spectrum.setRT(i);
    spectrum.setMSLevel(1);
    for (Real mz=500.0; mz<=900; mz+=100.0)
    { 
      Peak1D peak;
      peak.setMZ(mz+i);
      spectrum.push_back(peak);
    }
    exp.push_back(spectrum);
  } //end of creation

  for(PeakMap::AreaIterator it=exp.areaBegin(2.0, 3.0, 603.0, 802.0); it!=exp.areaEnd(); ++it)
  {
    cout << it.getRT() << " - " << it->getMZ() << endl;
  }

  for(PeakMap::Iterator s_it=exp.begin(); s_it!=exp.end(); ++s_it)
  {
    for (PeakSpectrum::Iterator p_it=s_it->begin(); p_it!=s_it->end(); ++p_it)
    {
      cout << s_it->getRT() << " - " << p_it->getMZ() << endl;
    }
  }

  return 0;
} //end of main
