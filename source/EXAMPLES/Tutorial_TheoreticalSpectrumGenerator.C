#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  TheoreticalSpectrumGenerator tsg;
  RichPeakSpectrum spec1, spec2;
  AASequence peptide("DFPIANGER");

  tsg.addPeaks(spec1, peptide, Residue::YIon, 1);

  tsg.getSpectrum(spec2, peptide, 2);

  cout << "Spectrum 1 has " << spec1.size() << " peaks. " << endl;
  cout << "Spectrum 2 has " << spec2.size() << " peaks. " << endl;

  return 0;
} //end of main
