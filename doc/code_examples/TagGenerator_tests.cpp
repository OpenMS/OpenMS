//
// Created by trapho on 10/27/23.
//

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  //MzMLFile reader;
  //PeakMap map;
  //reader.load("/home/trapho/test/OpenMS/doc/code_examples/data/Targeted_carbonic_anhydrase_CID12pt5V_deconv.mzML", map);
  //MSSpectrum spectrum = map.getSpectrum(0);

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  auto queryP = AASequence::fromString("LQSRPAAPPAPGPGQLTHALLIWASGEHT");
  tsg.getSpectrum(b_y_ions, queryP,1, 1);
  MSSpectrum spec;
  Precursor prec;
  prec.setMZ(queryP.getMZ(1));
  spec.setPrecursors({prec});
  spec.setMSLevel(2);
  for(Peak1D p: b_y_ions){
    spec.push_back(p);
  }

  TagGenerator tg(spec);
  cout << "starting global selection" << endl;
  tg.globalSelection();
  cout << "global done" << endl;
  tg.localSelection();
  cout << "local done" << endl;

  tg.generateAllNodes(3);


  tg.generateDirectedAcyclicGraph(0.2);


  vector<MultiPeak> allquadpeaks;
  tg.generateAllMultiPeaks(allquadpeaks, 2);

  return 0;
}