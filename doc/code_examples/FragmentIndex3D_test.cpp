//
// Created by trapho on 10/31/23.
//
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTD.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTDScorer.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex3D.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>

using namespace OpenMS;
using namespace std;



int main()
{
  std::vector<FASTAFile::FASTAEntry> entries {
    {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
    {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
    {"test3", "test3", "KVKLQSRPKAAPPAPGPGQLT"},

    {"test4", "test4", "MATEGMILTNHDHQIRVGV"},
  };


  //FragmentIndex3D sdb;
  //sdb.build(entries);


  //250
  FASTAFile fasta;
  vector<FASTAFile::FASTAEntry> entries2;
  fasta.load("/home/trapho/test/OpenMS/doc/code_examples/data/250_bovine.fasta", entries2);

  FragmentIndex3D sdb2;
  sdb2.build(entries2);

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  auto query = AASequence::fromString("SHHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK");
  tsg.getSpectrum(b_y_ions, query,1, 1);
  MSSpectrum spec;
  Precursor prec;
  prec.setMZ(query.getMZ(1));
  spec.setPrecursors({prec});
  spec.setMSLevel(2);
  for(Peak1D p: b_y_ions){
    spec.push_back(p);
  }
  TagGenerator tagGenerator(spec);

  tagGenerator.globalSelection();
  tagGenerator.localSelection();
  tagGenerator.generateDirectedAcyclicGraph(0.4);
  vector<MultiPeak> mPeaks;
  tagGenerator.generateAllMultiPeaks(mPeaks);


  auto range = sdb2.getPeptideRange(prec.getMZ(), {0,0});
  vector<FragmentIndexTD::Hit> hits;
  sdb2.query(hits, mPeaks[55], range, {0,0});


  for(FragmentIndexTD::Hit hit: hits){
    cout << hit.peptide_idx << endl;
    cout << sdb2.getFiPeptides().at(hit.peptide_idx).protein_idx << endl;
  }




}