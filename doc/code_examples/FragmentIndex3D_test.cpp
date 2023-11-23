//
// Created by trapho on 10/31/23.
//
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex3D.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexScorer.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <iostream>

using namespace OpenMS;
using namespace std;



int main()
{


 /* std::vector<FASTAFile::FASTAEntry> entries {
    {"test1", "test1", "LRLRACGLNFADLMARQGLY"},
    {"test2", "test2", "AAASPPLLRCLVLTGFGGYD"},
    {"test3", "test3", "KVKLQSRPKAAPPAPGPGQLT"},

    {"test4", "test4", "MATEGMILTNHDHQIRVGV"},
  };
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  auto query = AASequence::fromString("LRLRACGLNFADLMARQGLY");
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
  tagGenerator.generateDirectedAcyclicGraph(0.05);
  vector<MultiPeak> mPeaks;
  tagGenerator.generateAllMultiPeaks(mPeaks);

  FragmentIndex3D sdb;
  sdb.build(entries);

  auto range = sdb.getPeptideRange(prec.getMZ(), {0,0});
  vector<FragmentIndex::Hit> hits;
  sdb.query(hits, mPeaks[10], range, {0,0});


  for(FragmentIndex::Hit hit: hits){
    cout << hit.peptide_idx << endl;
    cout << sdb.getFiPeptides().at(hit.peptide_idx).protein_idx << endl;
  }
*/

  //250 on computational generated data
 FASTAFile fasta;
  vector<FASTAFile::FASTAEntry> entries2;
  fasta.load("/home/trapho/test/OpenMS/doc/code_examples/data/47128_bovine.fasta", entries2);
  cout << entries2[774].identifier << endl;
  auto e1 = entries2.begin();
  auto e2 = entries2.begin() + 10000;
  vector<FASTAFile::FASTAEntry> entries2_s(e1, e2);

  FragmentIndex3D sdb2;
  sdb2.build(entries2_s);
  cout << "DB build " << endl;

  // Real data
  MzMLFile reader;
  PeakMap map;
  reader.load("/home/trapho/test/OpenMS/doc/code_examples/data/Targeted_carbonic_anhydrase_CID12pt5V_deconv.mzML", map);
  MSSpectrum spectrum_exp = map.getSpectrum(0);

  //Comp data
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

  FragmentIndexScorer::InitHits hits;
  FragmentIndexScorer scorer;
  scorer.setDB(&sdb2);
  scorer.simpleScoring(spectrum_exp, hits);
  scorer.multiDimScoring(spectrum_exp, hits);
  for(auto h: hits.hits_){
    cout << h.peptide_idx_ << " " << sdb2.getFiPeptides().at(h.peptide_idx_).protein_idx << " " << h.precursor_charge_ << " " << h.num_matched_ <<endl;
  }

}