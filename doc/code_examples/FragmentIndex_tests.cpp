//
// Created by trapho on 10/27/23.
//
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <iostream>

using namespace OpenMS;
using namespace std;



int main()
{
  //spectrum Exp
  MzMLFile reader;
  PeakMap map;
  reader.load("/home/trapho/test/OpenMS/doc/code_examples/data/Targeted_carbonic_anhydrase_CID12pt5V_deconv.mzML", map);
  MSSpectrum spectrum = map.getSpectrum(0);


  // spectrum Theo

  FASTAFile fasta;
  vector<FASTAFile::FASTAEntry> entries;
  fasta.load("/home/trapho/test/OpenMS/doc/code_examples/data/47128_bovine.fasta", entries);

  cout << "loaded succesfully \n" << entries[0].sequence << endl;

  FragmentIndex sdb;
  auto params = sdb.getParameters();
  //params.setValue("modifications_fixed", std::vector<std::string>{});
  //params.setValue("modifications_variable", std::vector<std::string>{"Acetyl (N-term)"});

  //sdb.setParameters(params);

  cout << "params set" << endl;
  sdb.build(entries);

  cout << "build succes \n" << sdb.getFiPeptidesSequences()[0].toString() << endl;


  auto prec = spectrum.getPrecursors();
  cout << "precussors: " << prec.size() << " mass: " << prec[0].getMZ() << endl;


  FragmentIndexScorer scorer;
  FragmentIndexScorer::InitHits inithits;
  scorer.setDB(&sdb);
  scorer.simpleScoring(spectrum, inithits);



  cout << inithits.matched_peaks_ << " " << inithits.scored_candidates_ << endl;
  for(auto i: inithits.hits_){
    FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
    cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_ <<
      " Fasta Entry: " << entries[pep.protein_idx].identifier << " Mass: "<< pep.mass << endl;
  }

  auto scorer_params = scorer.getParameters();
  scorer_params.setValue("open_search", "false");
  scorer.setParameters(scorer_params);
  inithits.clear();
  scorer.simpleScoring(spectrum, inithits);
  cout << inithits.matched_peaks_ << " " << inithits.scored_candidates_ << endl;
  for(auto i: inithits.hits_)
  {
    FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
    cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_
         << " Fasta Entry: " << entries[pep.protein_idx].identifier << endl;
  }


  ///experiment with theoretical peptide

  AASequence peptide_query = AASequence::fromString("MATVPEPINEMMAYYSDENELLFEADGPKQMKSCIQHLDLGSMGDGNIQLQISHQFYNKSFRQVVSVIVAMEKLRNSAYAHVFHDDDLRSILSFIFEEEPVIFETSSDEFLCDAPVQSIKCKLQDREQKSLVLASPCVLKALHLLSQEMNREVVFCMSFVQGEERDNKIPVALGIKDKNLYLSCVKKGDTPTLQLEEVDPKVYPKRNMEKRFVFYKTEIKNTVEFESVLYPNWYISTSQIEERPVFLGHFRGGQDITDFRMETLSP");
  AASequence peptide_special = AASequence::fromString("MATVPEPINEMMAYYSDENELLFEADGPKQMKSCIQHLDLGSMGDGNIQLQISHQFYNKSFRQVVSVIVAMEKLRNSAYAHVFHDDDLRSILSFIFEEEPVIFETSSDEFLCDAPVQSIKCKLQDREQKSLVLASPCVLKALHLLSQEMNREVVFCMSFVQGEERDNKIPVALGIKDKNLYLSCVKKGDTPTLQLEEVDPKVYPKRNMEKRFVFYKTEIKNTVEFESVLYPNWYISTSQIEERPVFLGHFRGGQDITDFRMETLSP");
  peptide_special.setCTerminalModification("Acetyl (N-term)");
  peptide_special.setNTerminalModification("Acetyl (C-term)");
  //Acetyl (N-term)
  cout << "mass query " <<  peptide_query.getMZ(1) << " " << peptide_query.getMonoWeight() << " special " << peptide_special.getMZ(1) << endl;

  sdb.addSpecialPeptide(peptide_special, entries.size());
  sdb.build(entries);

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  tsg.getSpectrum(b_y_ions, peptide_query,1, 1);
  MSSpectrum spec_theo;
  Precursor prec_theo;
  prec_theo.setMZ(peptide_query.getMZ(1));
  spec_theo.setPrecursors({prec_theo});
  spec_theo.setMSLevel(2);
  for(Peak1D p: b_y_ions){
    spec_theo.push_back(p);
  }

  inithits.clear();
  scorer.simpleScoring(spec_theo, inithits);
  cout << inithits.matched_peaks_ << " " << inithits.scored_candidates_ << endl;
  for(auto i: inithits.hits_)
  {
    FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
    cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_
         << " Fasta Entry: " << entries[pep.protein_idx].identifier << endl;
  }

  return 0;
}