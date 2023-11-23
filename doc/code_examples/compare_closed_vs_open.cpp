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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <fstream>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{

  //values of interest:
  int tp53;
  int var;
  int next_best;
  double false_peaks;


  const ResidueModification* phosphoS = ModificationsDB::getInstance()->getModification("Phospho (S)");
  const ResidueModification* AcetylK = ModificationsDB::getInstance()->getModification("Acetyl (K)");
  AASequence peptide_query = AASequence::fromString("MEESQAELNVEPPLSQETFSDLWNLLPENNLLSSELSAPVDDLLPYTDVATWLDECPNEAPQMPEPSAPAAPPPATPAPATSWPLSSFVPSQKTYPGNYGFRLGFLQSGTAKSVTCTYSPSLNKLFCQLAKTCPVQLWVDSPPPPGTRVRAMAIYKKLEHMTEVVRRCPHHERSSDYSDGLAPPQHLIRVEGNLRAEYLDDRNTFRHSVVVPYESPEIDSECTTIHYNFMCNSSCMGGMNRRPILTIITLEDSCGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGQSCPEPPPRSTKRALPTNTSSSPQPKKKPLDGEYFTLQIRGFKRYEMFRELNDALELKDALDGREPGESRAHSSHLKSKKRPSPSCHKKPMLKREGPDSD");


  FASTAFile fasta;
  vector<FASTAFile::FASTAEntry> entries;
  fasta.load("/home/trapho/OpenMS/doc/code_examples/data/47128_bovine.fasta", entries);

  cout << "loaded succesfully \n" << entries[0].sequence << endl;
  FragmentIndex sdb;

  sdb.build(entries);

  cout << "build succes \n" << sdb.getFiPeptidesSequences()[0].toString() << endl;

  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;

  FragmentIndexScorer scorer;
  FragmentIndexScorer::InitHits inithits;
  scorer.setDB(&sdb);

  vector<pair<size_t ,const ResidueModification*> > modification_idx{make_pair(14, phosphoS),
                                                                     make_pair(384, phosphoS),
                                                                     make_pair(111, AcetylK),
                                                                     make_pair(374, AcetylK),
                                                                     make_pair(19, phosphoS),
                                                                     make_pair(313, AcetylK),
                                                                     make_pair(261, phosphoS)};

  std::ofstream experiment_one("/home/trapho/OpenMS/doc/code_examples/output/experiment_one.csv");

  experiment_one << "modi_pos, modi_mass, experiment, subexperiment, hit, variant, next_max, average_false_peaks\n";
  std::ofstream experiment_one_histo("/home/trapho/OpenMS/doc/code_examples/output/experiment_one_histo.csv");

  experiment_one_histo << "modi_pos, modi_mass, experiment, subexperiment, peptide, number_hits\n";
  for(size_t modi = 0; modi <= modification_idx.size(); modi++){
    if(modi != 0){
      peptide_query.setModification(modification_idx[modi-1].first, modification_idx[modi-1].second);
    }
    cout << "Current Peptide mass: " << peptide_query.getMZ(1) << endl;
    MSSpectrum spec_theo;
    Precursor prec_theo;
    tsg.getSpectrum(b_y_ions, peptide_query,1, 1);
    prec_theo.setMZ(peptide_query.getMZ(1));
    spec_theo.setPrecursors({prec_theo});
    spec_theo.setMSLevel(2);
    for(Peak1D peak: b_y_ions){
      spec_theo.push_back(peak);
    }

    scorer.simpleScoring(spec_theo, inithits);
    cout << "Closed search: total number of matched peaks: " << inithits.matched_peaks_ << " Number of scored candidates: " << inithits.scored_candidates_ << endl;
    bool found = false;
    bool found_var = false;
    int max = 0;
    int adjustor =3;
    for(auto i: inithits.hits_)
    {
      auto result_scorer = scorer.scoreCandidate(i, spec_theo);
      cout << result_scorer->summed_b_ << " " << result_scorer->summed_y_
           << " " <<result_scorer->matched_b_ << " " <<result_scorer->matched_y_
           << " " <<result_scorer->longest_b << " " <<result_scorer->longest_y << endl;


      experiment_one_histo << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                           << ",Standard, close," << entries[sdb.getFiPeptides()[i.peptide_idx_].protein_idx].identifier
                            << "," << i.num_matched_ << "\n";
      if(i.peptide_idx_ == 23481){
        tp53 = i.num_matched_;
        found = true;
      }
      else if(i.peptide_idx_ == 23473){
        var = i.num_matched_;
        found_var = true;
      }
      else if((i.num_matched_ > max) && i.peptide_idx_ != 23480){
        max = i.num_matched_;
      }

      FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
      cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_
           << " Fasta Entry: " << entries[pep.protein_idx].identifier << endl;
    }
    if(!found)
    {
      tp53 = 0;
      adjustor -= 2;
    }
    if(!found_var){
      var = 0;
      adjustor -= 1;
    }
    next_best = max;
    false_peaks = static_cast<double>(inithits.matched_peaks_ - tp53*2 - var) / static_cast<double>(inithits.scored_candidates_ - adjustor);
    experiment_one << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                                                                   << ",Standard, close," <<  tp53 << "," << var <<","<< next_best
                                                                           << "," << false_peaks << "\n";

    inithits.clear();

    auto sParams = scorer.getParameters();
    sParams.setValue("open_search", "true");
    scorer.setParameters(sParams);

    scorer.simpleScoring(spec_theo, inithits);
    cout << "OPENSEARCh: total number of matched peaks: " << inithits.matched_peaks_ << " Number of scored candidates: " << inithits.scored_candidates_ << endl;
    found = false;
    found_var = false;
    max = 0;
    adjustor =3;
    for(auto i: inithits.hits_)
    {
      experiment_one_histo << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                           << ",Standard, open," << entries[sdb.getFiPeptides()[i.peptide_idx_].protein_idx].identifier
                           << "," << i.num_matched_ << "\n";
      if(i.peptide_idx_ == 23481){
        tp53 = i.num_matched_;
        found = true;
      }
      else if(i.peptide_idx_ == 23473){
        var = i.num_matched_;
        found_var = true;
      }
      else if((i.num_matched_ > max) && i.peptide_idx_ != 23480){
        max = i.num_matched_;
      }

      FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
      cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_
           << " Fasta Entry: " << entries[pep.protein_idx].identifier << endl;
    }
    if(!found)
    {
      tp53 = 0;
      adjustor -= 2;
    }
    if(!found_var){
      var = 0;
      adjustor -= 1;
    }
    next_best = max;
    false_peaks = static_cast<double>(inithits.matched_peaks_ - tp53*2 - var) / static_cast<double>(inithits.scored_candidates_ - adjustor);
    experiment_one << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                   << ",Standard, open," <<  tp53 << "," << var <<","<< next_best
                   << "," << false_peaks << "\n";

    inithits.clear();
    sParams.setValue("open_search", "false");
    scorer.setParameters(sParams);

    b_y_ions.clear(true);
  }
  sdb.clear();
  peptide_query = AASequence::fromString("MEESQAELNVEPPLSQETFSDLWNLLPENNLLSSELSAPVDDLLPYTDVATWLDECPNEAPQMPEPSAPAAPPPATPAPATSWPLSSFVPSQKTYPGNYGFRLGFLQSGTAKSVTCTYSPSLNKLFCQLAKTCPVQLWVDSPPPPGTRVRAMAIYKKLEHMTEVVRRCPHHERSSDYSDGLAPPQHLIRVEGNLRAEYLDDRNTFRHSVVVPYESPEIDSECTTIHYNFMCNSSCMGGMNRRPILTIITLEDSCGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGQSCPEPPPRSTKRALPTNTSSSPQPKKKPLDGEYFTLQIRGFKRYEMFRELNDALELKDALDGREPGESRAHSSHLKSKKRPSPSCHKKPMLKREGPDSD");

  FragmentIndex3D index3D;
  index3D.build(entries);
  FragmentIndexScorer scorer2;
  scorer2.setDB(&index3D);

  for(size_t modi = 0; modi <= modification_idx.size(); modi++){
    if(modi != 0){
      peptide_query.setModification(modification_idx[modi-1].first, modification_idx[modi-1].second);
    }
    cout << "Current Peptide mass: " << peptide_query.getMZ(1) << endl;
    MSSpectrum spec_theo;
    Precursor prec_theo;
    tsg.getSpectrum(b_y_ions, peptide_query,1, 1);
    prec_theo.setMZ(peptide_query.getMZ(1));
    spec_theo.setPrecursors({prec_theo});
    spec_theo.setMSLevel(2);
    for(Peak1D peak: b_y_ions){
      spec_theo.push_back(peak);
    }

    for(double window : {100, 200, 300, 400, 500}){
      auto multiDimParams = scorer2.getParameters();
      multiDimParams.setValue("open_precursor_window", window);
      scorer2.setParameters(multiDimParams);

      scorer2.multiDimScoring(spec_theo, inithits);
      cout << window <<"-size: total number of matched peaks: " << inithits.matched_peaks_ << " Number of scored candidates: " << inithits.scored_candidates_ << endl;
      bool found = false;
      bool found_var = false;
      int max = 0;
      int adjustor =3;
      for(auto i: inithits.hits_)
      {
        experiment_one_histo << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                             << ",Multi_dim, " << window << "," << entries[sdb.getFiPeptides()[i.peptide_idx_].protein_idx].identifier
                             << "," << i.num_matched_ << "\n";
        if(i.peptide_idx_ == 23481){
          tp53 = i.num_matched_;
          found = true;
        }
        else if(i.peptide_idx_ == 23473){
          var = i.num_matched_;
          found_var = true;
        }
        else if((i.num_matched_ > max) && i.peptide_idx_ != 23480){
          max = i.num_matched_;
        }

        FragmentIndex::Peptide pep = scorer.getDb()->getFiPeptides()[i.peptide_idx_];
        cout << "#matched: " << i.num_matched_ << " isotope error:  " << i.isotope_error_ << " PepIdx: " << i.peptide_idx_ << " Charge: " << i.precursor_charge_
             << " Fasta Entry: " << entries[pep.protein_idx].identifier << endl;
      }
      if(!found)
      {
        tp53 = 0;
        adjustor -= 2;
      }
      if(!found_var){
        var = 0;
        adjustor -= 1;
      }
      next_best = max;
      false_peaks = static_cast<double>(inithits.matched_peaks_ - tp53*2 - var) / static_cast<double>(inithits.scored_candidates_ - adjustor);
      experiment_one << modification_idx[modi].first << "," << peptide_query.getMZ(1)
                     << ",Multi_dim," << window << "," <<  tp53 << "," << var <<","<< next_best
                     << "," << false_peaks << "\n";

      inithits.clear();
    }




    b_y_ions.clear(true);
  }

  experiment_one.close();
  experiment_one_histo.close();
}