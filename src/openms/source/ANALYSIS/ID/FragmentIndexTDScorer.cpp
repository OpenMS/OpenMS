//
// Created by trapho on 10/12/23.
//
#include <OpenMS/ANALYSIS/ID/FragmentIndexTDScorer.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTD.h>

#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>

#include <functional>


using namespace std;

namespace OpenMS
{

  void FragmentIndexTDScorer::setDB(FragmentIndexTD& db){
    db_ = db;
  }
  void FragmentIndexTDScorer::buildDB(const std::vector<FASTAFile::FASTAEntry> & fasta_entries){
    db_.build(fasta_entries);
  }



  void FragmentIndexTDScorer::simpleScoring(MSSpectrum& specturm, InitHits& initHits){

    if (!db_.isBuild()){
      OPENMS_LOG_WARN << "FragmentIndex not yet build \n";
      return;
    }

    if(specturm.empty() || (specturm.getMSLevel() != 2)) return;

    auto precursor = specturm.getPrecursors();
    if(precursor.size() != 1){
      OPENMS_LOG_WARN << "Number of precursors is not equal 1 \n";
      return;
    }
    //does OpenMs work with [M] or [MH+] !!!
    //auto mz = precursor[0].getMZ - Constants::PROTON_MASS;

    //two posible modes. Precursor has a charge or we test all possible charges
    vector<size_t> charges;
    cout << "precursor charge = " << precursor[0].getCharge() << endl;
    if(precursor[0].getCharge()){
      cout << "precursor charge found" << endl;
      charges.push_back(precursor[0].getCharge());
    }else{
      for(size_t i = min_precursor_charge_; i <= max_precursor_charge_; i++){
        charges.push_back(i);
      }
    }
    //loop over all PRECURSOR-charges


    for(size_t charge: charges){
      InitHits candidates_charge;

      double mz;
      if(precursor[0].getCharge())
        mz = precursor[0].getMZ();
      else
        mz = precursor[0].getMZ() * charge;

      if(open_search)
        openScoring(specturm, mz, candidates_charge, charge);
      else
        closedScoring(specturm, mz, candidates_charge, charge);

      trimHits(candidates_charge);
      initHits += candidates_charge;
    }
    trimHits(initHits);
  }

  void FragmentIndexTDScorer::closedScoring(OpenMS::MSSpectrum& spectrum, double mz, OpenMS::FragmentIndexTDScorer::InitHits& initHits, uint32_t charge)
  {
    for (int16_t isotope_error = min_isotope_error_; isotope_error <= max_isotope_error_; isotope_error++){
      InitHits candidates_iso_error;
      mz += isotope_error * Constants::NEUTRON_MASS;
      auto candidates_range = db_.getPeptideRange(mz, {0,0});  // for the simple search we do not apply any modification window!!
      candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first + 1);

      for(Peak1D peak: spectrum){

        //TODO: should we loop over the charges of the peaks?? or adjust the mass in general?
        for(auto hit: db_.query(peak, candidates_range, {0,0})){

          // the following part is 1:1 from sage
          size_t idx = hit.peptide_idx - candidates_range.first;

          auto source = &candidates_iso_error.hits_[idx];
          if (source->num_matched_ == 0){
            candidates_iso_error.scored_candidates_ +=1;
            source->precursor_charge_ = charge;
            source->peptide_idx_ = hit.peptide_idx;
            source->isotope_error_ = isotope_error;
          }
          source->num_matched_+= 1;
          candidates_iso_error.matched_peaks_+=1;
        }
      }
      //take only top 50 hits
      trimHits(candidates_iso_error);
      initHits += candidates_iso_error;
    }
  }


  void FragmentIndexTDScorer::openScoring(MSSpectrum& spectrum, double mz, InitHits& open_init_hits, uint32_t charge)
  {
    pair<double, double> precursor_window;
    if(open_precursor_window == 0)
      precursor_window = make_pair(-mz * 0.05, mz * 0.01);
    else
      precursor_window = make_pair(-open_precursor_window,  open_precursor_window);//sage took (-500, 100) fixed
    OPENMS_LOG_WARN << "The Precursor mass window was set to: " << precursor_window.first << " " << precursor_window.second<< endl;
    auto candidates_range = db_.getPeptideRange(mz, precursor_window);

    // adjust the window. with each decreasing peak the prob. that we have all modifications inside gets lower
    // the linear decrease of the window below is a very coarse approximation and must be improved
    pair<double, double> fragment_window;
    //TODO: update the adaption of the window
    for (double window_step = precursor_window.first; window_step <= precursor_window.second; window_step += open_fragment_window){
      InitHits hits_per_window;
      hits_per_window.hits_.resize(candidates_range.second - candidates_range.first + 1);
      fragment_window = make_pair(window_step, window_step + open_fragment_window > precursor_window.second ? precursor_window.second : window_step + open_fragment_window);

      for(auto peak = spectrum.end()-1; peak >= spectrum.begin(); peak--){
        //TODO: should we loop over the charges of the peaks?? or adjust the mass in general?
        for(auto hit: db_.query(*peak, candidates_range, fragment_window)){

          // the following part is 1:1 from sage
          size_t idx = hit.peptide_idx - candidates_range.first;
          auto source = &hits_per_window.hits_[idx];
          if (source->num_matched_ == 0){
            hits_per_window.scored_candidates_ +=1;
            source->precursor_charge_ = charge;
            source->peptide_idx_ = hit.peptide_idx;
            source->isotope_error_ = window_step;
          }
          source->num_matched_+= 1;
          hits_per_window.matched_peaks_+=1;

        }
      }
      trimHits(hits_per_window);
      open_init_hits += hits_per_window;
    }
  }

  void FragmentIndexTDScorer::trimHits(OpenMS::FragmentIndexTDScorer::InitHits& init_hits)
  {
    buildKMinHeap<PreHits, uint32_t>(init_hits.hits_,  max_processed_hits_, [](PreHits a){return a.num_matched_;});
    if(max_processed_hits_ < init_hits.hits_.size()){
      init_hits.hits_.resize(max_processed_hits_);
    }
  }



  template<class T, class A>
  void FragmentIndexTDScorer::buildKMinHeap(vector<T>& slice, uint32_t size, A (*access)(T))
  {
    if (slice.size() <= size) return;  // in such case we do not have to perform this step

    // build k-heap
    for (size_t idx = size ; idx > 0; idx--){
      minHeapify(slice, idx, size-1, access);
    }

    //insert elements from outside into the heap, if they are larger than the current min element

    for (size_t outer_idx = size; outer_idx < slice.size(); outer_idx++){
      if(access(slice[outer_idx]) > access(slice[0])){
        iter_swap(slice.begin(), slice.begin()+outer_idx);
        minHeapify(slice, 0, size, access);
      }
    }


  }
  template<class T, class A>
  void FragmentIndexTDScorer::minHeapify(vector<T>& slice, size_t idx, uint32_t size, A (*access) (T))
  {

    while(true){  // TODO: maybe build in a failsafe or something
      size_t min = idx;
      size_t left  = idx * 2 +1;
      size_t right = idx * 2 +2;

      if((left < size) && (access(slice[left]) < access(slice[min]))){
        min = left;
      }
      if((right < size) && (access(slice[right]) < access(slice[min]))){
        min = right;
      }
      if(idx != min){
        iter_swap(slice.begin()+idx, slice.begin() + min);
        idx = min;
      }
      else
        return;
    }
  }

  void FragmentIndexTDScorer::testHeapify()
  {
    vector<int> test{10,5,6,2,20,1,4,10,1, 3, 5, 6 ,1, 2, 45, 66, 89, 23, 4, 55, 7, 100, 3, 2, 44, 645, 58};

    buildKMinHeap<int, int>(test, 9, [](int a){return a;});

    for(int dfs: test){
      cout << dfs << endl;
    }
  }


  void FragmentIndexTDScorer::buildKMinHeapforTag(std::vector<TagGenerator::IdxAndIntensity>& slice, uint32_t size)
  {
      buildKMinHeap<TagGenerator::IdxAndIntensity, double>(slice, size, [](TagGenerator::IdxAndIntensity a){return a.intensity_;});
  }

  FragmentIndexTDScorer::FragmentIndexTDScorer() : DefaultParamHandler("FragmentIndexTDScorer")
  {
    defaults_.setValue("Min_matched_peaks", 5, "Minimal number of matched ions to report a PSM");
    defaults_.setValue("min_isotope_error", -1, "Precursor isotope error");
    defaults_.setValue("max_isotope_error", 1, "precursor isotope error");
    defaults_.setValue("min_precursor_charge", 1, "min precursor charge");
    defaults_.setValue("max_precursor_charge", 1, "max precursor charge");
    defaults_.setValue("max_fragment_charge", 1, "max fragment charge");
    defaults_.setValue("max_processed_hits", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue("open_search", "true", "Open or standard search");
    defaults_.setValue("open_precursor_window", 200.0, "Size of the open precursor window, if set to 0 it is set automatically");
    defaults_.setValue("open_fragment_window", 1.0, "Size of the (multiple) open fragment windows");

    defaultsToParam_();
  }

  void FragmentIndexTDScorer::updateMembers_()
  {
    min_matched_peaks_ = param_.getValue("Min_matched_peaks");
    min_isotope_error_ = param_.getValue("min_isotope_error");
    max_isotope_error_ = param_.getValue("max_isotope_error");
    min_precursor_charge_ = param_.getValue("min_precursor_charge");
    max_precursor_charge_ = param_.getValue("max_precursor_charge");
    max_processed_hits_ = param_.getValue("max_processed_hits");
    open_search = param_.getValue("open_search").toBool();
    open_precursor_window = param_.getValue("open_precursor_window");
    open_fragment_window = param_.getValue("open_fragment_window");
  }
  const FragmentIndexTD& FragmentIndexTDScorer::getDb() const
  {
    return db_;
  }
}