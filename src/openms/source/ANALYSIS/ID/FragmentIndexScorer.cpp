// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
//#include <OpenMS/ANALYSIS/ID/FragmentIndex3D.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexScorer.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGenerator.h>
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
#include <functional>


using namespace std;

namespace OpenMS
{

  void FragmentIndexScorer::buildDB(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    db_.build(fasta_entries);
  }

  void FragmentIndexScorer::querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms)
  {
    /*
    auto* ptrBase = dynamic_cast<FragmentIndex3D*>(db_);

    if (ptrBase)
    {
      OPENMS_LOG_WARN << "The Database has the wrong format" << endl;
      return;
    }
    */
    if (!db_.isBuild())
    {
      OPENMS_LOG_WARN << "FragmentIndex not yet build \n";
      return;
    }

/* TODO:
    if(db_.getIonTypes() != )
    {
      OPENMS_LOG_WARN << "FragmentIndex was not build with ions set by the Scorer object" << endl;
      return;
    }
*/    
    if (spectrum.empty() || (spectrum.getMSLevel() != 2))
    {
      return;
    }
      
    const auto& precursor = spectrum.getPrecursors();
    if (precursor.size() != 1)
    {
      OPENMS_LOG_WARN << "Number of precursors is not equal 1 \n";
      return;
    }
    // does OpenMs work with [M] or [MH+] !!!
    // auto mz = precursor[0].getMZ - Constants::PROTON_MASS;

    // two posible modes. Precursor has a charge or we test all possible charges
    vector<size_t> charges;
    //cout << "precursor charge = " << precursor[0].getCharge() << endl;
    if (precursor[0].getCharge())
    {
      //cout << "precursor charge found" << endl;
      charges.push_back(precursor[0].getCharge());
    }
    else
    {
      for (uint16_t i = min_precursor_charge_; i <= max_precursor_charge_; i++)
      {
        charges.push_back(i);
      }
    }
    // loop over all PRECURSOR-charges

    for (uint16_t charge : charges)
    {
      SpectrumMatchesTopN candidates_charge;

      //cout << "mz" << precursor[0].getMZ() << " uw " << precursor[0].getUnchargedMass() << endl;
      float mz;
      if (precursor[0].getCharge())
        mz = precursor[0].getMZ(); // TODO: What does this do? precursor[0].getUnchargedMass()
      else
        mz = precursor[0].getMZ() * charge;

      if (open_search)
        openSearch(spectrum, mz, candidates_charge, charge);
      else
        closedSearch(spectrum, mz, candidates_charge, charge);

      trimHits(candidates_charge);
      sms += candidates_charge;
    }
    trimHits(sms);
  }

  void FragmentIndexScorer::queryPeak(SpectrumMatchesTopN& candidates, const OpenMS::Peak1D& peak, pair<size_t, size_t> candidates_range, int16_t isotope_error, uint16_t precursor_charge)
  {
    for (uint16_t fragment_charge = 1; fragment_charge <= max_fragment_charge_; fragment_charge++)
    {
      for (auto hit : db_.query(peak, candidates_range, fragment_charge))
      {
        // the following part is 1:1 from sage
        size_t idx = hit.peptide_idx - candidates_range.first;

        auto& source = candidates.hits_[idx];
        if (source.num_matched_ == 0)
        {
          ++candidates.scored_candidates_;
          source.precursor_charge_ = precursor_charge;
          source.peptide_idx_ = hit.peptide_idx;
          source.isotope_error_ = isotope_error;
        }
        ++source.num_matched_;
        ++candidates.matched_peaks_;
      }
    }
  }

  void FragmentIndexScorer::closedSearch(const MSSpectrum& spectrum, float mz, SpectrumMatchesTopN& sms, uint16_t charge)
  {
    for (int16_t isotope_error = min_isotope_error_; isotope_error <= max_isotope_error_; isotope_error++)
    {
      SpectrumMatchesTopN candidates_iso_error;
      mz += isotope_error * Constants::NEUTRON_MASS;
      auto candidates_range = db_.getPeptidesInPrecursorRange(mz, {0, 0}); // for the simple search we do not apply any modification window!!
      candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first + 1);

      for (const Peak1D& peak : spectrum)
      {
        queryPeak(candidates_iso_error, peak, candidates_range, isotope_error, charge);
      }
      // take only top 50 hits
      trimHits(candidates_iso_error);
      sms += candidates_iso_error;
    }
  }


  void FragmentIndexScorer::openSearch(const MSSpectrum& spectrum, float precursor_mass, SpectrumMatchesTopN& sms, uint16_t charge)
  {
    pair<float, float> precursor_window;
    if (open_precursor_window == 0)
      precursor_window = make_pair(-precursor_mass * 0.05, precursor_mass * 0.01);
    else
      precursor_window = make_pair(-open_precursor_window, open_precursor_window); // sage took (-500, 100) fixed
    OPENMS_LOG_WARN << "The Precursor mass window was set to: " << precursor_window.first << " " << precursor_window.second << endl;
    auto candidates_range = db_.getPeptidesInPrecursorRange(precursor_mass, precursor_window);

    // adjust the window. with each decreasing peak the prob. that we have all modifications inside gets lower
    // the linear decrease of the window below is a very coarse approximation and must be improved
    pair<float, float> fragment_window;

    SpectrumMatchesTopN hits_per_window;
    hits_per_window.hits_.resize(candidates_range.second - candidates_range.first + 1);

    for (auto peak = spectrum.end() - 1; peak >= spectrum.begin(); peak--)
    {
      queryPeak(hits_per_window, *peak, candidates_range, 0, charge);
    }
    trimHits(hits_per_window);
    sms += hits_per_window;
  }

  void FragmentIndexScorer::extractHits(OpenMS::FragmentIndexScorer::SpectrumMatchesTopN& candidates, const vector<FragmentIndex::Hit>& hits, uint32_t charge, int16_t isotope_error,
                                          std::pair<size_t, size_t> peptide_range)
  {
    for (FragmentIndex::Hit hit : hits)
    {
      // the following part is 1:1 from sage
      size_t idx = hit.peptide_idx - peptide_range.first;

      auto source = &candidates.hits_[idx];
      if (source->num_matched_ == 0)
      {
        candidates.scored_candidates_ += 1;
        source->precursor_charge_ = charge;
        source->peptide_idx_ = hit.peptide_idx;
        source->isotope_error_ = isotope_error;
      }
      source->num_matched_ += 1;
      candidates.matched_peaks_ += 1;
    }
  }

/*
  void FragmentIndexScorer::multiDimScoring(const OpenMS::MSSpectrum& spectrum, OpenMS::FragmentIndexScorer::SpectrumMatchesTopN& SpectrumMatchesTopN)
  {
    // 1.) First check all requirements for the function to actually work
    // a) We have selected the correct database type
    auto* ptrDerived = dynamic_cast<FragmentIndex3D*>(db_);
    if (!ptrDerived)
    {
      OPENMS_LOG_WARN << "The Database has the wrong format" << endl;
      return;
    }
    // b) The database was build
    if (!db_.isBuild())
    {
      OPENMS_LOG_WARN << "FragmentIndex not yet build \n";
      return;
    }
    // c) The query spectrum has the correct MS level and contains data
    if (spectrum.empty() || (spectrum.getMSLevel() != 2))
      return;

    // d) The number of precursors is correct
    auto precursor = spectrum.getPrecursors();
    if (precursor.size() != 1)
    {
      OPENMS_LOG_WARN << "Number of precursors is not equal 1 \n";
      return;
    }


    // 2.) Generate all MultiPeaks (Tags) (this should introduce all possible fragment charges)
    TagGenerator tagGenerator(spectrum);
    tagGenerator.globalSelection();
    tagGenerator.localSelection();
    cout << "DEBUG: frag mz tol: " << db_.getParameters().getValue("fragment_mz_tolerance") << endl;
    tagGenerator.generateDirectedAcyclicGraph(db_.getParameters().getValue("fragment_mz_tolerance"));
    vector<MultiPeak> mPeaks;
    tagGenerator.generateAllMultiPeaks(mPeaks, db_.getParameters().getValue("depth"));

    // 3.) Loop over all precursor charges
    vector<size_t> charges;
    cout << "precursor charge = " << precursor[0].getCharge() << endl;
    if (precursor[0].getCharge())
    {
      cout << "precursor charge found" << endl;
      charges.push_back(precursor[0].getCharge());
    }
    else
    {
      for (size_t i = min_precursor_charge_; i <= max_precursor_charge_; i++)
      {
        charges.push_back(i);
      }
    }
    // loop over all PRECURSOR-charges
    for (size_t charge : charges)
    {
      SpectrumMatchesTopN candidates_charge;
      vector<FragmentIndex::Hit> hits_charge;
      float mz = precursor[0].getMZ() * static_cast<float>(charge);
      auto range = ptrDerived->getPeptidesInPrecursorRange(mz, {-open_precursor_window, open_precursor_window});
      candidates_charge.hits_.resize(range.second - range.first + 1);
      for (const MultiPeak& mp : mPeaks)
      {
        ptrDerived->query(hits_charge, mp, range, {-open_precursor_window, open_precursor_window});
        extractHits(candidates_charge, hits_charge, charge, 0, range);
        hits_charge.clear();
      }

      trimHits(candidates_charge);
      SpectrumMatchesTopN += candidates_charge;
    }
    trimHits(SpectrumMatchesTopN);
  }
*/

  void FragmentIndexScorer::trimHits(OpenMS::FragmentIndexScorer::SpectrumMatchesTopN& init_hits)
  {
    buildKMinHeap<SpectrumMatch, uint32_t>(init_hits.hits_, max_processed_hits_, [](SpectrumMatch a) { return a.num_matched_; });
    if (max_processed_hits_ < init_hits.hits_.size())
    {
      init_hits.hits_.resize(max_processed_hits_);
    }
  }


  template<class T, class A>
  void FragmentIndexScorer::buildKMinHeap(vector<T>& slice, uint32_t size, A (*access)(T))
  {
    if (slice.size() <= size)
      return; // in such case we do not have to perform this step

    // build k-heap
    for (size_t idx = size; idx > 0; idx--)
    {
      minHeapify(slice, idx, size - 1, access);
    }

    // insert elements from outside into the heap, if they are larger than the current min element

    for (size_t outer_idx = size; outer_idx < slice.size(); outer_idx++)
    {
      if (access(slice[outer_idx]) > access(slice[0]))
      {
        iter_swap(slice.begin(), slice.begin() + outer_idx);
        minHeapify(slice, 0, size, access);
      }
    }
  }
  template<class T, class A>
  void FragmentIndexScorer::minHeapify(vector<T>& slice, size_t idx, uint32_t size, A (*access)(T))
  {
    while (true)
    { // TODO: maybe build in a failsafe or something
      size_t min = idx;
      size_t left = idx * 2 + 1;
      size_t right = idx * 2 + 2;

      if ((left < size) && (access(slice[left]) < access(slice[min])))
      {
        min = left;
      }
      if ((right < size) && (access(slice[right]) < access(slice[min])))
      {
        min = right;
      }
      if (idx != min)
      {
        iter_swap(slice.begin() + idx, slice.begin() + min);
        idx = min;
      }
      else
        return;
    }
  }

  void FragmentIndexScorer::testHeapify()
  {
    vector<int> test {10, 5, 6, 2, 20, 1, 4, 10, 1, 3, 5, 6, 1, 2, 45, 66, 89, 23, 4, 55, 7, 100, 3, 2, 44, 645, 58};

    buildKMinHeap<int, int>(test, 9, [](int a) { return a; });

    for (int dfs : test)
    {
      cout << dfs << endl;
    }
  }

  void FragmentIndexScorer::buildKMinHeapforTag(std::vector<TagGenerator::IdxAndIntensity>& slice, uint32_t size)
  {
    buildKMinHeap<TagGenerator::IdxAndIntensity, float>(slice, size, [](TagGenerator::IdxAndIntensity a) { return a.intensity_; });
  }

  FragmentIndexScorer::FragmentIndexScorer() : DefaultParamHandler("FragmentIndexScorer")
  {
    defaults_.setValue("Min_matched_peaks", 5, "Minimal number of matched ions to report a PSM");
    defaults_.setValue("min_isotope_error", 0, "Precursor isotope error");
    defaults_.setValue("max_isotope_error", 0, "precursor isotope error");
    defaults_.setValue("min_precursor_charge", 1, "min precursor charge");
    defaults_.setValue("max_precursor_charge", 1, "max precursor charge");
    defaults_.setValue("max_fragment_charge", 1, "max fragment charge");
    defaults_.setValue("max_processed_hits", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue("open_search", "false", "Open or standard search");
    defaults_.setValue("open_precursor_window", 200.0, "Size of the open precursor window, if set to 0 it is set automatically");
    defaults_.setValue("open_fragment_window", 0, "Size of the (multiple) open fragment windows");
    vector<string> frag{"a_x", "b_y", "c_z"};
    defaults_.setValue("fragmentation", "b_y");
    defaults_.setValidStrings("fragmentation", frag);


    defaultsToParam_();
  }

  void FragmentIndexScorer::updateMembers_()
  {
    min_matched_peaks_ = param_.getValue("Min_matched_peaks");
    min_isotope_error_ = param_.getValue("min_isotope_error");
    max_isotope_error_ = param_.getValue("max_isotope_error");
    min_precursor_charge_ = param_.getValue("min_precursor_charge");
    max_precursor_charge_ = param_.getValue("max_precursor_charge");
    max_fragment_charge_ = param_.getValue("max_fragment_charge");
    max_processed_hits_ = param_.getValue("max_processed_hits");
    open_search = param_.getValue("open_search").toBool();
    open_precursor_window = param_.getValue("open_precursor_window");
    open_fragment_window = param_.getValue("open_fragment_window");
    fragmentation_ = param_.getValue("fragmentation").toString();
  }
  const FragmentIndex& FragmentIndexScorer::getDB() const
  {
    return db_;
  }
} // namespace OpenMS