// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndex3D.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGenerator.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#include <functional>


using namespace std;


namespace OpenMS
{
  FragmentIndex3D::FragmentIndex3D() : DefaultParamHandler("FragmentIndex3D")
  {
    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", {"true", "false"});

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", {"true", "false"});

    defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", {"true", "false"});

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", {"true", "false"});

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", {"true", "false"});

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", {"true", "false"});

    defaults_.setValue("depth", 3, "The number of adjacent peaks that are taken into account");

    vector<string> tolerance_units {"DA", "PPM"}; // TODO: check if same string in other OpenMS files
    defaults_.setValue("fragment_min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment_max_mz", 500000, "Maximal fragment mz for database");


    defaults_.setValue("fragment_min_mz", 0, "Minimal fragment mz for database");
    defaults_.setValue("fragment_max_mz", 5000000, "Maximal fragment mz for database");
    defaults_.setValue("precursor_mz_tolerance", 2.0, "Tolerance for precursor-m/z in search");
    defaults_.setValue("fragment_mz_tolerance", 0.05, "Tolerance for fragment-m/z in search");
    defaults_.setValue("precursor_mz_tolerance_unit", "DA", "Unit of tolerance for precursor-m/z");
    defaults_.setValidStrings("precursor_mz_tolerance_unit", tolerance_units);
    defaults_.setValue("fragment_mz_tolerance_unit", "DA", "Unit of tolerance for fragment-m/z");
    defaults_.setValidStrings("fragment_mz_tolerance_unit", tolerance_units);
    is_build_ = false;

    defaults_.setValue("min_matched_peaks", 5, "Minimal number of matched ions to report a PSM");
    defaults_.setValue("min_precursor_charge", 1, "min precursor charge");
    defaults_.setValue("max_precursor_charge", 4, "max precursor charge");
    defaults_.setValue("max_fragment_charge", 1, "max fragment charge");
    defaults_.setValue("max_processed_hits", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue("open_precursor_window_lower", -500.0, "lower bound of the open precursor window");
    defaults_.setValue("open_precursor_window_upper", 500.0, "upper bound of the open precursor window");
    defaults_.setValue("open_fragment_window_lower", -200.0, "lower bound of the open fragment window");
    defaults_.setValue("open_fragment_window_upper", 500.0, "upper bound of the open fragment window");

    defaultsToParam_();
  }
  void FragmentIndex3D::updateMembers_()
  {
    depth_ = param_.getValue("depth");
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    fragment_min_mz_ = param_.getValue("fragment_min_mz");
    fragment_max_mz_ = param_.getValue("fragment_max_mz");
    precursor_mz_tolerance_ = param_.getValue("precursor_mz_tolerance");
    fragment_mz_tolerance_ = param_.getValue("fragment_mz_tolerance");
    precursor_mz_tolerance_unit_ppm_ = param_.getValue("precursor_mz_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment_mz_tolerance_unit").toString() == "ppm";
    min_matched_peaks_ = param_.getValue("min_matched_peaks");
    min_precursor_charge_ = param_.getValue("min_precursor_charge");
    max_precursor_charge_ = param_.getValue("max_precursor_charge");
    max_fragment_charge_ = param_.getValue("max_fragment_charge");
    max_processed_hits_ = param_.getValue("max_processed_hits");
    open_precursor_window_lower_ = param_.getValue("open_precursor_window_lower");
    open_precursor_window_upper_ = param_.getValue("open_precursor_window_upper");
    open_fragment_window_lower_ = param_.getValue("open_fragment_window_lower");
    open_fragment_window_upper_ = param_.getValue("open_fragment_window_upper");
  }


  void FragmentIndex3D::generate_peptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    size_t skipped_peptides = 0;


#pragma omp parallel
    for (SignedSize i = 0; i < fasta_entries.size(); ++i)
    {
      const FASTAFile::FASTAEntry& protein = fasta_entries[i];
      /// DIGEST (if bottom-up)

      // remove peptides containing unknown AA
      if (protein.sequence.find('X') != string::npos)
      {
        skipped_peptides++;
        continue;
      }

      /// MODIFY (if modifications are specified)
      AASequence unmod_peptide = AASequence::fromString(protein.sequence);
      float unmodified_mz = unmod_peptide.getMZ(1);


#pragma omp critical(FIIndex)
      {
        fi_peptides_.push_back({static_cast<UInt32>(i), unmodified_mz});
      }
    }
    if (skipped_peptides > 0)
    {
      OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unkown AA \n";
    }
    // sort the peptide vector, critical for following steps
    sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b) { return std::tie(a.precursor_mz, a.protein_idx_) < std::tie(b.precursor_mz, b.protein_idx_); });
  }

  void FragmentIndex3D::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    TheoreticalSpectrumGenerator tsg;     /// here we actually need the TSG because we need the meta info
    // store ion types for each peak
    Param tsg_settings = tsg.getParameters();
    tsg_settings.setValue("add_metainfo", "true");
    tsg.setParameters(tsg_settings);

    PeakSpectrum b_y_ions;

    /// generate all Peptides
    generate_peptides(fasta_entries);

    size_t peptide_idx = 0;

    /// For each Peptides get all theoretical b and y ions // TODO: include other fragmentation methods
    for (const Peptide& pep : fi_peptides_)
    {
      AASequence pep_sequence = AASequence::fromString(fasta_entries[pep.protein_idx_].sequence);

      tsg.getSpectrum(b_y_ions, pep_sequence, 1, 1);
      TagGenerator tag_generator {b_y_ions};

      tag_generator.generateAllMultiFragments(fi_fragments_, depth_, peptide_idx, fragment_min_mz_, fragment_max_mz_);

      peptide_idx++;
      b_y_ions.clear(true);
    }

    /// Calculate bucket size
    bucketsize_ = static_cast<size_t>(pow(static_cast<float>(fi_fragments_.size()), (1.0 / static_cast<float>(depth_ + 2)))); // The bucketsize depends on the dimensionallity, thus on the depth
    vector<size_t> iterBucketsize;
    iterBucketsize.push_back(fi_fragments_.size());
    cout << "DB size: " << fi_fragments_.size() << "\n Bucket sizes: ";
    for (uint16_t d = depth_ + 1; d >= 1; d--)
    {
      iterBucketsize.push_back(pow(bucketsize_, d));
      cout << iterBucketsize[iterBucketsize.size() - 1] << endl;
    }


    /// 1.) First we sort According to follow up peaks. Keep in mind our data set is (x+2)-D, where x is the depth.
    for (uint16_t d = 0; d < depth_; d++)
    {
      if (d > 0)
        follow_up_peaks_buckets_min_mz.emplace_back();
      for (size_t i = 0; i < fi_fragments_.size(); i += iterBucketsize[d])
      {
        if (d > 0)
          follow_up_peaks_buckets_min_mz[follow_up_peaks_buckets_min_mz.size() - 1].emplace_back(fi_fragments_[i].getFollowUpPeaks()[follow_up_peaks_buckets_min_mz.size() - 1]);
        auto sec_bucket_start = fi_fragments_.begin() + i;
        auto sec_bucket_end = (i + (iterBucketsize[d])) > fi_fragments_.size() ? fi_fragments_.end() : (sec_bucket_start + iterBucketsize[d]);
        sort(sec_bucket_start, sec_bucket_end, [d](const MultiFragment& a, const MultiFragment& b) { return a.getFollowUpPeaks()[d] < b.getFollowUpPeaks()[d]; });
      }
    }


    /// 2.) Sort by Fragment mass
    follow_up_peaks_buckets_min_mz.emplace_back();
    for (size_t i = 0; i < fi_fragments_.size(); i += (iterBucketsize[depth_]))
    {
      follow_up_peaks_buckets_min_mz[follow_up_peaks_buckets_min_mz.size() - 1].emplace_back(fi_fragments_[i].getFollowUpPeaks()[follow_up_peaks_buckets_min_mz.size() - 1]);
      auto sec_bucket_start = fi_fragments_.begin() + i;
      auto sec_bucket_end = (i + (iterBucketsize[depth_])) > fi_fragments_.size() ? fi_fragments_.end() : sec_bucket_start + iterBucketsize[depth_];

      sort(sec_bucket_start, sec_bucket_end, [](const MultiFragment& a, const MultiFragment& b) { return a.getFragmentMz() < b.getFragmentMz(); });
    }


    /// 3.) Sort by Precursor mass
    for (size_t j = 0; j < fi_fragments_.size(); j += iterBucketsize[depth_ + 1])
    {
      bucket_min_mz_.emplace_back(fi_fragments_[j].getFragmentMz());

      auto bucket_start = fi_fragments_.begin() + j;
      auto bucket_end = (j + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

      sort(bucket_start, bucket_end, [](const MultiFragment& a, const MultiFragment& b) { return a.getPeptideIdx() < b.getPeptideIdx(); });
    }

    is_build_ = true;
  }


  bool FragmentIndex3D::inRange(float hit, float query, float tolerance, std::pair<float, float> window)
  {
    return (query >= (hit - tolerance + window.first)) && (query <= (hit + tolerance + window.second));
  }

  bool FragmentIndex3D::inRangeFollowUpPeaks(std::vector<float> hit, std::vector<float> query, float tolerance)
  {
    if (hit.size() != query.size())
    {
      OPENMS_LOG_WARN << "The Query has a different Neighborhood depth than the Database!";
      return false;
    }
    for (size_t i = 0; i < hit.size(); i++)
    {
      if ((hit[i] <= (query[i] - tolerance)) || (hit[i] >= (query[i] + tolerance)))
        return false;
    }
    return true;
  }


  std::pair<size_t, size_t> FragmentIndex3D::getPeptidesInPrecursorRange(float precursor_mass, std::pair<float, float> window)
  {
    float prec_tol = precursor_mz_tolerance_unit_ppm_ ? Math::ppmToMass(precursor_mz_tolerance_, precursor_mass) : precursor_mz_tolerance_;

    auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass - prec_tol + window.first, [](const Peptide& a, float b) { return a.precursor_mz < b; });
    auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass + prec_tol + window.second, [](float b, const Peptide& a) { return b < a.precursor_mz; });
    return make_pair(std::distance(fi_peptides_.begin(), left_it), std::distance(fi_peptides_.begin(), right_it));
  }

  void FragmentIndex3D::query(vector<FragmentIndex3D::Hit>& hits, const MultiPeak& peak, std::pair<size_t, size_t> peptide_idx_range, std::pair<float, float> window)
  {
    float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, (float)peak.getPeak().getMZ()) : fragment_mz_tolerance_; // TODO??? apply ppm to mass * charge or not

    recursiveQuery(hits, peak, peptide_idx_range, window, depth_ + 1, 0, frag_tol);
  }

  void FragmentIndex3D::recursiveQuery(vector<OpenMS::FragmentIndex3D::Hit>& hits, const MultiPeak& peak, std::pair<size_t, size_t> peptide_idx_range, std::pair<float, float> window,
                                       size_t recursion_step,
                                       size_t current_slice, // From the last recursiv step. Holds the info in which branch of the tree we are in
                                       float fragment_tolerance)
  {
    vector<float>* current_level;
    float current_query;
    std::pair<float, float> applied_window;

    if (recursion_step == 0)
    { // last (precursor mz) level of the tree. Push hits into the hits vector
      auto last_slice_begin = fi_fragments_.begin() + current_slice * bucketsize_;
      auto last_slice_end = ((current_slice + 1) * bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : fi_fragments_.begin() + (current_slice + 1) * bucketsize_ + 1;

      auto left_it = std::lower_bound(last_slice_begin, last_slice_end, peptide_idx_range.first, [](MultiFragment a, UInt32 b) { return a.getPeptideIdx() < b; });

      while (left_it != last_slice_end && left_it->getPeptideIdx() <= peptide_idx_range.second)
      {
        if (inRange(left_it->getFragmentMz(), peak.getPeak().getMZ(), fragment_tolerance, window) && inRangeFollowUpPeaks(left_it->getFollowUpPeaks(), peak.getFollowUpPeaks(), fragment_tolerance))
        {
          hits.push_back({left_it->getPeptideIdx(), left_it->getFragmentMz()});
        }
        ++left_it;
      }
      return;
    }
    if (recursion_step == 1)
    { // Fragment mz level
      current_level = &bucket_min_mz_;
      current_query = peak.getPeak().getMZ();
      applied_window = window;
    }
    if (recursion_step > 1)
    {                                                                                    // All follow up peak levels
      current_level = &(follow_up_peaks_buckets_min_mz.at(depth_ + 1 - recursion_step)); // get the current vector of interests
      current_query = peak.getFollowUpPeaks().at(depth_ + 1 - recursion_step);           // current query value of interst
      applied_window = make_pair<float, float>(0, 0);                                    // For the follow up peaks we do not have any window
    }
    auto slice_begin = (*current_level).begin() + current_slice * bucketsize_;
    auto slice_end = ((current_slice + 1) * bucketsize_) > (*current_level).size() ? (*current_level).end() : (*current_level).begin() + (current_slice + 1) * bucketsize_ + 1;
    if (recursion_step == depth_ + 1) // edge case for the very first tree layer
      slice_end = (*current_level).end();

    auto left_it = std::lower_bound(slice_begin, slice_end, current_query - fragment_tolerance + applied_window.first);
    auto right_it = std::upper_bound(slice_begin, slice_end, current_query + fragment_tolerance + applied_window.second);

    auto next_slices = make_pair(std::distance(current_level->begin(), left_it), distance(current_level->begin(), right_it));
    for (size_t next_slice = next_slices.first; next_slice < next_slices.second; next_slice++)
    {
      recursiveQuery(hits, peak, peptide_idx_range, window, recursion_step - 1, next_slice, fragment_tolerance);
    }
  }


  void FragmentIndex3D::multiDimScoring(const OpenMS::MSSpectrum& spectrum, OpenMS::FragmentIndex3D::SpectrumMatchesTopN& SpectrumMatchesTopN)
  {
    // b) The database was build
    if (!is_build_)
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
    cout << "DEBUG: frag mz tol: " << getParameters().getValue("fragment_mz_tolerance") << endl;
    tagGenerator.generateDirectedAcyclicGraph(getParameters().getValue("fragment_mz_tolerance"));
    vector<MultiPeak> mPeaks;
    tagGenerator.generateAllMultiPeaks(mPeaks, getParameters().getValue("depth"));

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
      FragmentIndex3D::SpectrumMatchesTopN candidates_charge;
      vector<FragmentIndex3D::Hit> hits_charge;
      float mz = precursor[0].getMZ() * static_cast<float>(charge);
      auto range = getPeptidesInPrecursorRange(mz, {-open_precursor_window_lower_, open_precursor_window_upper_});
      candidates_charge.hits_.resize(range.second - range.first + 1);
      for (const MultiPeak& mp : mPeaks)
      {
        query(hits_charge, mp, range, {-open_precursor_window_lower_, open_precursor_window_upper_});
        {
          for (FragmentIndex3D::Hit hit : hits_charge)
          {
            // the following part is 1:1 from sage
            size_t idx = hit.peptide_idx - range.first;

            auto source = &candidates_charge.hits_[idx];
            if (source->num_matched_ == 0)
            {
              candidates_charge.scored_candidates_ += 1;
              source->precursor_charge_ = charge;
              source->peptide_idx_ = hit.peptide_idx;
              source->delta_precursor_mass = mz - fi_peptides_[hit.peptide_idx].precursor_mz;
            }
            source->num_matched_ += 1;
            candidates_charge.matched_peaks_ += 1;
          }
        }
        hits_charge.clear();
      }

      trimHits(candidates_charge);
      SpectrumMatchesTopN += candidates_charge;
    }
    trimHits(SpectrumMatchesTopN);
  }


  void FragmentIndex3D::trimHits(OpenMS::FragmentIndex3D::SpectrumMatchesTopN& init_hits) const
  {
    if (init_hits.hits_.size() > max_processed_hits_)
    {
      std::partial_sort(init_hits.hits_.begin(), init_hits.hits_.begin() + max_processed_hits_, init_hits.hits_.end(), [](const SpectrumMatch& a, const SpectrumMatch& b) {
        if (a.num_matched_ != b.num_matched_)
        {
          return a.num_matched_ > b.num_matched_;
        }
        else
        {
          return std::tie(a.delta_precursor_mass, a.precursor_charge_) < std::tie(b.delta_precursor_mass, b.precursor_charge_);
        }
      });

      init_hits.hits_.resize(max_processed_hits_);
    }
    for (auto hit_iter = init_hits.hits_.rbegin(); hit_iter != init_hits.hits_.rend(); hit_iter++)
    {
      if (hit_iter->num_matched_ >= min_matched_peaks_) // search for the first element that should be included
      {
        init_hits.hits_.resize(init_hits.hits_.size() - (distance(init_hits.hits_.rbegin(), hit_iter)));
        break;
      }
      if (hit_iter == init_hits.hits_.rend() - 1) // we reached the last element without activating the previous if statement -> no hits at all
      {
        init_hits.hits_.resize(0);
      }
    }
  }



} // namespace OpenMS