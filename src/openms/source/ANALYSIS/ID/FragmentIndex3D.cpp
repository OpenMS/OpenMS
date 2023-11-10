//
// Created by trapho on 10/31/23.
//
#include <OpenMS/ANALYSIS/ID/FragmentIndexTD.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex3D.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>

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
  FragmentIndex3D::FragmentIndex3D() : FragmentIndexTD()
  {
    defaults_.setValue("depth", 3, "The number of adjacent peaks that are taken into account");
    defaultsToParam_();
  }
  void FragmentIndex3D::updateMembers_()
  {
    depth_ = param_.getValue("depth");
    FragmentIndexTD::updateMembers_();
  }

  void FragmentIndex3D::build (const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
    TheoreticalSpectrumGenerator tsg;
    //store ion types for each peak
    Param tsg_settings = tsg.getParameters();
    tsg_settings.setValue("add_metainfo", "true");
    tsg.setParameters(tsg_settings);

    PeakSpectrum b_y_ions;
    std::vector<Fragment> all_frags;
    /// generate all Peptides
    generate_peptides(fasta_entries);

    size_t peptide_idx = 0;
    /// For each Peptides get all theoretical b and y ions // TODO: include other fragmentation methods
    for (Peptide pep: fi_peptides_)
    {
      tsg.getSpectrum(b_y_ions, pep.sequence, 1, 1);
      TagGenerator tag_generator{b_y_ions};

      tag_generator.generateAllMultiFragments(fi_fragments_, depth_, peptide_idx, fragment_min_mz_, fragment_max_mz_);

      peptide_idx++;
      b_y_ions.clear(true);
    }

    ///Calculate bucket size
    bucketsize_ = static_cast<size_t>(pow(static_cast<double>(fi_fragments_.size()), (1.0/static_cast<double>(depth_+2))));  // The bucketsize depends on the dimensionallity, thus on the depth
    vector<size_t> iterBucketsize;
    iterBucketsize.push_back(fi_fragments_.size());
    cout << "DB size: " << fi_fragments_.size() << "\n Bucket sizes: " ;
    for(uint16_t d = depth_+1; d >= 1; d--){
      iterBucketsize.push_back(pow(bucketsize_, d));
      cout << iterBucketsize[iterBucketsize.size()-1] << endl;
    }


    ///1.) First we sort According to follow up peaks. Keep in mind our data set is (x+2)-D, where x is the depth.
    for(uint16_t d = 0; d< depth_; d++){
      if(d > 0)
        follow_up_peaks_buckets_min_mz.emplace_back();
      for(size_t i = 0; i < fi_fragments_.size(); i+= iterBucketsize[d]){
        if(d > 0)
          follow_up_peaks_buckets_min_mz[follow_up_peaks_buckets_min_mz.size()-1].emplace_back(fi_fragments_[i].getFollowUpPeaks()[follow_up_peaks_buckets_min_mz.size()-1]);
        auto sec_bucket_start = fi_fragments_.begin()+i;
        auto sec_bucket_end = (i+ (iterBucketsize[d])) > fi_fragments_.size() ? fi_fragments_.end() : (sec_bucket_start + iterBucketsize[d]);
        sort(sec_bucket_start, sec_bucket_end, [d](const MultiFragment& a, const MultiFragment& b){
          return a.getFollowUpPeaks()[d] < b.getFollowUpPeaks()[d];  });

      }

    }


    ///2.) Sort by Fragment mass
    follow_up_peaks_buckets_min_mz.emplace_back();
    for(size_t i = 0; i < fi_fragments_.size(); i += (iterBucketsize[depth_])){
      follow_up_peaks_buckets_min_mz[follow_up_peaks_buckets_min_mz.size()-1].emplace_back(fi_fragments_[i].getFollowUpPeaks()[follow_up_peaks_buckets_min_mz.size()-1]);
      auto sec_bucket_start = fi_fragments_.begin()+i;
      auto sec_bucket_end = (i+ (iterBucketsize[depth_])) > fi_fragments_.size() ? fi_fragments_.end() : sec_bucket_start + iterBucketsize[depth_];

      sort(sec_bucket_start, sec_bucket_end, [](const MultiFragment& a, const MultiFragment& b){
        return a.getFragmentMz() < b.getFragmentMz();
      });
    }


    ///3.) Sort by Precursor mass
    for (size_t j = 0; j < fi_fragments_.size(); j += iterBucketsize[depth_+1]){
      bucket_min_mz_.emplace_back(fi_fragments_[j].getFragmentMz());

      auto bucket_start = fi_fragments_.begin()+j;
      auto bucket_end = (j + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

      sort(bucket_start, bucket_end, [](const MultiFragment& a, const MultiFragment& b) {
        return a.getPeptideIdx() < b.getPeptideIdx();
      });
    }

    is_build_ = true;


  }

  bool FragmentIndex3D::inRange(double hit, double query, double tolerance, std::pair<double, double> window)
  {
    return (query >= (hit -tolerance +window.first)) && (query <= (hit + tolerance + window.second));
  }

  bool FragmentIndex3D::inRangeFollowUpPeaks(std::vector<double> hit, std::vector<double> query, double tolerance)
  {
    if(hit.size() != query.size()){
      OPENMS_LOG_WARN << "The Query has a different Neighborhood depth than the Database!";
      return false;
    }
    for(size_t i = 0; i < hit.size(); i++){
      if((hit[i] <= (query[i] - tolerance)) || (hit[i] >= (query[i] + tolerance)))
        return false;
    }
    return true;
  }

  void FragmentIndex3D::query(vector<FragmentIndexTD::Hit>& hits,
                              const MultiPeak& peak, std::pair<size_t, size_t> peptide_idx_range,
                              std::pair<double, double> window)
  {

    double frag_tol = (fragment_mz_tolerance_unit_ == "DA") ? fragment_mz_tolerance_ : Math::ppmToMass(fragment_mz_tolerance_, peak.getPeak().getMZ());

    recursiveQuery(hits, peak, peptide_idx_range, window, depth_+1, 0, frag_tol);


  }

  void FragmentIndex3D::recursiveQuery(vector<OpenMS::FragmentIndexTD::Hit>& hits,
                                       const OpenMS::MultiPeak& peak,
                                       std::pair<size_t, size_t> peptide_idx_range,
                                       std::pair<double, double> window,
                                       size_t recursion_step,
                                       size_t current_slice,              // From the last recursiv step. Holds the info in which branch of the tree we are in
                                       double fragment_tolerance)
  {
    vector<double>* current_level;
    double current_query;
    std::pair<double, double> applied_window;

    if(recursion_step == 0){ // last (precursor mz) level of the tree. Push hits into the hits vector
      auto last_slice_start = fi_fragments_.begin() + current_slice * bucketsize_;
      auto last_slice_end = ((current_slice +1) * bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : fi_fragments_.begin() + (current_slice+1)* bucketsize_;
      vector<MultiFragment> last_slice(last_slice_start, last_slice_end);
      auto hits_slice = FragmentIndexTD::binary_search_slice_mf(last_slice,
                                                                   peptide_idx_range.first,
                                                                   peptide_idx_range.second,
                                                                   [](MultiFragment a){return a.getPeptideIdx();},
                                                                    false);
      for(size_t i = hits_slice.first; i <= hits_slice.second; i++)
      {
        if (inRange(last_slice[i].getFragmentMz(), peak.getPeak().getMZ(), fragment_tolerance, window) &&
            inRangeFollowUpPeaks(last_slice[i].getFollowUpPeaks(), peak.getFollowUpPeaks(), fragment_tolerance)) {
          hits.push_back({last_slice[i].getPeptideIdx(), last_slice[i].getFragmentMz()});
        }
      }
      return;
    }
    if(recursion_step == 1){  // Fragment mz level
      current_level = &bucket_min_mz_;
      current_query = peak.getPeak().getMZ();
      applied_window = window;
    }
    if(recursion_step > 1){  // All follow up peak levels
      current_level = &(follow_up_peaks_buckets_min_mz.at(depth_ +1 - recursion_step));   // get the current vector of interests
      current_query = peak.getFollowUpPeaks().at(depth_ +1 -recursion_step);              // current query value of interst
      applied_window = make_pair<double, double>(0, 0);                                 // For the follow up peaks we do not have any window
    }
    auto slice_start = (*current_level).begin() + current_slice * bucketsize_;
    auto slice_end = ((current_slice +1) * bucketsize_) > (*current_level).size() ? (*current_level).end() : (*current_level).begin() + (current_slice+1)* bucketsize_;
    if(recursion_step == depth_+1)            //edge case for the very first tree layer
      slice_end = (*current_level).end();
    vector<double> slice(slice_start, slice_end);
    auto next_slices = FragmentIndexTD::binary_search_slice_double(slice,
                                                           current_query -fragment_tolerance + applied_window.first,
                                                           current_query + fragment_tolerance + applied_window.second,
                                                            true);
    for(size_t next_slice = next_slices.first; next_slice <= next_slices.second; next_slice++){
      recursiveQuery(hits,
                     peak,
                     peptide_idx_range,
                     window,
                     recursion_step-1,
                     current_slice * bucketsize_ + next_slice,                    // first slide to the start position of the window, than add the actual idx of the entry found
                     fragment_tolerance);
    }

  }


}