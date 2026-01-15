// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <OpenMS/ML/INTERPOLATION/LinearInterpolation.h>

// #define Debug_PoseClusteringShiftSuperimposer
#ifdef Debug_PoseClusteringShiftSuperimposer
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla)
#endif
#define VV_(bla) V_("" # bla ": " << bla)

namespace OpenMS
{

  PoseClusteringShiftSuperimposer::PoseClusteringShiftSuperimposer() :
    BaseSuperimposer()
  {
    setName("PoseClusteringShiftSuperimposer");

    defaults_.setValue("mz_pair_max_distance", 0.5, "Maximum of m/z deviation of corresponding elements in different maps.  "
                                                    "This condition applies to the pairs considered in hashing.");
    defaults_.setMinFloat("mz_pair_max_distance", 0.);

    defaults_.setValue("num_used_points", 2000, "Maximum number of elements considered in each map "
                                                "(selected by intensity).  Use this to reduce the running time "
                                                "and to disregard weak signals during alignment.  For using all points, set this to -1.");
    defaults_.setMinInt("num_used_points", -1);

    defaults_.setValue("shift_bucket_size", 3.0, "The shift of the retention time "
                                                 "interval is being hashed into buckets of this size during pose "
                                                 "clustering.  A good choice for this would be about "
                                                 "the time between consecutive MS scans.");
    defaults_.setMinFloat("shift_bucket_size", 0.);

    defaults_.setValue("max_shift", 1000.0, "Maximal shift which is considered during histogramming.  "
                                            "This applies for both directions.", {"advanced"});
    defaults_.setMinFloat("max_shift", 0.);

    defaults_.setValue("dump_buckets", "", "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  "
                                           "A serial number for each invocation will be appended automatically.", {"advanced"});

    defaults_.setValue("dump_pairs", "", "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  "
                                         "A serial number for each invocation will be appended automatically.", {"advanced"});

    defaultsToParam_();
    return;
  }

  void PoseClusteringShiftSuperimposer::run(const ConsensusMap & map_model, const ConsensusMap & map_scene, TransformationDescription & transformation)
  {
    typedef ConstRefVector<ConsensusMap> PeakPointerArray_;
    typedef Math::LinearInterpolation<double, double> LinearInterpolationType_;

    LinearInterpolationType_ shift_hash_;

    // OLD STUFF
    //    LinearInterpolationType_ scaling_hash_1;
    //    LinearInterpolationType_ scaling_hash_2;
    //    LinearInterpolationType_ shift_hash_;
    //    LinearInterpolationType_ rt_high_hash_;

    /// Maximum deviation in mz of two partner points
    const double mz_pair_max_distance = param_.getValue("mz_pair_max_distance");

    /// Size of each shift bucket
    const double shift_bucket_size = param_.getValue("shift_bucket_size");

    const UInt struc_elem_length_datapoints = 21; // MAGIC ALERT: number of data points in structuring element for tophat filter, which removes baseline from histogram
    const double scaling_histogram_crossing_slope = 3.0; // MAGIC ALERT: used when distinguishing noise level and enriched histogram bins
    const double scaling_cutoff_stdev_multiplier = 1.5; // MAGIC ALERT: multiplier for stdev in cutoff for outliers
    const UInt loops_mean_stdev_cutoff = 3; // MAGIC ALERT: number of loops in stdev cutoff for outliers

    startProgress(0, 100, "shift pose clustering");
    UInt actual_progress = 0;
    setProgress(++actual_progress);

    // Optionally, we will write dumps of the hash table buckets.
    bool do_dump_buckets = false;
    String dump_buckets_basename;
    if (param_.getValue("dump_buckets") != "")
    {
      do_dump_buckets = true;
      dump_buckets_basename = param_.getValue("dump_buckets").toString();
    }
    setProgress(++actual_progress);

    // Even more optionally, we will write dumps of the hashed pairs.
    bool do_dump_pairs = false;
    String dump_pairs_basename;
    if (param_.getValue("dump_pairs") != "")
    {
      do_dump_pairs = true;
      dump_pairs_basename = param_.getValue("dump_pairs").toString();
    }
    setProgress(++actual_progress);

    //**************************************************************************
    // Select the most abundant data points only.  After that, disallow modifications
    // (we tend to have annoying issues with const_iterator versus iterator).
    PeakPointerArray_ model_map_ini(map_model.begin(), map_model.end());
    const PeakPointerArray_ & model_map(model_map_ini);
    PeakPointerArray_ scene_map_ini(map_scene.begin(), map_scene.end());
    const PeakPointerArray_ & scene_map(scene_map_ini);
    {
      // truncate the data as necessary
      // casting to SignedSize is done on PURPOSE here! (num_used_points will be maximal if -1 is used)
      const Size num_used_points = (SignedSize) param_.getValue("num_used_points");
      if (model_map_ini.size() > num_used_points)
      {
        model_map_ini.sortByIntensity(true);
        model_map_ini.resize(num_used_points);
      }
      model_map_ini.sortByComparator(Peak2D::MZLess());
      setProgress(++actual_progress);
      if (scene_map_ini.size() > num_used_points)
      {
        scene_map_ini.sortByIntensity(true);
        scene_map_ini.resize(num_used_points);
      }
      scene_map_ini.sortByComparator(Peak2D::MZLess());
      setProgress(++actual_progress);
      // Note: model_map_ini and scene_map_ini will not be used further below
    }
    setProgress((actual_progress = 10));

    //**************************************************************************
    // Preprocessing

    // get RT ranges (NOTE: we trust that min and max have been updated in the
    // ConsensusMap::convert() method !)

    const double model_low = map_model.getMinRT();
    const double scene_low = map_scene.getMinRT();
    const double model_high = map_model.getMaxRT();
    const double scene_high = map_scene.getMaxRT();

    // OLD STUFF
    //    const double rt_low = (maps[0].getMin()[ConsensusFeature::RT] + maps[1].getMin()[ConsensusFeature::RT]) / 2.;
    //    const double rt_high = (maps[0].getMax()[ConsensusFeature::RT] + maps[1].getMax()[ConsensusFeature::RT]) / 2.;

    // Initialize the hash tables: shift_hash_
    // OLD STUFF: was:  rt_scaling_hash_, rt_low_hash_, and rt_high_hash_
    {
      // (over)estimate the required number of buckets for shifting
      double max_shift = param_.getValue("max_shift");
      // actually the largest possible shift can be much smaller, depending on the data
      do
      {
        if (max_shift < 0)
          max_shift = -max_shift;
        //     ...ml@@@mh........    ,    ........ml@@@mh...
        //     ........sl@@@sh...    ,    ...sl@@@sh........
        double diff;
        diff = model_high - scene_low;
        if (diff < 0)
          diff = -diff;
        if (max_shift > diff)
          max_shift = diff;
        diff = model_low - scene_high;
        if (diff < 0)
          diff = -diff;
        if (max_shift > diff)
          max_shift = diff;
      } while (false);

      const Int shift_buckets_num_half = 4 + (Int) ceil((max_shift) / shift_bucket_size);
      const Int shift_buckets_num = 1 + 2 * shift_buckets_num_half;

      shift_hash_.getData().clear();
      shift_hash_.getData().resize(shift_buckets_num);
      shift_hash_.setMapping(shift_bucket_size, shift_buckets_num_half, 0);
    }
    setProgress(++actual_progress);

    //**************************************************************************
    // compute the ratio of the total intensities of both maps, for normalization
    double total_intensity_ratio;
    do
    {
      double total_int_model_map = 0;
      for (Size i = 0; i < model_map.size(); ++i)
      {
        total_int_model_map += model_map[i].getIntensity();
      }
      setProgress(++actual_progress);
      double total_int_scene_map = 0;
      for (Size i = 0; i < scene_map.size(); ++i)
      {
        total_int_scene_map += scene_map[i].getIntensity();
      }
      setProgress(++actual_progress);
      // ... and finally ...
      total_intensity_ratio = total_int_model_map / total_int_scene_map;
    } while (false);   // (the extra syntax helps with code folding in eclipse!)
    setProgress((actual_progress = 20));

    /// The serial number is incremented for each invocation of this, to avoid overwriting of hash table dumps.
    static Int dump_buckets_serial = 0;
    ++dump_buckets_serial;

    //**************************************************************************
    // Hashing

    // Compute the transformations between each point pair in the model map
    // and each point pair in the scene map and hash the shift
    // transformation.

    // To speed up the calculation of the final transformation, we confine the number of
    // considered point pairs.  We match a point p in the model map only onto those points p'
    // in the scene map that lie in a certain mz interval.

    Size const model_map_size = model_map.size(); // i  /* OLD STUFF: also: j */
    Size const scene_map_size = scene_map.size(); // k  /* OLD STUFF: also: l */

    const double winlength_factor_baseline = 0.1; // MAGIC ALERT: Each window is given unit weight.  If there are too many pairs for a window, the individual contributions will be very small, but running time will be high, so we provide a cutoff for this.  Typically this will exclude compounds which elute over the whole retention time range from consideration.


    ///////////////////////////////////////////////////////////////////
    // Hashing:  Estimate the shift

    do // begin of hashing (the extra syntax helps with code folding in eclipse!)
    {
      String dump_pairs_filename;
      std::ofstream dump_pairs_file;
      if (do_dump_pairs)
      {
        dump_pairs_filename = dump_pairs_basename + String(dump_buckets_serial);
        dump_pairs_file.open(dump_pairs_filename.c_str());
        dump_pairs_file << "#" << ' ' << "i" << ' ' << "k" << std::endl;
      }
      setProgress(++actual_progress);

      // first point in model map
      for (Size i = 0, i_low = 0, i_high = 0, k_low = 0, k_high = 0; i < model_map_size - 1; ++i)
      {
        setProgress(actual_progress + float(i) / model_map_size * 10.f);

        // Adjust window around i in model map
        while (i_low < model_map_size && model_map[i_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance)
          ++i_low;
        while (i_high < model_map_size && model_map[i_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance)
          ++i_high;
        double i_winlength_factor = 1. / (i_high - i_low);
        i_winlength_factor -= winlength_factor_baseline;
        if (i_winlength_factor <= 0)
          continue;

        // Adjust window around k in scene map
        while (k_low < scene_map_size && scene_map[k_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance)
          ++k_low;
        while (k_high < scene_map_size && scene_map[k_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance)
          ++k_high;

        // first point in scene map
        for (Size k = k_low; k < k_high; ++k)
        {
          double k_winlength_factor = 1. / (k_high - k_low);
          k_winlength_factor -= winlength_factor_baseline;
          if (k_winlength_factor <= 0)
            continue;

          // compute similarity of intensities i k
          double similarity_ik;
          {
            const double int_i = model_map[i].getIntensity();
            const double int_k = scene_map[k].getIntensity() * total_intensity_ratio;
            similarity_ik = (int_i < int_k) ? int_i / int_k : int_k / int_i;
            // weight is inverse proportional to number of elements with similar mz
            similarity_ik *= i_winlength_factor;
            similarity_ik *= k_winlength_factor;
            // VV_(int_i<<' '<<int_k<<' '<<int_similarity_ik);
          }

          // compute the transformation (i) -> (k)
          double shift = model_map[i].getRT() - scene_map[k].getRT();

          // hash the images of scaling, rt_low and rt_high into their respective hash tables
          shift_hash_.addValue(shift, similarity_ik);

          if (do_dump_pairs)
          {
            dump_pairs_file << i << ' ' << model_map[i].getRT() << ' ' << model_map[i].getMZ() << ' ' << k << ' ' << scene_map[k].getRT() << ' '
                            << scene_map[k].getMZ() << ' ' << similarity_ik << ' ' << std::endl;
          }

        } // k
      } // i
    } while (false);   // end of hashing (the extra syntax helps with code folding in eclipse!)

    setProgress((actual_progress = 30));

    ///////////////////////////////////////////////////////////////////
    // work on shift_hash_
    //   double shift_low;
    //   double shift_centroid;
    //   double shift_high;

    // OLD STUFF
    // double shift_low;
    double shift_centroid(0.0);
    // double shift_high;
    do
    {

      UInt filtering_stage = 0;

      // optionally, dump before filtering
      String dump_buckets_filename;
      std::ofstream dump_buckets_file;
      if (do_dump_buckets)
      {
        dump_buckets_filename = dump_buckets_basename + "_" + String(dump_buckets_serial);
        dump_buckets_file.open(dump_buckets_filename.c_str());
        VV_(dump_buckets_filename);

        dump_buckets_file << "# shift hash table buckets dump ( scale, height ) : " << dump_buckets_filename << std::endl;
        dump_buckets_file << "# unfiltered hash data\n";
        for (Size index = 0; index < shift_hash_.getData().size(); ++index)
        {
          const double image = shift_hash_.index2key(index);
          const double height = shift_hash_.getData()[index];
          dump_buckets_file << filtering_stage << '\t' << index << '\t' << image << '\t' << height << '\n';
        }
        dump_buckets_file << '\n';
      }

      ++filtering_stage;
      setProgress(++actual_progress);

      // apply tophat filter to histogram
      MorphologicalFilter morph_filter;
      Param morph_filter_param;
      morph_filter_param.setValue("struc_elem_unit", "DataPoints");
      morph_filter_param.setValue("struc_elem_length", double(struc_elem_length_datapoints));
      morph_filter_param.setValue("method", "tophat");
      morph_filter.setParameters(morph_filter_param);

      LinearInterpolationType_::container_type buffer(shift_hash_.getData().size());
      morph_filter.filterRange(shift_hash_.getData().begin(), shift_hash_.getData().end(), buffer.begin());
      shift_hash_.getData().swap(buffer);

      // optionally, dump after filtering
      if (do_dump_buckets)
      {
        dump_buckets_file << "# tophat filtered hash data\n";
        for (Size index = 0; index < shift_hash_.getData().size(); ++index)
        {
          const double image = shift_hash_.index2key(index);
          const double height = shift_hash_.getData()[index];
          dump_buckets_file << filtering_stage << '\t' << index << '\t' << image << '\t' << height << '\n';
        }
        dump_buckets_file << '\n';
      }
      setProgress(++actual_progress);

      ++filtering_stage;

      // compute freq_cutoff using a fancy criterion to distinguish between the noise level of the histogram and enriched histogram bins
      double freq_cutoff_low;
      do
      {
        {
          std::copy(shift_hash_.getData().begin(), shift_hash_.getData().end(), buffer.begin());
          std::sort(buffer.begin(), buffer.end(), std::greater<double>());
          double freq_intercept = shift_hash_.getData().front();
          double freq_slope = (shift_hash_.getData().back() - shift_hash_.getData().front()) / double(buffer.size())
                                  / scaling_histogram_crossing_slope;
          if (!freq_slope || buffer.empty())
          {
            // in fact these conditions are actually impossible, but let's be really sure ;-)
            freq_cutoff_low = 0;
          }
          else
          {
            Size index = 1; // not 0 (!)
            while (buffer[index] >= freq_intercept + freq_slope * double(index))
            {
              ++index;
            }
            freq_cutoff_low = buffer[--index]; // note that we have index >= 1
          }
        }
      } while (false);
      setProgress(++actual_progress);

      // apply freq_cutoff, setting smaller values to zero
      for (Size index = 0; index < shift_hash_.getData().size(); ++index)
      {
        if (shift_hash_.getData()[index] < freq_cutoff_low)
        {
          shift_hash_.getData()[index] = 0;
        }
      }
      setProgress(++actual_progress);

      // optionally, dump after noise filtering using freq_cutoff
      if (do_dump_buckets)
      {
        dump_buckets_file << "# after freq_cutoff, which is: " << freq_cutoff_low << '\n';
        for (Size index = 0; index < shift_hash_.getData().size(); ++index)
        {
          const double image = shift_hash_.index2key(index);
          const double height = shift_hash_.getData()[index];
          dump_buckets_file << filtering_stage << '\t' << index << '\t' << image << '\t' << height << '\n';
        }
        dump_buckets_file << '\n';
      }
      setProgress(++actual_progress);

      // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
      {
        Math::BasicStatistics<double> statistics;
        std::vector<double>::const_iterator data_begin = shift_hash_.getData().begin();
        const Size data_size = shift_hash_.getData().size();
        Size data_range_begin = 0;
        Size data_range_end = data_size;
        for (UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop)   // MAGIC ALERT: number of loops
        {
          statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
          double mean = statistics.mean() + data_range_begin;
          double stdev = sqrt(statistics.variance());
          data_range_begin = floor(std::max<double>(mean - scaling_cutoff_stdev_multiplier * stdev, 0));
          data_range_end = ceil(std::min<double>(mean + scaling_cutoff_stdev_multiplier * stdev + 1, data_size));
          const double outside_mean = shift_hash_.index2key(mean);
          const double outside_stdev = stdev * shift_hash_.getScale();
          // shift_low = (outside_mean - outside_stdev);
          shift_centroid = (outside_mean);
          // shift_high = (outside_mean + outside_stdev);
          if (do_dump_buckets)
          {
            dump_buckets_file << "# loop: " << loop << "  mean: " << outside_mean << "  stdev: " << outside_stdev << "  (mean-stdev): "
                              << outside_mean - outside_stdev << "  (mean+stdev): " << outside_mean + outside_stdev
                              << "  data_range_begin: " << data_range_begin << "  data_range_end: "
                              << data_range_end << std::endl;
          }
        }
        setProgress(++actual_progress);
      }
      if (do_dump_buckets)
      {
        dump_buckets_file << "# EOF" << std::endl;
        dump_buckets_file.close();
      }
      setProgress(80);

    } while (false);

    //************************************************************************************
    // Estimate transform

    // Compute the shifts at the low and high ends by looking at (around) the fullest bins.
    double intercept;
#if 1 // yes of course, use centroids for images of rt_low and rt_high
    intercept = shift_centroid;
#else // ooh, use maximum bins instead (Note: this is a fossil which would disregard most of the above computations!  The code is left here for developers/debugging only.)
    const Size rt_low_max_index = std::distance(shift_hash_.getData().begin(),
                                                std::max_element(shift_hash_.getData().begin(), shift_hash_.getData().end()));
    intercept = shift_hash_.index2key(rt_low_max_index);
#endif

    VV_(intercept);

    setProgress(++actual_progress);

    // set trafo
    {
      Param params;
      params.setValue("slope", 1.0);
      params.setValue("intercept", intercept);

      TransformationDescription trafo;
      trafo.fitModel("linear", params);
      transformation = trafo;
    }

    setProgress(++actual_progress);
    endProgress();

    return;
  } // run()

} // namespace OpenMS
