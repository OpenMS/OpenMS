// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/ML/INTERPOLATION/LinearInterpolation.h>

#include <boost/math/special_functions/fpclassify.hpp> // isnan

// #define Debug_PoseClusteringAffineSuperimposer

namespace OpenMS
{

  PoseClusteringAffineSuperimposer::PoseClusteringAffineSuperimposer() :
    BaseSuperimposer()
  {
    setName("PoseClusteringAffineSuperimposer");

    defaults_.setValue("mz_pair_max_distance", 0.5, "Maximum of m/z deviation of corresponding elements in different maps.  "
                                                    "This condition applies to the pairs considered in hashing.");
    defaults_.setMinFloat("mz_pair_max_distance", 0.);

    defaults_.setValue("rt_pair_distance_fraction", 0.1, "Within each of the two maps, the pairs considered for pose clustering "
                                                         "must be separated by at least this fraction of the total elution time "
                                                         "interval (i.e., max - min).  ", {"advanced"});
    defaults_.setMinFloat("rt_pair_distance_fraction", 0.);
    defaults_.setMaxFloat("rt_pair_distance_fraction", 1.);

    defaults_.setValue("num_used_points", 2000, "Maximum number of elements considered in each map "
                                                "(selected by intensity).  Use this to reduce the running time "
                                                "and to disregard weak signals during alignment.  For using all points, set this to -1.");
    defaults_.setMinInt("num_used_points", -1);

    defaults_.setValue("scaling_bucket_size", 0.005, "The scaling of the retention time "
                                                     "interval is being hashed into buckets of this size during pose "
                                                     "clustering.  A good choice for this would be a bit smaller than the "
                                                     "error you would expect from repeated runs.");
    defaults_.setMinFloat("scaling_bucket_size", 0.);

    defaults_.setValue("shift_bucket_size", 3.0, "The shift at the lower (respectively, higher) end of the retention time "
                                                 "interval is being hashed into buckets of this size during pose "
                                                 "clustering.  A good choice for this would be about "
                                                 "the time between consecutive MS scans.");
    defaults_.setMinFloat("shift_bucket_size", 0.);

    defaults_.setValue("max_shift", 1000.0, "Maximal shift which is considered during histogramming (in seconds).  "
                                            "This applies for both directions.", {"advanced"});
    defaults_.setMinFloat("max_shift", 0.);

    defaults_.setValue("max_scaling", 2.0, "Maximal scaling which is considered during histogramming.  "
                                           "The minimal scaling is the reciprocal of this.", {"advanced"});
    defaults_.setMinFloat("max_scaling", 1.);

    defaults_.setValue("dump_buckets", "", "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  "
                                           "A serial number for each invocation will be appended automatically.", {"advanced"});

    defaults_.setValue("dump_pairs", "", "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  "
                                         "A serial number for each invocation will be appended automatically.", {"advanced"});

    defaultsToParam_();
  }

  /**
    @brief Initialize hash maps for the algorithm.

    The hash maps will contain a histogram of the values to be estimated.

  */
  void initializeHashTables(
    Math::LinearInterpolation<double, double>& scaling_hash_1,
    Math::LinearInterpolation<double, double>& scaling_hash_2,
    Math::LinearInterpolation<double, double>& rt_low_hash_,
    Math::LinearInterpolation<double, double>& rt_high_hash_,
    const double max_scaling, const double max_shift,
    const double scaling_bucket_size, const double shift_bucket_size,
    const double rt_low, const double rt_high)
  {
    const Int scaling_buckets_num_half = (Int) ceil(log(max_scaling) / scaling_bucket_size) + 1;

    // set scale to scaling_bucket_size and establish initial mapping of scaling_buckets_num_half to zero

    scaling_hash_1.getData().clear();
    scaling_hash_1.getData().resize(2 * scaling_buckets_num_half + 1);
    scaling_hash_1.setMapping(scaling_bucket_size, scaling_buckets_num_half, 0.);

    scaling_hash_2.getData().clear();
    scaling_hash_2.getData().resize(2 * scaling_buckets_num_half + 1);
    scaling_hash_2.setMapping(scaling_bucket_size, scaling_buckets_num_half, 0.); // map scaling_buckets_num_half to zero

    // (over)estimate the required number of buckets for shifting
    const Int rt_buckets_num_half = 4 + 2 * (Int) ceil((max_shift * max_scaling) / shift_bucket_size);
    const Int rt_buckets_num = 1 + 2 * rt_buckets_num_half;

    // set scale to shift_bucket_size and establish initial mapping of rt_buckets_num_half to rt_low/rt_high
    rt_low_hash_.getData().clear();
    rt_low_hash_.getData().resize(rt_buckets_num);
    rt_low_hash_.setMapping(shift_bucket_size, rt_buckets_num_half, rt_low);

    rt_high_hash_.getData().clear();
    rt_high_hash_.getData().resize(rt_buckets_num);
    rt_high_hash_.setMapping(shift_bucket_size, rt_buckets_num_half, rt_high);
  }

  /**
    @brief Estimates scaling by trying different (weighted) affine transformations.

    Basically try all combinations of two pairs from map model (i,j) and two
    pairs from map scene (k,l) and compute shift and scale based on these
    four points. The computed value is weighed by the intensity of all
    points, thus this is a density-based approach.

    In the first round, compute and store every combination. In the second
    round, only consider quadruplets where the scaling factor matches the
    estimated bounds of (scale_low_1,scale_high_1), discard all other data.

  */
  void affineTransformationHashing(const bool do_dump_pairs,
                                   const std::vector<Peak2D> & model_map,
                                   const std::vector<Peak2D> & scene_map,
                                   Math::LinearInterpolation<double, double>& scaling_hash_1,
                                   Math::LinearInterpolation<double, double>& scaling_hash_2,
                                   Math::LinearInterpolation<double, double>& rt_low_hash_,
                                   Math::LinearInterpolation<double, double>& rt_high_hash_,
                                   const int hashing_round,
                                   const double rt_pair_min_distance,
                                   const String& dump_pairs_basename,
                                   const Int dump_buckets_serial,
                                   const double mz_pair_max_distance,
                                   const double winlength_factor_baseline,
                                   const double total_intensity_ratio,
                                   const double scale_low_1,
                                   const double scale_high_1,
                                   const double rt_low, const double rt_high)
  {
    Size const model_map_size = model_map.size();   // i j
    Size const scene_map_size = scene_map.size();   // k l

    String dump_pairs_filename;
    std::ofstream dump_pairs_file;
    if (do_dump_pairs)
    {
      dump_pairs_filename = dump_pairs_basename + "_phase_two_" + String(dump_buckets_serial);
      dump_pairs_file.open(dump_pairs_filename.c_str());
      dump_pairs_file << "#" << ' ' << "i" << ' ' << "j" << ' ' << "k" << ' ' << "l" << ' ' << std::endl;
    }

    // first point in model map (i)
    for (Size i = 0, i_low = 0, i_high = 0, k_low = 0, k_high = 0; i < model_map_size - 1; ++i)
    {
      // Adjust window around i in model map (get all features in a m/z range of item i in the model map)
      while (i_low < model_map_size && model_map[i_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance)
        ++i_low;
      while (i_high < model_map_size && model_map[i_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance)
        ++i_high;
      // stop if there are too many features are in our window
      double i_winlength_factor = 1. / (i_high - i_low);
      i_winlength_factor -= winlength_factor_baseline;
      if (i_winlength_factor <= 0)
        continue;

      // Adjust window around k in scene map (get all features in a m/z range of item i in the scene map)
      while (k_low < scene_map_size && scene_map[k_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance)
        ++k_low;
      while (k_high < scene_map_size && scene_map[k_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance)
        ++k_high;

      // Iterate through all matching features in the scene map that are
      // within the m/z distance of item i from the model map.
      // first point in scene map (k)
      for (Size k = k_low; k < k_high; ++k)
      {
        // stop if there are too many features are in our window
        double k_winlength_factor = 1. / (k_high - k_low);
        k_winlength_factor -= winlength_factor_baseline;
        if (k_winlength_factor <= 0)
          continue;

        // compute similarity of intensities i k by taking the ratio of the two intensities
        double similarity_ik;
        {
          const double int_i = model_map[i].getIntensity();
          const double int_k = scene_map[k].getIntensity() * total_intensity_ratio;
          similarity_ik = (int_i < int_k) ? int_i / int_k : int_k / int_i;
          // weight is inverse proportional to number of elements with similar mz
          similarity_ik *= i_winlength_factor;
          similarity_ik *= k_winlength_factor;
        }

        // second point in model map (j)
        for (Size j = i + 1, j_low = i_low, j_high = i_low, l_low = k_low, l_high = k_high; j < model_map_size; ++j)
        {
          // diff in model map -> skip features that are too far away in RT
          double diff_model = model_map[j].getRT() - model_map[i].getRT();
          if (fabs(diff_model) < rt_pair_min_distance)
            continue;

          // Adjust window around j in model map
          while (j_low < model_map_size && model_map[j_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance)
            ++j_low;
          while (j_high < model_map_size && model_map[j_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance)
            ++j_high;
          double j_winlength_factor = 1. / (j_high - j_low);
          j_winlength_factor -= winlength_factor_baseline;
          if (j_winlength_factor <= 0)
            continue;

          // Adjust window around l in scene map
          while (l_low < scene_map_size && scene_map[l_low].getMZ() < model_map[j].getMZ() - mz_pair_max_distance)
            ++l_low;
          while (l_high < scene_map_size && scene_map[l_high].getMZ() <= model_map[j].getMZ() + mz_pair_max_distance)
            ++l_high;

          // second point in scene map (l)
          for (Size l = l_low; l < l_high; ++l)
          {
            double l_winlength_factor = 1. / (l_high - l_low);
            l_winlength_factor -= winlength_factor_baseline;
            if (l_winlength_factor <= 0)
              continue;

            // diff in scene map -> skip features that are too far away in RT
            double diff_scene = scene_map[l].getRT() - scene_map[k].getRT();

            // avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
            // and point pairs with equal retention times (e.g. i_rt == j_rt)
            if (fabs(diff_scene) < rt_pair_min_distance || ((diff_model > 0) != (diff_scene > 0)))
              continue;

            // compute the transformation (i,j) -> (k,l)
            double scaling = diff_model / diff_scene;
            double shift = model_map[i].getRT() - scene_map[k].getRT() * scaling;

            // compute similarity of intensities i k j l
            double similarity_ik_jl;
            {
              // compute similarity of intensities j l
              const double int_j = model_map[j].getIntensity();
              const double int_l = scene_map[l].getIntensity() * total_intensity_ratio;
              double similarity_jl = (int_j < int_l) ? int_j / int_l : int_l / int_j;
              // weight is inverse proportional to number of elements with similar mz
              similarity_jl *= j_winlength_factor;
              similarity_jl *= l_winlength_factor;
              similarity_ik_jl = similarity_ik * similarity_jl;
            }

            // hash the images of scaling, rt_low and rt_high into their respective hash tables
            // store the scaling parameter and the (estimated) transformation of start/end of the maps in hashes
            //   -> in round 2, discard values outside of scale_low_1 and
            //   scale_high_1 (estimated before in scalingEstimate)
            if (hashing_round == 1)
            {
              // hashing round 1 (estimate the scaling only)
              scaling_hash_1.addValue(log(scaling), similarity_ik_jl);
            }
            else if (scaling >= scale_low_1 && scaling <= scale_high_1)
            {
              // hashing round 2 (estimate scaling and shift)
              scaling_hash_2.addValue(log(scaling), similarity_ik_jl);

              const double rt_low_image = shift + rt_low * scaling;
              rt_low_hash_.addValue(rt_low_image, similarity_ik_jl);
              const double rt_high_image = shift + rt_high * scaling;
              rt_high_hash_.addValue(rt_high_image, similarity_ik_jl);

              if (do_dump_pairs)
              {
                dump_pairs_file << i << ' ' << model_map[i].getRT() << ' ' << model_map[i].getMZ() << ' ' << j << ' ' << model_map[j].getRT() << ' '
                                << model_map[j].getMZ() << ' ' << k << ' ' << scene_map[k].getRT() << ' ' << scene_map[k].getMZ() << ' ' << l << ' '
                                << scene_map[l].getRT() << ' ' << scene_map[l].getMZ() << ' ' << similarity_ik_jl << ' ' << std::endl;
              }
            }
          }   // l
        }   // j
      }   // k
    }   // i
  }

  /**
    @brief Estimates likely position of the scale factor based on scaling_hash_1.

    Uses the histogram given in scaling_hash_1 to perform some filtering and
    correction of the data and then estimate the mean of the scaling factor
    and standard deviation, thus returning a likely range for the scaling
    factor. Outliers are removed in an iterative process.

    Returns scale_centroid_1 (mean), scale_low_1 (lower bound) and
    scale_high_1 (upper bound) of the scaling factor.

  */
  void scalingEstimate(
    Math::LinearInterpolation<double, double>& scaling_hash_1,
    const bool do_dump_buckets,
    const UInt struc_elem_length_datapoints,
    const String& dump_buckets_basename,
    const Int dump_buckets_serial,
    const double scaling_histogram_crossing_slope,
    const double scaling_cutoff_stdev_multiplier,
    const UInt loops_mean_stdev_cutoff,
    double& scale_low_1,
    double& scale_high_1,
    double& scale_centroid_1)
  {
    typedef Math::LinearInterpolation<double, double> LinearInterpolationType_;
    UInt filtering_stage = 0;

    // optionally, dump before filtering
    String dump_buckets_filename;
    std::ofstream dump_buckets_file;
    if (do_dump_buckets)
    {
      dump_buckets_filename = dump_buckets_basename + "_scale_" + String(dump_buckets_serial);
      dump_buckets_file.open(dump_buckets_filename.c_str());

      dump_buckets_file << "# rt scale hash table buckets dump ( scale, height ) : " << dump_buckets_filename << std::endl;
      dump_buckets_file << "# unfiltered hash data\n";
      for (Size index = 0; index < scaling_hash_1.getData().size(); ++index)
      {
        const double log_of_scale = scaling_hash_1.index2key(index);
        const double height = scaling_hash_1.getData()[index];
        dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_file << '\n';
    }

    ++filtering_stage;

    // ***************************************************************************
    // Data filtering: apply tophat filter to histogram of different scales
    // ***************************************************************************
    MorphologicalFilter morph_filter;
    Param morph_filter_param;
    morph_filter_param.setValue("struc_elem_unit", "DataPoints");
    morph_filter_param.setValue("struc_elem_length", double(struc_elem_length_datapoints));
    morph_filter_param.setValue("method", "tophat");
    morph_filter.setParameters(morph_filter_param);
    LinearInterpolationType_::container_type buffer(scaling_hash_1.getData().size());
    morph_filter.filterRange(scaling_hash_1.getData().begin(), scaling_hash_1.getData().end(), buffer.begin());
    scaling_hash_1.getData().swap(buffer);

    // optionally, dump after filtering
    if (do_dump_buckets)
    {
      dump_buckets_file << "# tophat filtered hash data\n";
      for (Size index = 0; index < scaling_hash_1.getData().size(); ++index)
      {
        const double log_of_scale = scaling_hash_1.index2key(index);
        const double height = scaling_hash_1.getData()[index];
        dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_file << '\n';
    }

    ++filtering_stage;

    // ***************************************************************************
    // Data cutoff: estimate cutoff for filtered histogram
    // compute freq_cutoff using a fancy criterion to distinguish between the
    // noise level of the histogram and enriched histogram bins
    // ***************************************************************************
    double freq_cutoff;
    do
    {
      std::copy(scaling_hash_1.getData().begin(), scaling_hash_1.getData().end(), buffer.begin());
      std::sort(buffer.begin(), buffer.end(), std::greater<double>());
      double freq_intercept = scaling_hash_1.getData().front();
      double freq_slope = (scaling_hash_1.getData().back() - scaling_hash_1.getData().front()) / double(buffer.size())
                          / scaling_histogram_crossing_slope;
      if (!freq_slope || buffer.empty())
      {
        // in fact these conditions are actually impossible, but let's be really sure ;-)
        freq_cutoff = 0;
      }
      else
      {
        // -> basically trying to find the intersection where sorted values fall
        // below fitted line with slop "freq_slope"
        Size index = 1;   // not 0 (!)
        while (buffer[index] >= freq_intercept + freq_slope * double(index))
        {
          ++index;
        }
        freq_cutoff = buffer[--index];   // note that we have index >= 1
      }
    } while (false);

    // ***************************************************************************
    // apply freq_cutoff, setting smaller values to zero
    for (Size index = 0; index < scaling_hash_1.getData().size(); ++index)
    {
      if (scaling_hash_1.getData()[index] < freq_cutoff)
      {
        scaling_hash_1.getData()[index] = 0;
      }
    }

    // optionally, dump after noise filtering using freq_cutoff
    if (do_dump_buckets)
    {
      dump_buckets_file << "# after freq_cutoff, which is: " << freq_cutoff << '\n';
      for (Size index = 0; index < scaling_hash_1.getData().size(); ++index)
      {
        const double log_of_scale = scaling_hash_1.index2key(index);
        const double height = scaling_hash_1.getData()[index];
        dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_file << '\n';
    }

    // ***************************************************************************
    // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
    // ***************************************************************************
    Math::BasicStatistics<double> statistics;
    std::vector<double>::const_iterator data_begin = scaling_hash_1.getData().begin();
    const Size data_size = scaling_hash_1.getData().size();
    Size data_range_begin = 0;
    Size data_range_end = data_size;
    for (UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop)     // MAGIC ALERT: number of loops
    {
      statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
      double mean = statistics.mean() + data_range_begin;
      double stdev = sqrt(statistics.variance());
      data_range_begin = floor(std::max<double>(mean - scaling_cutoff_stdev_multiplier * stdev, 0));
      data_range_end = ceil(std::min<double>(mean + scaling_cutoff_stdev_multiplier * stdev + 1, data_size));
      const double log_outside_mean = scaling_hash_1.index2key(mean);
      const double log_outside_stdev = stdev * scaling_hash_1.getScale();
      scale_low_1 = exp(log_outside_mean - log_outside_stdev);
      scale_centroid_1 = exp(log_outside_mean);
      scale_high_1 = exp(log_outside_mean + log_outside_stdev);
      if (do_dump_buckets)
      {
        dump_buckets_file << "# loop: " << loop << "  mean: " << log_outside_mean << " [" << exp(log_outside_mean) << "]  stdev: " << log_outside_stdev
                          << " [" << scale_centroid_1 << "]  (mean-stdev): " << log_outside_mean - log_outside_stdev << " [" << scale_low_1 << "]  (mean+stdev): "
                          << log_outside_mean + log_outside_stdev << " [" << scale_high_1 << "]  data_range_begin: " << data_range_begin << "  data_range_end: "
                          << data_range_end << std::endl;
      }
    }

    if (do_dump_buckets)
    {
      dump_buckets_file << "# EOF" << std::endl;
      dump_buckets_file.close();
    }
  }

  void shiftEstimate(
    const bool do_dump_buckets,
    Math::LinearInterpolation<double, double>& rt_low_hash_,
    Math::LinearInterpolation<double, double>& rt_high_hash_,
    const Int dump_buckets_serial,
    const UInt struc_elem_length_datapoints,
    const double scaling_histogram_crossing_slope,
    const double scaling_cutoff_stdev_multiplier,
    const UInt loops_mean_stdev_cutoff,
    const String& dump_buckets_basename,
    double& rt_low_centroid,
    double& rt_high_centroid)
  {
    UInt filtering_stage = 0;

    // optionally, dump before filtering
    String dump_buckets_low_filename;
    std::ofstream dump_buckets_low_file;
    String dump_buckets_high_filename;
    std::ofstream dump_buckets_high_file;
    if (do_dump_buckets)
    {
      dump_buckets_low_filename = dump_buckets_basename + "_low_" + String(dump_buckets_serial);
      dump_buckets_low_file.open(dump_buckets_low_filename.c_str());

      dump_buckets_low_file << "# rt low hash table buckets dump ( scale, height ) : " << dump_buckets_low_filename << std::endl;
      dump_buckets_low_file << "# unfiltered hash data\n";
      for (Size index = 0; index < rt_low_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_low_hash_.index2key(index);
        const double height = rt_low_hash_.getData()[index];
        dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_low_file << '\n';

      dump_buckets_high_filename = dump_buckets_basename + "_high_" + String(dump_buckets_serial);
      dump_buckets_high_file.open(dump_buckets_high_filename.c_str());

      dump_buckets_high_file << "# rt high hash table buckets dump ( scale, height ) : " << dump_buckets_high_filename << std::endl;
      dump_buckets_high_file << "# unfiltered hash data\n";
      for (Size index = 0; index < rt_high_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_high_hash_.index2key(index);
        const double height = rt_high_hash_.getData()[index];
        dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_high_file << '\n';
    }

    ++filtering_stage;

    // apply tophat filter to histogram
    MorphologicalFilter morph_filter;
    Param morph_filter_param;
    morph_filter_param.setValue("struc_elem_unit", "DataPoints");
    morph_filter_param.setValue("struc_elem_length", double(struc_elem_length_datapoints));
    morph_filter_param.setValue("method", "tophat");
    morph_filter.setParameters(morph_filter_param);

    typedef Math::LinearInterpolation<double, double> LinearInterpolationType_;
    LinearInterpolationType_::container_type buffer(rt_low_hash_.getData().size());
    morph_filter.filterRange(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end(), buffer.begin());
    rt_low_hash_.getData().swap(buffer);
    morph_filter.filterRange(rt_high_hash_.getData().begin(), rt_high_hash_.getData().end(), buffer.begin());
    rt_high_hash_.getData().swap(buffer);

    // optionally, dump after filtering
    if (do_dump_buckets)
    {
      dump_buckets_low_file << "# tophat filtered hash data\n";
      for (Size index = 0; index < rt_low_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_low_hash_.index2key(index);
        const double height = rt_low_hash_.getData()[index];
        dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_low_file << '\n';

      dump_buckets_high_file << "# tophat filtered hash data\n";
      for (Size index = 0; index < rt_high_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_high_hash_.index2key(index);
        const double height = rt_high_hash_.getData()[index];
        dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_high_file << '\n';
    }

    ++filtering_stage;

    // compute freq_cutoff using a fancy criterion to distinguish between the noise level of the histogram and enriched histogram bins
    double freq_cutoff_low;
    double freq_cutoff_high;
    do
    {
      {
        std::copy(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end(), buffer.begin());
        std::sort(buffer.begin(), buffer.end(), std::greater<double>());
        double freq_intercept = rt_low_hash_.getData().front();
        double freq_slope = (rt_low_hash_.getData().back() - rt_low_hash_.getData().front()) / double(buffer.size())
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
      {
        std::copy(rt_high_hash_.getData().begin(), rt_high_hash_.getData().end(), buffer.begin());
        std::sort(buffer.begin(), buffer.end(), std::greater<double>());
        double freq_intercept = rt_high_hash_.getData().front();
        double freq_slope = (rt_high_hash_.getData().back() - rt_high_hash_.getData().front()) / double(buffer.size())
                            / scaling_histogram_crossing_slope;
        if (!freq_slope || buffer.empty())
        {
          // in fact these conditions are actually impossible, but let's be really sure ;-)
          freq_cutoff_high = 0;
        }
        else
        {
          Size index = 1; // not 0 (!)
          while (buffer[index] >= freq_intercept + freq_slope * double(index))
          {
            ++index;
          }
          freq_cutoff_high = buffer[--index]; // note that we have index >= 1
        }
      }
    } while (false);

    // apply freq_cutoff, setting smaller values to zero
    for (Size index = 0; index < rt_low_hash_.getData().size(); ++index)
    {
      if (rt_low_hash_.getData()[index] < freq_cutoff_low)
      {
        rt_low_hash_.getData()[index] = 0;
      }
    }
    for (Size index = 0; index < rt_high_hash_.getData().size(); ++index)
    {
      if (rt_high_hash_.getData()[index] < freq_cutoff_high)
      {
        rt_high_hash_.getData()[index] = 0;
      }
    }

    // optionally, dump after noise filtering using freq_cutoff
    if (do_dump_buckets)
    {
      dump_buckets_low_file << "# after freq_cutoff, which is: " << freq_cutoff_low << '\n';
      for (Size index = 0; index < rt_low_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_low_hash_.index2key(index);
        const double height = rt_low_hash_.getData()[index];
        dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_low_file << '\n';

      dump_buckets_high_file << "# after freq_cutoff, which is: " << freq_cutoff_high << '\n';
      for (Size index = 0; index < rt_high_hash_.getData().size(); ++index)
      {
        const double rt_image = rt_high_hash_.index2key(index);
        const double height = rt_high_hash_.getData()[index];
        dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
      }
      dump_buckets_high_file << '\n';
    }

    // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
    {
      Math::BasicStatistics<double> statistics;
      std::vector<double>::const_iterator data_begin = rt_low_hash_.getData().begin();
      const Size data_size = rt_low_hash_.getData().size();
      Size data_range_begin = 0;
      Size data_range_end = data_size;
      for (UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop)   // MAGIC ALERT: number of loops
      {
        statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
        double mean = statistics.mean() + data_range_begin;
        double stdev = sqrt(statistics.variance());
        data_range_begin = floor(std::max<double>(mean - scaling_cutoff_stdev_multiplier * stdev, 0));
        data_range_end = ceil(std::min<double>(mean + scaling_cutoff_stdev_multiplier * stdev + 1, data_size));
        const double outside_mean = rt_low_hash_.index2key(mean);
        const double outside_stdev = stdev * rt_low_hash_.getScale();
        // rt_low_low = (outside_mean - outside_stdev);
        rt_low_centroid = (outside_mean);
        // rt_low_high = (outside_mean + outside_stdev);
        if (do_dump_buckets)
        {
          dump_buckets_low_file << "# loop: " << loop << "  mean: " << outside_mean << "  stdev: " << outside_stdev << "  (mean-stdev): " << outside_mean
          - outside_stdev << "  (mean+stdev): " << outside_mean + outside_stdev << "  data_range_begin: " << data_range_begin << "  data_range_end: "
                                << data_range_end << std::endl;
        }
      }
    }

    // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
    {
      Math::BasicStatistics<double> statistics;
      std::vector<double>::const_iterator data_begin = rt_high_hash_.getData().begin();
      const Size data_size = rt_high_hash_.getData().size();
      Size data_range_begin = 0;
      Size data_range_end = data_size;
      for (UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop)   // MAGIC ALERT: number of loops
      {
        statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
        double mean = statistics.mean() + data_range_begin;
        double stdev = sqrt(statistics.variance());
        data_range_begin = floor(std::max<double>(mean - scaling_cutoff_stdev_multiplier * stdev - 1, 0));
        data_range_end = ceil(std::min<double>(mean + scaling_cutoff_stdev_multiplier * stdev + 2, data_size));
        const double outside_mean = rt_high_hash_.index2key(mean);
        const double outside_stdev = stdev * rt_high_hash_.getScale();
        // rt_high_low = (outside_mean - outside_stdev);
        rt_high_centroid = (outside_mean);
        // rt_high_high = (outside_mean + outside_stdev);
        if (do_dump_buckets)
        {
          dump_buckets_high_file << "# loop: " << loop << "  mean: " << outside_mean << "  stdev: " << outside_stdev << "  (mean-stdev): " << outside_mean
          - outside_stdev << "  (mean+stdev): " << outside_mean + outside_stdev << "  data_range_begin: " << data_range_begin << "  data_range_end: "
                                 << data_range_end << std::endl;
        }
      }
    }
    if (do_dump_buckets)
    {
      dump_buckets_low_file << "# EOF" << std::endl;
      dump_buckets_low_file.close();
      dump_buckets_high_file << "# EOF" << std::endl;
      dump_buckets_high_file.close();
    }
  }

  double computeIntensityRatio(const std::vector<Peak2D> & model_map, const std::vector<Peak2D> & scene_map)
  {
    double total_int_model_map = 0;
    for (Size i = 0; i < model_map.size(); ++i)
    {
      total_int_model_map += model_map[i].getIntensity();
    }

    double total_int_scene_map = 0;
    for (Size i = 0; i < scene_map.size(); ++i)
    {
      total_int_scene_map += scene_map[i].getIntensity();
    }

    return total_int_model_map / total_int_scene_map;
  }

  void PoseClusteringAffineSuperimposer::run(const std::vector<Peak2D> & map_model,
                                             const std::vector<Peak2D> & map_scene, 
                                             TransformationDescription & transformation)

  {
    if (map_model.empty() || map_scene.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "One of the input maps is empty! This is not allowed!");
    }

    //**************************************************************************
    // Parameters
    //**************************************************************************

    // number of data points in structuring element for tophat filter, which removes baseline from histogram
    const UInt struc_elem_length_datapoints = 21;
    // used when distinguishing noise level and enriched histogram bins
    const double scaling_histogram_crossing_slope = 3.0;
    // multiplier for stdev in cutoff for outliers
    const double scaling_cutoff_stdev_multiplier = 1.5;
    // number of loops in stdev cutoff for outliers
    const UInt loops_mean_stdev_cutoff = 3;
    // Each m/z window is given unit weight. If there are too many pairs for a
    // window, the individual contributions will be very small, but running
    // time will be high, so we provide a cutoff for this.  Typically this will
    // exclude compounds which elute over the whole retention time range from
    // consideration.
    // This may lead to the exclusion of certain m/z windows that are very
    // crowded.
    const double winlength_factor_baseline = 0.1;

    /// Maximum deviation in mz of two partner points
    const double mz_pair_max_distance = param_.getValue("mz_pair_max_distance");

    //**************************************************************************
    // Working variables
    //**************************************************************************
    typedef Math::LinearInterpolation<double, double> LinearInterpolationType_;
    // these are a set of hashes that transform bins to actual RT values ...
    LinearInterpolationType_ scaling_hash_1; //scaling estimate from round 1 hashing
    LinearInterpolationType_ scaling_hash_2; //scaling estimate from round 2 hashing
    LinearInterpolationType_ rt_low_hash_; // rt shift estimate of map start
    LinearInterpolationType_ rt_high_hash_; // rt shift estimate of map end
    UInt actual_progress = 0;

    startProgress(0, 100, "affine pose clustering");
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
    // Step 1: Select the most abundant data points only.
    //**************************************************************************
    // use copy to truncate
    std::vector<Peak2D> model_map(map_model);
    std::vector<Peak2D> scene_map(map_scene);
    {
      // truncate the data as necessary
      const Size num_used_points = (Int) param_.getValue("num_used_points");

      // sort the last data points by ascending intensity (from the right, using reverse iterators)
      //  -> linear in complexity, should be faster than sorting and then taking cutoff
      if (model_map.size() > num_used_points)
      {
        std::nth_element(model_map.rbegin(), model_map.rbegin() + (model_map.size() - num_used_points),
            model_map.rend(), Peak2D::IntensityLess());
        model_map.resize(num_used_points);
      }
      setProgress(++actual_progress);
      if (scene_map.size() > num_used_points)
      {
        std::nth_element(scene_map.rbegin(), scene_map.rbegin() + (scene_map.size() - num_used_points),
            scene_map.rend(), Peak2D::IntensityLess());
        scene_map.resize(num_used_points);
      }
      setProgress(++actual_progress);
    }
    // sort by ascending m/z
    std::sort(model_map.begin(), model_map.end(), Peak2D::MZLess());
    std::sort(scene_map.begin(), scene_map.end(), Peak2D::MZLess());
    setProgress((actual_progress = 10));

    //**************************************************************************
    // Preprocessing
    //**************************************************************************
    // take estimates of the minimal / maximal element from both maps
    // possible improvement: use the truncated map from above which should be
    // more reliable (one outlier of low intensity could derail the estimate
    // below)
    const double model_minrt = std::min_element(map_model.begin(), map_model.end(), Peak2D::RTLess())->getRT();
    const double scene_minrt = std::min_element(map_scene.begin(), map_scene.end(), Peak2D::RTLess())->getRT();
    const double model_maxrt = std::max_element(map_model.begin(), map_model.end(), Peak2D::RTLess())->getRT();
    const double scene_maxrt = std::max_element(map_scene.begin(), map_scene.end(), Peak2D::RTLess())->getRT();
    const double rt_low =  (model_minrt + scene_minrt) / 2.;
    const double rt_high = (model_maxrt + scene_maxrt) / 2.;

    //**************************************************************************
    // Sanity check
    //**************************************************************************
    {
      // crude estimate of the shift and slope
      double shift = std::fabs(model_minrt - scene_minrt);
      double slope = (model_maxrt - model_minrt) / (scene_maxrt - scene_minrt);

      if ( (double)param_.getValue("max_scaling") < slope * 1.2 || 
           1.0 / (double)param_.getValue("max_scaling") > slope / 1.2)
      {
        std::cout << "WARNING: your map likely has a scaling around " << slope
          << " but your parameters only allow for a maximal scaling of " <<
          param_.getValue("max_scaling") << std::endl;
        std::cout << "It is strongly advised to adjust your max_scaling factor" << std::endl;
      }

      if ( (double)param_.getValue("max_shift") < shift * 1.2)
      {
        std::cout << "WARNING: your map likely has a shift around " << shift
          << " but your parameters only allow for a maximal shift of " <<
          param_.getValue("max_shift") << std::endl;
        std::cout << "It is strongly advised to adjust your max_shift factor" << std::endl;
      }

    }
    

    // Distance in RT two points need to have at most to be considered for clustering
    const double rt_pair_min_distance = (double) param_.getValue("rt_pair_distance_fraction") * (rt_high - rt_low);

    //**************************************************************************
    // Step 2: Initialize the hash tables: rt_scaling_hash_, rt_low_hash_, and
    //         rt_high_hash_. (over)estimate the required number of buckets
    //         for scaling
    //         Note: the user-specified bucket size only applies to scales
    //         around 1.  The hashing uses a log transformation because we do
    //         not like skewed distributions.
    //**************************************************************************
    initializeHashTables(scaling_hash_1, scaling_hash_2, rt_low_hash_, rt_high_hash_,
                         param_.getValue("max_scaling"), param_.getValue("max_shift"),
                         param_.getValue("scaling_bucket_size"), param_.getValue("shift_bucket_size"),
                         rt_low, rt_high);

    setProgress(++actual_progress);

    //**************************************************************************
    // Step 3: compute the ratio of the total intensities of both maps, for
    //         normalization
    //**************************************************************************
    double total_intensity_ratio = computeIntensityRatio(model_map, scene_map);
    setProgress((actual_progress = 20));

    // The serial number is incremented for each invocation of this, to avoid
    // overwriting of hash table dumps.
    static Int dump_buckets_serial = 0;
    ++dump_buckets_serial;

    //**************************************************************************
    // Step 4: Hashing
    //         Compute the transformations between each point pair in the model
    //         map and each point pair in the scene map and hash the affine
    //         transformation.
    //
    //         To speed up the calculation of the final transformation, we
    //         confine the number of considered point pairs.  We match a point
    //         p in the model map only onto those points p' in the scene map
    //         that lie in a certain mz interval.
    //**************************************************************************

    ///////////////////////////////////////////////////////////////////
    // Step 4.1 First round of hashing: Estimate the scaling
    affineTransformationHashing(
      do_dump_pairs,
      model_map, scene_map,
      scaling_hash_1, scaling_hash_2, rt_low_hash_, rt_high_hash_,
      1,
      rt_pair_min_distance,
      dump_pairs_basename,
      dump_buckets_serial,
      mz_pair_max_distance,
      winlength_factor_baseline,
      total_intensity_ratio,
      -1, // only used in 2nd round of hashing
      -1, // only used in 2nd round of hashing
      rt_low, rt_high);
    setProgress((actual_progress = 30));

    ///////////////////////////////////////////////////////////////////
    // Step 4.2 Estimate the scaling factor (and potential bounds) based on the
    // histogram work on rt_scaling_hash_
    double scale_low_1;
    double scale_centroid_1;
    double scale_high_1;
    scalingEstimate(
      scaling_hash_1,
      do_dump_buckets,
      struc_elem_length_datapoints,
      dump_buckets_basename,
      dump_buckets_serial,
      scaling_histogram_crossing_slope,
      scaling_cutoff_stdev_multiplier,
      loops_mean_stdev_cutoff,
      scale_low_1,
      scale_high_1,
      scale_centroid_1);
    setProgress((actual_progress = 40));

    ///////////////////////////////////////////////////////////////////
    // Step 4.3 Second round of hashing: Estimate the shift at both ends and
    // thereby re-estimate the scaling. This uses the first guess of the
    // scaling to reduce noise in the histograms.
    affineTransformationHashing(
      do_dump_pairs,
      model_map, scene_map,
      scaling_hash_1, scaling_hash_2, rt_low_hash_, rt_high_hash_,
      2,
      rt_pair_min_distance,
      dump_pairs_basename,
      dump_buckets_serial,
      mz_pair_max_distance,
      winlength_factor_baseline,
      total_intensity_ratio,
      scale_low_1,
      scale_high_1,
      rt_low, rt_high);
    setProgress((actual_progress = 50));

    ///////////////////////////////////////////////////////////////////
    // Step 4.4 Estimate the shift factor at start/end of the map based on the
    // histogram work on rt_low_hash_ and rt_high_hash_
    double rt_low_centroid;
    double rt_high_centroid;
    shiftEstimate(
      do_dump_buckets,
      rt_low_hash_,
      rt_high_hash_,
      dump_buckets_serial,
      struc_elem_length_datapoints,
      scaling_histogram_crossing_slope,
      scaling_cutoff_stdev_multiplier,
      loops_mean_stdev_cutoff,
      dump_buckets_basename,
      rt_low_centroid,
      rt_high_centroid);
    setProgress(80);

    //**************************************************************************
    // Step 5: Estimate transform
    //         Compute the shifts at the low and high ends by either using the
    //         estimated centroids from the distribution (new method) or by
    //         looking at (around) the fullest bins (old method).
    //**************************************************************************

    double rt_low_image;
    double rt_high_image;

    // 5.1 use centroids for images of rt_low and rt_high
#if 1
    rt_low_image = rt_low_centroid;
    rt_high_image = rt_high_centroid;
#else
    // Alternative: use maximum bins instead (i.e. most likely shift)
    // (Note: this is a fossil which would disregard most of the above
    // computations! The code is left here for developers/debugging only.)
    // This does not fully take into account the shape of the histogram and may
    // be potentially not be as robust as working on histogram data

    const Size rt_low_max_index = std::distance(rt_low_hash_.getData().begin(),
                                                std::max_element(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end()));
    rt_low_image = rt_low_hash_.index2key(rt_low_max_index);

    const Size rt_high_max_index = std::distance(rt_high_hash_.getData().begin(), std::max_element(rt_high_hash_.getData().begin(),
                                                                                                   rt_high_hash_.getData().end()));
    rt_high_image = rt_high_hash_.index2key(rt_high_max_index);
#endif

    setProgress(++actual_progress);

    // 5.2 compute slope and intercept from matching high/low retention times
    {
      Param params;
      const double slope = ((rt_high_image - rt_low_image) / (rt_high - rt_low));
      params.setValue("slope", slope);

      const double intercept = rt_low_image - rt_low * slope;
      params.setValue("intercept", intercept);

      if (std::isinf(slope) || std::isnan(slope) || std::isinf(intercept) || std::isnan(intercept))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Superimposer could not compute an initial transformation!") +
                                      "You can try to increase 'max_num_peaks_considered' to solve this.", String(intercept * slope));
      }

      transformation.fitModel("linear", params);       // no data, but explicit parameters
    }

    setProgress(++actual_progress);
    endProgress();
  }

  void PoseClusteringAffineSuperimposer::run(const ConsensusMap& map_model,
                                             const ConsensusMap& map_scene,
                                             TransformationDescription& transformation)
  {
    std::vector<Peak2D> c_map_model, c_map_scene;
    for (ConsensusMap::const_iterator it = map_model.begin(); it != map_model.end(); ++it)
    {
      Peak2D c;
      c.setIntensity( it->getIntensity() );
      c.setRT( it->getRT() );
      c.setMZ( it->getMZ() );
      c_map_model.push_back(c);
    }
    for (ConsensusMap::const_iterator it = map_scene.begin(); it != map_scene.end(); ++it)
    {
      Peak2D c;
      c.setIntensity( it->getIntensity() );
      c.setRT( it->getRT() );
      c.setMZ( it->getMZ() );
      c_map_scene.push_back(c);
    }

    run(c_map_model, c_map_scene, transformation);
  }

} // namespace OpenMS
