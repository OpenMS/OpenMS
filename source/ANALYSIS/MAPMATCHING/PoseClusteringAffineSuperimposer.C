// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <boost/math/special_functions/fpclassify.hpp>

// #define Debug_PoseClusteringAffineSuperimposer
#ifdef Debug_PoseClusteringAffineSuperimposer
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla)
#endif
#define VV_(bla) V_(""#bla": " << bla)

namespace OpenMS
{

  PoseClusteringAffineSuperimposer::PoseClusteringAffineSuperimposer() :
    BaseSuperimposer()
  {
    setName(getProductName());

    defaults_.setValue("mz_pair_max_distance", 0.5, "Maximum of m/z deviation of corresponding elements in different maps.  "
      "This condition applies to the pairs considered in hashing.");
    defaults_.setMinFloat("mz_pair_max_distance",0.);

    defaults_.setValue("rt_pair_distance_fraction", 0.1, "Within each of the two maps, the pairs considered for pose clustering "
      "must be separated by at least this fraction of the total elution time "
      "interval (i.e., max - min).  ", StringList::create("advanced"));
    defaults_.setMinFloat("rt_pair_distance_fraction",0.);
    defaults_.setMaxFloat("rt_pair_distance_fraction",1.);

    defaults_.setValue("num_used_points", 2000, "Maximum number of elements considered in each map "
      "(selected by intensity).  Use this to reduce the running time "
      "and to disregard weak signals during alignment.  For using all points, set this to -1.");
    defaults_.setMinInt("num_used_points", -1);

    defaults_.setValue("scaling_bucket_size", 0.005, "The scaling of the retention time "
      "interval is being hashed into buckets of this size during pose "
      "clustering.  A good choice for this would be a bit smaller than the "
      "error you would expect from repeated runs.");
    defaults_.setMinFloat("scaling_bucket_size",0.);

    defaults_.setValue("shift_bucket_size", 3.0, "The shift at the lower (respectively, higher) end of the retention time "
      "interval is being hashed into buckets of this size during pose "
      "clustering.  A good choice for this would be about "
      "the time between consecutive MS scans.");
    defaults_.setMinFloat("shift_bucket_size",0.);

    defaults_.setValue("max_shift", 1000.0, "Maximal shift which is considered during histogramming.  "
      "This applies for both directions.", StringList::create("advanced"));
    defaults_.setMinFloat("max_shift",0.);

    defaults_.setValue("max_scaling", 2.0, "Maximal scaling which is considered during histogramming.  "
      "The minimal scaling is the reciprocal of this.", StringList::create("advanced"));
    defaults_.setMinFloat("max_scaling",1.);

    defaults_.setValue("dump_buckets", "", "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  "
      "A serial number for each invocation will be appended automatically.", StringList::create("advanced"));

    defaults_.setValue("dump_pairs", "", "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  "
      "A serial number for each invocation will be appended automatically.", StringList::create("advanced"));

    defaultsToParam_();
    return;
  }

  void
  PoseClusteringAffineSuperimposer::run(const std::vector<ConsensusMap>& maps, std::vector<TransformationDescription>& transformations)
  {

    if ( maps.size() != 2 )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Exactly two input maps are required");
    }
    
    if ( maps[0].empty() || maps[1].empty() )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "One of the input maps is empty! This is not allowed!");
    }

    typedef ConstRefVector<ConsensusMap> PeakPointerArray_;
    typedef Math::LinearInterpolation<DoubleReal, DoubleReal> LinearInterpolationType_;

    LinearInterpolationType_ scaling_hash_1;
    LinearInterpolationType_ scaling_hash_2;
    LinearInterpolationType_ rt_low_hash_;
    LinearInterpolationType_ rt_high_hash_;

    /// Maximum deviation in mz of two partner points
    const DoubleReal mz_pair_max_distance = param_.getValue("mz_pair_max_distance");

    /// Size of each shift bucket
    const DoubleReal shift_bucket_size = param_.getValue("shift_bucket_size");

    const UInt struc_elem_length_datapoints = 21; // MAGIC ALERT: number of data points in structuring element for tophat filter, which removes baseline from histogram
    const DoubleReal scaling_histogram_crossing_slope = 3.0; // MAGIC ALERT: used when distinguishing noise level and enriched histogram bins
    const DoubleReal scaling_cutoff_stdev_multiplier = 1.5; // MAGIC ALERT: multiplier for stdev in cutoff for outliers
    const UInt loops_mean_stdev_cutoff = 3; // MAGIC ALERT: number of loops in stdev cutoff for outliers

    startProgress(0, 100, "affine pose clustering");
    UInt actual_progress = 0;
    setProgress(++actual_progress);

    // Optionally, we will write dumps of the hash table buckets.
    bool do_dump_buckets = false;
    String dump_buckets_basename;
    if ( param_.getValue("dump_buckets") != "" )
    {
      do_dump_buckets = true;
      dump_buckets_basename = param_.getValue("dump_buckets");
    }
    setProgress(++actual_progress);

    // Even more optionally, we will write dumps of the hashed pairs.
    bool do_dump_pairs = false;
    String dump_pairs_basename;
    if ( param_.getValue("dump_pairs") != "" )
    {
      do_dump_pairs = true;
      dump_pairs_basename = param_.getValue("dump_pairs");
    }
    setProgress(++actual_progress);

    //**************************************************************************
    // Select the most abundant data points only.  After that, disallow modifications
    // (we tend to have annoying issues with const_iterator versus iterator).
    PeakPointerArray_ model_map_ini(maps[0].begin(), maps[0].end());
    const PeakPointerArray_ & model_map(model_map_ini);
    PeakPointerArray_ scene_map_ini(maps[1].begin(), maps[1].end());
    const PeakPointerArray_ & scene_map(scene_map_ini);
    {
      // truncate the data as necessary
      const Size num_used_points = (Int) param_.getValue("num_used_points");
      if ( model_map_ini.size() > num_used_points )
      {
        model_map_ini.sortByIntensity(true);
        model_map_ini.resize(num_used_points);
      }
      model_map_ini.sortByComparator(Peak2D::MZLess());
      setProgress(++actual_progress);
      if ( scene_map_ini.size() > num_used_points )
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
    const DoubleReal rt_low = (maps[0].getMin()[ConsensusFeature::RT] + maps[1].getMin()[ConsensusFeature::RT]) / 2.;
    const DoubleReal rt_high = (maps[0].getMax()[ConsensusFeature::RT] + maps[1].getMax()[ConsensusFeature::RT]) / 2.;

    const DoubleReal rt_pair_min_distance = (DoubleReal) param_.getValue("rt_pair_distance_fraction") * (rt_high - rt_low);

    // Initialize the hash tables: rt_scaling_hash_, rt_low_hash_, and rt_high_hash_
    {
      // (over)estimate the required number of buckets for scaling
      const DoubleReal max_scaling = param_.getValue("max_scaling");
      const DoubleReal max_shift = param_.getValue("max_shift");
      // Note: the user-specified bucket size only applies to scales around 1.  The hashing uses a log transformation because we do not like skewed distributions.
      const DoubleReal scaling_bucket_size = param_.getValue("scaling_bucket_size");
      const Int scaling_buckets_num_half = (Int) ceil(log(max_scaling) / scaling_bucket_size) + 1;

      scaling_hash_1.getData().clear();
      scaling_hash_1.getData().resize(2* scaling_buckets_num_half + 1);
      scaling_hash_1.setMapping(scaling_bucket_size, scaling_buckets_num_half, 0.);

      scaling_hash_2.getData().clear();
      scaling_hash_2.getData().resize(2* scaling_buckets_num_half + 1);
      scaling_hash_2.setMapping(scaling_bucket_size, scaling_buckets_num_half, 0.);

      // (over)estimate the required number of buckets for shifting
      const Int rt_buckets_num_half = 4 + 2 * (Int) ceil((max_shift*max_scaling) / shift_bucket_size);
      const Int rt_buckets_num = 1 + 2 * rt_buckets_num_half;

      rt_low_hash_.getData().clear();
      rt_low_hash_.getData().resize(rt_buckets_num);
      rt_low_hash_.setMapping(shift_bucket_size, rt_buckets_num_half, rt_low);

      rt_high_hash_.getData().clear();
      rt_high_hash_.getData().resize(rt_buckets_num);
      rt_high_hash_.setMapping(shift_bucket_size, rt_buckets_num_half, rt_high);
    }
    setProgress(++actual_progress);

    //**************************************************************************
    // compute the ratio of the total intensities of both maps, for normalization
    DoubleReal total_intensity_ratio;
    do
    {
      DoubleReal total_int_model_map = 0;
      for ( Size i = 0; i < model_map.size(); ++i )
      {
        total_int_model_map += model_map[i].getIntensity();
      }
      setProgress(++actual_progress);
      DoubleReal total_int_scene_map = 0;
      for ( Size i = 0; i < scene_map.size(); ++i )
      {
        total_int_scene_map += scene_map[i].getIntensity();
      }
      setProgress(++actual_progress);
      // ... and finally ...
      total_intensity_ratio = total_int_model_map / total_int_scene_map;
    }
    while ( 0 ); // (the extra syntax helps with code folding in eclipse!)
    setProgress((actual_progress = 20));

    /// The serial number is incremented for each invocation of this, to avoid overwriting of hash table dumps.
    static Int dump_buckets_serial = 0;
    ++dump_buckets_serial;

    //**************************************************************************
    // Hashing

    // Compute the transformations between each point pair in the model map
    // and each point pair in the scene map and hash the affine
    // transformation.

    // To speed up the calculation of the final transformation, we confine the number of
    // considered point pairs.  We match a point p in the model map only onto those points p'
    // in the scene map that lie in a certain mz interval.

    VV_(rt_pair_min_distance);

    Size const model_map_size = model_map.size(); // i j
    Size const scene_map_size = scene_map.size(); // k l

    const DoubleReal winlength_factor_baseline = 0.1; // MAGIC ALERT: Each window is given unit weight.  If there are too many pairs for a window, the individual contributions will be very small, but running time will be high, so we provide a cutoff for this.  Typically this will exclude compounds which elute over the whole retention time range from consideration.


    ///////////////////////////////////////////////////////////////////
    // First round of hashing:  Estimate the scaling

    do // begin of hashing (the extra syntax helps with code folding in eclipse!)
    {
      String dump_pairs_filename;
      std::ofstream dump_pairs_file;
      if ( do_dump_pairs )
      {
        dump_pairs_filename = dump_pairs_basename + "_phase_one_" + String(dump_buckets_serial);
        dump_pairs_file.open(dump_pairs_filename.c_str());
        dump_pairs_file << "#" << ' ' << "i" << ' ' << "j" << ' ' << "k" << ' ' << "l" << ' ' << std::endl;
      }
      setProgress(++actual_progress);

      // first point in model map
      for ( Size i = 0, i_low = 0, i_high = 0, k_low = 0, k_high = 0; i < model_map_size - 1; ++i )
      {
        setProgress(actual_progress + Real(i) / model_map_size * 10.f);

        // Adjust window around i in model map
        while ( i_low < model_map_size && model_map[i_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
          ++i_low;
        while ( i_high < model_map_size && model_map[i_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
          ++i_high;
        DoubleReal i_winlength_factor = 1. / (i_high - i_low);
        i_winlength_factor -= winlength_factor_baseline;
        if ( i_winlength_factor <= 0 )
          continue;

        // Adjust window around k in scene map
        while ( k_low < scene_map_size && scene_map[k_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
          ++k_low;
        while ( k_high < scene_map_size && scene_map[k_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
          ++k_high;

        // first point in scene map
        for ( Size k = k_low; k < k_high; ++k )
        {
          DoubleReal k_winlength_factor = 1. / (k_high - k_low);
          k_winlength_factor -= winlength_factor_baseline;
          if ( k_winlength_factor <= 0 )
            continue;

          // compute similarity of intensities i k
          DoubleReal similarity_ik;
          {
            const DoubleReal int_i = model_map[i].getIntensity();
            const DoubleReal int_k = scene_map[k].getIntensity() * total_intensity_ratio;
            similarity_ik = (int_i < int_k) ? int_i / int_k : int_k / int_i;
            // weight is inverse proportional to number of elements with similar mz
            similarity_ik *= i_winlength_factor;
            similarity_ik *= k_winlength_factor;
            // VV_(int_i<<' '<<int_k<<' '<<int_similarity_ik);
          }

          // second point in model map
          for ( Size j = i + 1, j_low = i_low, j_high = i_low, l_low = k_low, l_high = k_high; j < model_map_size; ++j )
          {
            // diff in model map
            DoubleReal diff_model = model_map[j].getRT() - model_map[i].getRT();
            if ( fabs(diff_model) < rt_pair_min_distance )
              continue;

            // Adjust window around j in model map
            while ( j_low < model_map_size && model_map[j_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
              ++j_low;
            while ( j_high < model_map_size && model_map[j_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
              ++j_high;
            DoubleReal j_winlength_factor = 1. / (j_high - j_low);
            j_winlength_factor -= winlength_factor_baseline;
            if ( j_winlength_factor <= 0 )
              continue;

            // Adjust window in scene map
            while ( l_low < scene_map_size && scene_map[l_low].getMZ() < model_map[j].getMZ() - mz_pair_max_distance )
              ++l_low;
            while ( l_high < scene_map_size && scene_map[l_high].getMZ() <= model_map[j].getMZ() + mz_pair_max_distance )
              ++l_high;

            // second point in scene map
            for ( Size l = l_low; l < l_high; ++l )
            {
              DoubleReal l_winlength_factor = 1. / (l_high - l_low);
              l_winlength_factor -= winlength_factor_baseline;
              if ( l_winlength_factor <= 0 )
                continue;

              // diff in scene map
              DoubleReal diff_scene = scene_map[l].getRT() - scene_map[k].getRT();

              // avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
              // and point pairs with equal retention times (e.g. i_rt == j_rt)
              if ( fabs(diff_scene) < rt_pair_min_distance || ((diff_model > 0) != (diff_scene > 0)) )
                continue;

              // compute the transformation (i,j) -> (k,l)
              DoubleReal scaling = diff_model / diff_scene;
              // DoubleReal shift = model_map[i].getRT() - scene_map[k].getRT() * scaling;

              // compute similarity of intensities i k j l
              DoubleReal similarity_ik_jl;
              {
                // compute similarity of intensities j l
                const DoubleReal int_j = model_map[j].getIntensity();
                const DoubleReal int_l = scene_map[l].getIntensity() * total_intensity_ratio;
                DoubleReal similarity_jl = (int_j < int_l) ? int_j / int_l : int_l / int_j;
                // weight is inverse proportional to number of elements with similar mz
                similarity_jl *= j_winlength_factor;
                similarity_jl *= l_winlength_factor;

                // ... and finally ...
                similarity_ik_jl = similarity_ik * similarity_jl;
                // VV_(int_j<<' '<<int_l<<' '<<int_similarity_ik<<' '<<int_similarity_jl<<' '<<int_similarity);
              }

              // hash the images of scaling, rt_low and rt_high into their respective hash tables
              {
                scaling_hash_1.addValue(log(scaling), similarity_ik_jl);

                ///// This will take place in the second round of hashing!
                //  const DoubleReal rt_low_image = shift + rt_low * scaling;
                //  rt_low_hash_.addValue(rt_low_image, similarity_ik_jl);
                //  const DoubleReal rt_high_image = shift + rt_high * scaling;
                //  rt_high_hash_.addValue(rt_high_image, similarity_ik_jl);
              }

              if ( do_dump_pairs )
              {
                dump_pairs_file << i << ' ' << model_map[i].getRT() << ' ' << model_map[i].getMZ() << ' ' << j << ' ' << model_map[j].getRT() << ' '
                    << model_map[j].getMZ() << ' ' << k << ' ' << scene_map[k].getRT() << ' ' << scene_map[k].getMZ() << ' ' << l << ' '
                    << scene_map[l].getRT() << ' ' << scene_map[l].getMZ() << ' ' << similarity_ik_jl << ' ' << std::endl;
              }
            } // l
          } // j
        } // k
      } // i
    }
    while ( 0 ); // end of hashing (the extra syntax helps with code folding in eclipse!)

    setProgress((actual_progress = 30));

    ///////////////////////////////////////////////////////////////////
    // work on rt_scaling_hash_
    DoubleReal scale_low_1;
    DoubleReal scale_centroid_1;
    DoubleReal scale_high_1;
    do
    {

      UInt filtering_stage = 0;

      // optionally, dump before filtering
      String dump_buckets_filename;
      std::ofstream dump_buckets_file;
      if ( do_dump_buckets )
      {
        dump_buckets_filename = dump_buckets_basename + "_scale_" + String(dump_buckets_serial);
        dump_buckets_file.open(dump_buckets_filename.c_str());
        VV_(dump_buckets_filename);

        dump_buckets_file << "# rt scale hash table buckets dump ( scale, height ) : " << dump_buckets_filename << std::endl;
        dump_buckets_file << "# unfiltered hash data\n";
        for ( Size index = 0; index < scaling_hash_1.getData().size(); ++index )
        {
          const DoubleReal log_of_scale = scaling_hash_1.index2key(index);
          const DoubleReal height = scaling_hash_1.getData()[index];
          dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_file << '\n';
      }

      ++filtering_stage;
      setProgress(++actual_progress);

      // apply tophat filter to histogram
      MorphologicalFilter morph_filter;
      Param morph_filter_param;
      morph_filter_param.setValue("struc_elem_unit", "DataPoints");
      morph_filter_param.setValue("struc_elem_length", DoubleReal(struc_elem_length_datapoints));
      morph_filter_param.setValue("method", "tophat");
      morph_filter.setParameters(morph_filter_param);
      LinearInterpolationType_::container_type buffer(scaling_hash_1.getData().size());
      morph_filter.filterRange(scaling_hash_1.getData().begin(), scaling_hash_1.getData().end(), buffer.begin());
      scaling_hash_1.getData().swap(buffer);
      setProgress(++actual_progress);

      // optionally, dump after filtering
      if ( do_dump_buckets )
      {
        dump_buckets_file << "# tophat filtered hash data\n";
        for ( Size index = 0; index < scaling_hash_1.getData().size(); ++index )
        {
          const DoubleReal log_of_scale = scaling_hash_1.index2key(index);
          const DoubleReal height = scaling_hash_1.getData()[index];
          dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_file << '\n';
      }
      setProgress(++actual_progress);

      ++filtering_stage;

      // compute freq_cutoff using a fancy criterion to distinguish between the noise level of the histogram and enriched histogram bins
      DoubleReal freq_cutoff;
      do
      {
        std::copy(scaling_hash_1.getData().begin(), scaling_hash_1.getData().end(), buffer.begin());
        std::sort(buffer.begin(), buffer.end(), std::greater<DoubleReal>());
        DoubleReal freq_intercept = scaling_hash_1.getData().front();
        DoubleReal freq_slope = (scaling_hash_1.getData().back() - scaling_hash_1.getData().front()) / DoubleReal(buffer.size())
            / scaling_histogram_crossing_slope;
        if ( !freq_slope || !buffer.size() )
        {
          // in fact these conditions are actually impossible, but let's be really sure ;-)
          freq_cutoff = 0;
        }
        else
        {
          Size index = 1; // not 0 (!)
          while ( buffer[index] >= freq_intercept + freq_slope * DoubleReal(index) )
          {
            ++index;
          }
          freq_cutoff = buffer[--index]; // note that we have index >= 1
        }
      }
      while ( 0 );
      setProgress(++actual_progress);

      // apply freq_cutoff, setting smaller values to zero
      for ( Size index = 0; index < scaling_hash_1.getData().size(); ++index )
      {
        if ( scaling_hash_1.getData()[index] < freq_cutoff )
        {
          scaling_hash_1.getData()[index] = 0;
        }
      }
      setProgress(++actual_progress);

      // optionally, dump after noise filtering using freq_cutoff
      if ( do_dump_buckets )
      {
        dump_buckets_file << "# after freq_cutoff, which is: " << freq_cutoff << '\n';
        for ( Size index = 0; index < scaling_hash_1.getData().size(); ++index )
        {
          const DoubleReal log_of_scale = scaling_hash_1.index2key(index);
          const DoubleReal height = scaling_hash_1.getData()[index];
          dump_buckets_file << log_of_scale << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_file << '\n';
      }
      setProgress(++actual_progress);

      // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
      Math::BasicStatistics<DoubleReal> statistics;
      std::vector<DoubleReal>::const_iterator data_begin = scaling_hash_1.getData().begin();
      const Size data_size = scaling_hash_1.getData().size();
      Size data_range_begin = 0;
      Size data_range_end = data_size;
      DoubleReal mean;
      DoubleReal stdev;
      for ( UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop ) // MAGIC ALERT: number of loops
      {
        statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
        mean = statistics.mean() + data_range_begin;
        stdev = sqrt(statistics.variance());
        data_range_begin = floor(std::max<DoubleReal>(mean - scaling_cutoff_stdev_multiplier * stdev, 0));
        data_range_end = ceil(std::min<DoubleReal>(mean + scaling_cutoff_stdev_multiplier * stdev + 1, data_size));
        const DoubleReal log_outside_mean = scaling_hash_1.index2key(mean);
        const DoubleReal log_outside_stdev = stdev * scaling_hash_1.getScale();
        scale_low_1 = exp(log_outside_mean - log_outside_stdev);
        scale_centroid_1 = exp(log_outside_mean);
        scale_high_1 = exp(log_outside_mean + log_outside_stdev);
        if ( do_dump_buckets )
        {
          dump_buckets_file << "# loop: " << loop << "  mean: " << log_outside_mean << " [" << exp(log_outside_mean) << "]  stdev: " << log_outside_stdev
              << " [" << scale_centroid_1 << "]  (mean-stdev): " << log_outside_mean - log_outside_stdev << " [" << scale_low_1 << "]  (mean+stdev): "
              << log_outside_mean + log_outside_stdev << " [" << scale_high_1 << "]  data_range_begin: " << data_range_begin << "  data_range_end: "
              << data_range_end << std::endl;
        }
        setProgress(++actual_progress);
      }

      if ( do_dump_buckets )
      {
        dump_buckets_file << "# EOF" << std::endl;
        dump_buckets_file.close();
      }
    }
    while ( 0 );
    setProgress((actual_progress = 40));

    ///////////////////////////////////////////////////////////////////
    // Second round of hashing:  Estimate the shift at both ends and thereby re-estimate the scaling.  This uses the first guess of the scaling to reduce noise in the histograms.

    do // begin of hashing (the extra syntax helps with code folding in eclipse!)
    {
      String dump_pairs_filename;
      std::ofstream dump_pairs_file;
      if ( do_dump_pairs )
      {
        dump_pairs_filename = dump_pairs_basename + "_phase_two_" + String(dump_buckets_serial);
        dump_pairs_file.open(dump_pairs_filename.c_str());
        dump_pairs_file << "#" << ' ' << "i" << ' ' << "j" << ' ' << "k" << ' ' << "l" << ' ' << std::endl;
      }
      setProgress(++actual_progress);

      // first point in model map
      for ( Size i = 0, i_low = 0, i_high = 0, k_low = 0, k_high = 0; i < model_map_size - 1; ++i )
      {
        setProgress(actual_progress + Real(i) / model_map_size * 10.f);

        // Adjust window around i in model map
        while ( i_low < model_map_size && model_map[i_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
          ++i_low;
        while ( i_high < model_map_size && model_map[i_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
          ++i_high;
        DoubleReal i_winlength_factor = 1. / (i_high - i_low);
        i_winlength_factor -= winlength_factor_baseline;
        if ( i_winlength_factor <= 0 )
          continue;

        // Adjust window around k in scene map
        while ( k_low < scene_map_size && scene_map[k_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
          ++k_low;
        while ( k_high < scene_map_size && scene_map[k_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
          ++k_high;

        // first point in scene map
        for ( Size k = k_low; k < k_high; ++k )
        {
          DoubleReal k_winlength_factor = 1. / (k_high - k_low);
          k_winlength_factor -= winlength_factor_baseline;
          if ( k_winlength_factor <= 0 )
            continue;

          // compute similarity of intensities i k
          DoubleReal similarity_ik;
          {
            const DoubleReal int_i = model_map[i].getIntensity();
            const DoubleReal int_k = scene_map[k].getIntensity() * total_intensity_ratio;
            similarity_ik = (int_i < int_k) ? int_i / int_k : int_k / int_i;
            // weight is inverse proportional to number of elements with similar mz
            similarity_ik *= i_winlength_factor;
            similarity_ik *= k_winlength_factor;
            // VV_(int_i<<' '<<int_k<<' '<<int_similarity_ik);
          }

          // second point in model map
          for ( Size j = i + 1, j_low = i_low, j_high = i_low, l_low = k_low, l_high = k_high; j < model_map_size; ++j )
          {
            // diff in model map
            DoubleReal diff_model = model_map[j].getRT() - model_map[i].getRT();
            if ( fabs(diff_model) < rt_pair_min_distance )
              continue;

            // Adjust window around j in model map
            while ( j_low < model_map_size && model_map[j_low].getMZ() < model_map[i].getMZ() - mz_pair_max_distance )
              ++j_low;
            while ( j_high < model_map_size && model_map[j_high].getMZ() <= model_map[i].getMZ() + mz_pair_max_distance )
              ++j_high;
            DoubleReal j_winlength_factor = 1. / (j_high - j_low);
            j_winlength_factor -= winlength_factor_baseline;
            if ( j_winlength_factor <= 0 )
              continue;

            // Adjust window in scene map
            while ( l_low < scene_map_size && scene_map[l_low].getMZ() < model_map[j].getMZ() - mz_pair_max_distance )
              ++l_low;
            while ( l_high < scene_map_size && scene_map[l_high].getMZ() <= model_map[j].getMZ() + mz_pair_max_distance )
              ++l_high;

            // second point in scene map
            for ( Size l = l_low; l < l_high; ++l )
            {
              DoubleReal l_winlength_factor = 1. / (l_high - l_low);
              l_winlength_factor -= winlength_factor_baseline;
              if ( l_winlength_factor <= 0 )
                continue;

              // diff in scene map
              DoubleReal diff_scene = scene_map[l].getRT() - scene_map[k].getRT();

              // avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
              // and point pairs with equal retention times (e.g. i_rt == j_rt)
              if ( fabs(diff_scene) < rt_pair_min_distance || ((diff_model > 0) != (diff_scene > 0)) )
                continue;

              // compute the transformation (i,j) -> (k,l)
              DoubleReal scaling = diff_model / diff_scene;
              DoubleReal shift = model_map[i].getRT() - scene_map[k].getRT() * scaling;

              // compute similarity of intensities i k j l
              DoubleReal similarity_ik_jl;
              {
                // compute similarity of intensities j l
                const DoubleReal int_j = model_map[j].getIntensity();
                const DoubleReal int_l = scene_map[l].getIntensity() * total_intensity_ratio;
                DoubleReal similarity_jl = (int_j < int_l) ? int_j / int_l : int_l / int_j;
                // weight is inverse proportional to number of elements with similar mz
                similarity_jl *= j_winlength_factor;
                similarity_jl *= l_winlength_factor;

                // ... and finally ...
                similarity_ik_jl = similarity_ik * similarity_jl;
                // VV_(int_j<<' '<<int_l<<' '<<int_similarity_ik<<' '<<int_similarity_jl<<' '<<int_similarity);
              }

              // hash the images of scaling, rt_low and rt_high into their respective hash tables
              if ( scaling >= scale_low_1 && scaling <= scale_high_1 )
              {
                scaling_hash_2.addValue(log(scaling), similarity_ik_jl);

                const DoubleReal rt_low_image = shift + rt_low * scaling;
                rt_low_hash_.addValue(rt_low_image, similarity_ik_jl);
                const DoubleReal rt_high_image = shift + rt_high * scaling;
                rt_high_hash_.addValue(rt_high_image, similarity_ik_jl);

                if ( do_dump_pairs )
                {
                  dump_pairs_file << i << ' ' << model_map[i].getRT() << ' ' << model_map[i].getMZ() << ' ' << j << ' ' << model_map[j].getRT() << ' '
                      << model_map[j].getMZ() << ' ' << k << ' ' << scene_map[k].getRT() << ' ' << scene_map[k].getMZ() << ' ' << l << ' '
                      << scene_map[l].getRT() << ' ' << scene_map[l].getMZ() << ' ' << similarity_ik_jl << ' ' << std::endl;
                }
              }
            } // l
          } // j
        } // k
      } // i
    }
    while ( 0 ); // end of hashing (the extra syntax helps with code folding in eclipse!)

    setProgress((actual_progress = 50));

    ///////////////////////////////////////////////////////////////////
    // work on rt_low_hash_ and rt_high_hash_
    DoubleReal rt_low_low;
    DoubleReal rt_low_centroid;
    DoubleReal rt_low_high;
    DoubleReal rt_high_low;
    DoubleReal rt_high_centroid;
    DoubleReal rt_high_high;
    do
    {

      UInt filtering_stage = 0;

      // optionally, dump before filtering
      String dump_buckets_low_filename;
      std::ofstream dump_buckets_low_file;
      String dump_buckets_high_filename;
      std::ofstream dump_buckets_high_file;
      if ( do_dump_buckets )
      {
        dump_buckets_low_filename = dump_buckets_basename + "_low_" + String(dump_buckets_serial);
        dump_buckets_low_file.open(dump_buckets_low_filename.c_str());
        VV_(dump_buckets_low_filename);

        dump_buckets_low_file << "# rt low hash table buckets dump ( scale, height ) : " << dump_buckets_low_filename << std::endl;
        dump_buckets_low_file << "# unfiltered hash data\n";
        for ( Size index = 0; index < rt_low_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_low_hash_.index2key(index);
          const DoubleReal height = rt_low_hash_.getData()[index];
          dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_low_file << '\n';

        dump_buckets_high_filename = dump_buckets_basename + "_high_" + String(dump_buckets_serial);
        dump_buckets_high_file.open(dump_buckets_high_filename.c_str());
        VV_(dump_buckets_high_filename);

        dump_buckets_high_file << "# rt high hash table buckets dump ( scale, height ) : " << dump_buckets_high_filename << std::endl;
        dump_buckets_high_file << "# unfiltered hash data\n";
        for ( Size index = 0; index < rt_high_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_high_hash_.index2key(index);
          const DoubleReal height = rt_high_hash_.getData()[index];
          dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_high_file << '\n';
      }

      ++filtering_stage;
      setProgress(++actual_progress);

      // apply tophat filter to histogram
      MorphologicalFilter morph_filter;
      Param morph_filter_param;
      morph_filter_param.setValue("struc_elem_unit", "DataPoints");
      morph_filter_param.setValue("struc_elem_length", DoubleReal(struc_elem_length_datapoints));
      morph_filter_param.setValue("method", "tophat");
      morph_filter.setParameters(morph_filter_param);

      LinearInterpolationType_::container_type buffer(rt_low_hash_.getData().size());
      morph_filter.filterRange(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end(), buffer.begin());
      rt_low_hash_.getData().swap(buffer);
      morph_filter.filterRange(rt_high_hash_.getData().begin(), rt_high_hash_.getData().end(), buffer.begin());
      rt_high_hash_.getData().swap(buffer);

      // optionally, dump after filtering
      if ( do_dump_buckets )
      {
        dump_buckets_low_file << "# tophat filtered hash data\n";
        for ( Size index = 0; index < rt_low_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_low_hash_.index2key(index);
          const DoubleReal height = rt_low_hash_.getData()[index];
          dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_low_file << '\n';

        dump_buckets_high_file << "# tophat filtered hash data\n";
        for ( Size index = 0; index < rt_high_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_high_hash_.index2key(index);
          const DoubleReal height = rt_high_hash_.getData()[index];
          dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_high_file << '\n';
      }
      setProgress(++actual_progress);

      ++filtering_stage;

      // compute freq_cutoff using a fancy criterion to distinguish between the noise level of the histogram and enriched histogram bins
      DoubleReal freq_cutoff_low;
      DoubleReal freq_cutoff_high;
      do
      {
        {
          std::copy(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end(), buffer.begin());
          std::sort(buffer.begin(), buffer.end(), std::greater<DoubleReal>());
          DoubleReal freq_intercept = rt_low_hash_.getData().front();
          DoubleReal freq_slope = (rt_low_hash_.getData().back() - rt_low_hash_.getData().front()) / DoubleReal(buffer.size())
              / scaling_histogram_crossing_slope;
          if ( !freq_slope || !buffer.size() )
          {
            // in fact these conditions are actually impossible, but let's be really sure ;-)
            freq_cutoff_low = 0;
          }
          else
          {
            Size index = 1; // not 0 (!)
            while ( buffer[index] >= freq_intercept + freq_slope * DoubleReal(index) )
            {
              ++index;
            }
            freq_cutoff_low = buffer[--index]; // note that we have index >= 1
          }
        }
        setProgress(++actual_progress);
        {
          std::copy(rt_high_hash_.getData().begin(), rt_high_hash_.getData().end(), buffer.begin());
          std::sort(buffer.begin(), buffer.end(), std::greater<DoubleReal>());
          DoubleReal freq_intercept = rt_high_hash_.getData().front();
          DoubleReal freq_slope = (rt_high_hash_.getData().back() - rt_high_hash_.getData().front()) / DoubleReal(buffer.size())
              / scaling_histogram_crossing_slope;
          if ( !freq_slope || !buffer.size() )
          {
            // in fact these conditions are actually impossible, but let's be really sure ;-)
            freq_cutoff_high = 0;
          }
          else
          {
            Size index = 1; // not 0 (!)
            while ( buffer[index] >= freq_intercept + freq_slope * DoubleReal(index) )
            {
              ++index;
            }
            freq_cutoff_high = buffer[--index]; // note that we have index >= 1
          }
        }
      }
      while ( 0 );
      setProgress(++actual_progress);

      // apply freq_cutoff, setting smaller values to zero
      for ( Size index = 0; index < rt_low_hash_.getData().size(); ++index )
      {
        if ( rt_low_hash_.getData()[index] < freq_cutoff_low )
        {
          rt_low_hash_.getData()[index] = 0;
        }
      }
      setProgress(++actual_progress);
      for ( Size index = 0; index < rt_high_hash_.getData().size(); ++index )
      {
        if ( rt_high_hash_.getData()[index] < freq_cutoff_high )
        {
          rt_high_hash_.getData()[index] = 0;
        }
      }
      setProgress(++actual_progress);

      // optionally, dump after noise filtering using freq_cutoff
      if ( do_dump_buckets )
      {
        dump_buckets_low_file << "# after freq_cutoff, which is: " << freq_cutoff_low << '\n';
        for ( Size index = 0; index < rt_low_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_low_hash_.index2key(index);
          const DoubleReal height = rt_low_hash_.getData()[index];
          dump_buckets_low_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_low_file << '\n';
        setProgress(++actual_progress);

        dump_buckets_high_file << "# after freq_cutoff, which is: " << freq_cutoff_high << '\n';
        for ( Size index = 0; index < rt_high_hash_.getData().size(); ++index )
        {
          const DoubleReal rt_image = rt_high_hash_.index2key(index);
          const DoubleReal height = rt_high_hash_.getData()[index];
          dump_buckets_high_file << rt_image << '\t' << height << '\t' << filtering_stage << '\n';
        }
        dump_buckets_high_file << '\n';
      }
      setProgress(++actual_progress);

      // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
      {
        Math::BasicStatistics<DoubleReal> statistics;
        std::vector<DoubleReal>::const_iterator data_begin = rt_low_hash_.getData().begin();
        const Size data_size = rt_low_hash_.getData().size();
        Size data_range_begin = 0;
        Size data_range_end = data_size;
        DoubleReal mean;
        DoubleReal stdev;
        for ( UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop ) // MAGIC ALERT: number of loops
        {
          statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
          mean = statistics.mean() + data_range_begin;
          stdev = sqrt(statistics.variance());
          data_range_begin = floor(std::max<DoubleReal>(mean - scaling_cutoff_stdev_multiplier * stdev, 0));
          data_range_end = ceil(std::min<DoubleReal>(mean + scaling_cutoff_stdev_multiplier * stdev + 1, data_size));
          const DoubleReal outside_mean = rt_low_hash_.index2key(mean);
          const DoubleReal outside_stdev = stdev * rt_low_hash_.getScale();
          rt_low_low = (outside_mean - outside_stdev);
          rt_low_centroid = (outside_mean);
          rt_low_high = (outside_mean + outside_stdev);
          if ( do_dump_buckets )
          {
            dump_buckets_low_file << "# loop: " << loop << "  mean: " << outside_mean << "  stdev: " << outside_stdev << "  (mean-stdev): " << outside_mean
                - outside_stdev << "  (mean+stdev): " << outside_mean + outside_stdev << "  data_range_begin: " << data_range_begin << "  data_range_end: "
                << data_range_end << std::endl;
          }
        }
        setProgress(++actual_progress);
      }

      // iterative cut-off based on mean and stdev - relies upon scaling_cutoff_stdev_multiplier which is a bit hard to set right.
      {
        Math::BasicStatistics<DoubleReal> statistics;
        std::vector<DoubleReal>::const_iterator data_begin = rt_high_hash_.getData().begin();
        const Size data_size = rt_high_hash_.getData().size();
        Size data_range_begin = 0;
        Size data_range_end = data_size;
        DoubleReal mean;
        DoubleReal stdev;
        for ( UInt loop = 0; loop < loops_mean_stdev_cutoff; ++loop ) // MAGIC ALERT: number of loops
        {
          statistics.update(data_begin + data_range_begin, data_begin + data_range_end);
          mean = statistics.mean() + data_range_begin;
          stdev = sqrt(statistics.variance());
          data_range_begin = floor(std::max<DoubleReal>(mean - scaling_cutoff_stdev_multiplier * stdev -1, 0));
          data_range_end = ceil(std::min<DoubleReal>(mean + scaling_cutoff_stdev_multiplier * stdev + 2, data_size));
          const DoubleReal outside_mean = rt_high_hash_.index2key(mean);
          const DoubleReal outside_stdev = stdev * rt_high_hash_.getScale();
          rt_high_low = (outside_mean - outside_stdev);
          rt_high_centroid = (outside_mean);
          rt_high_high = (outside_mean + outside_stdev);
          if ( do_dump_buckets )
          {
            dump_buckets_high_file << "# loop: " << loop << "  mean: " << outside_mean << "  stdev: " << outside_stdev << "  (mean-stdev): " << outside_mean
                - outside_stdev << "  (mean+stdev): " << outside_mean + outside_stdev << "  data_range_begin: " << data_range_begin << "  data_range_end: "
                << data_range_end << std::endl;
          }
        }
        setProgress(++actual_progress);
      }
      if ( do_dump_buckets )
      {
        dump_buckets_low_file << "# EOF" << std::endl;
        dump_buckets_low_file.close();
        dump_buckets_high_file << "# EOF" << std::endl;
        dump_buckets_high_file.close();
      }
      setProgress(80);

    }
    while ( 0 );

    //************************************************************************************
    // Estimate transform

    // Compute the shifts at the low and high ends by looking at (around) the fullest bins.
    DoubleReal rt_low_image;
    DoubleReal rt_high_image;
#if 1 // yes of course, use centroids for images of rt_low and rt_high
    rt_low_image = rt_low_centroid;
    rt_high_image = rt_high_centroid;
#else // ooh, use maximum bins instead (Note: this is a fossil which would disregard most of the above computations!  The code is left here for developers/debugging only.)
    const Size rt_low_max_index = std::distance(rt_low_hash_.getData().begin(),
        std::max_element(rt_low_hash_.getData().begin(), rt_low_hash_.getData().end()));
    rt_low_image = rt_low_hash_.index2key(rt_low_max_index);

    const Size rt_high_max_index = std::distance(rt_high_hash_.getData().begin(), std::max_element(rt_high_hash_.getData().begin(),
            rt_high_hash_.getData().end()));
    rt_high_image = rt_high_hash_.index2key(rt_high_max_index);
#endif

    VV_(rt_low); VV_(rt_low_image); VV_(rt_high); VV_(rt_high_image);

    setProgress(++actual_progress);

    // set trafo
    {
      transformations.clear();

			Param params;
      const DoubleReal slope = ((rt_high_image - rt_low_image) / 
																(rt_high - rt_low));
      params.setValue("slope", slope);

      const DoubleReal intercept = rt_low_image - rt_low * slope;
      params.setValue("intercept", intercept);

      if (boost::math::isinf(slope) || boost::math::isnan(slope) || boost::math::isinf(intercept) || boost::math::isnan(intercept))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Superimposer could not compute an initial transformation! You can try to increase 'max_num_peaks_considered' to solve this.", String(intercept*slope));
      }

			TransformationDescription trafo;
			trafo.fitModel("linear", params); // no data, but explicit parameters
			transformations.push_back(trafo);
    }

    setProgress(++actual_progress);

    endProgress();

    return;
  } // run()

} // namespace OpenMS
