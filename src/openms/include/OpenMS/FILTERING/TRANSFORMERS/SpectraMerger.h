// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch, Lars Nilse $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
#define OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <vector>

namespace OpenMS
{

  /**
  @brief Merges blocks of MS or MS2 spectra

  Parameter's are accessible via the DefaultParamHandler.

      @htmlinclude OpenMS_SpectraMerger.parameters

*/
  class OPENMS_DLLAPI SpectraMerger :
    public DefaultParamHandler, public ProgressLogger
  {

protected:

    /* Determine distance between two spectra

      Distance is determined as

        (d_rt/rt_max_ + d_mz/mz_max_) / 2

    */
    class SpectraDistance_ :
      public DefaultParamHandler
    {
public:
      SpectraDistance_() :
        DefaultParamHandler("SpectraDistance")
      {
        defaults_.setValue("rt_tolerance", 10.0, "Maximal RT distance (in [s]) for two spectra's precursors.");
        defaults_.setValue("mz_tolerance", 1.0, "Maximal m/z distance (in Da) for two spectra's precursors.");
        defaultsToParam_();
      }

      void updateMembers_() override
      {
        rt_max_ = (double) param_.getValue("rt_tolerance");
        mz_max_ = (double) param_.getValue("mz_tolerance");

        return;
      }

      double getSimilarity(const double d_rt, const double d_mz) const
      {
        //     1 - distance
        return 1 - ((d_rt / rt_max_ + d_mz / mz_max_) / 2);
      }

      // measure of SIMILARITY (not distance, i.e. 1-distance)!!
      double operator()(const BaseFeature& first, const BaseFeature& second) const
      {
        // get RT distance:
        double d_rt = fabs(first.getRT() - second.getRT());
        double d_mz = fabs(first.getMZ() - second.getMZ());

        if (d_rt > rt_max_ || d_mz > mz_max_)
        {
          return 0;
        }

        // calculate similarity (0-1):
        double sim = getSimilarity(d_rt, d_mz);

        return sim;
      }

protected:
      double rt_max_;
      double mz_max_;

    }; // end of SpectraDistance

public:

    /// blocks of spectra (master-spectrum index to sacrifice-spectra(the ones being merged into the master-spectrum))
    typedef Map<Size, std::vector<Size> > MergeBlocks;

    /// blocks of spectra (master-spectrum index to update to spectra to average over)
    typedef Map<Size, std::vector<std::pair<Size, double> > > AverageBlocks;

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectraMerger();

    /// copy constructor
    SpectraMerger(const SpectraMerger& source);

    /// destructor
    ~SpectraMerger() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    SpectraMerger& operator=(const SpectraMerger& source);
    // @}

    // @name Merging functions
    // @{
    ///
    template <typename MapType>
    void mergeSpectraBlockWise(MapType& exp)
    {
      IntList ms_levels = param_.getValue("block_method:ms_levels");
      Int rt_block_size(param_.getValue("block_method:rt_block_size"));
      double rt_max_length = (param_.getValue("block_method:rt_max_length"));

      if (rt_max_length == 0) // no rt restriction set?
      {
        rt_max_length = (std::numeric_limits<double>::max)(); // set max rt span to very large value
      }

      for (IntList::iterator it_mslevel = ms_levels.begin(); it_mslevel < ms_levels.end(); ++it_mslevel)
      {
        MergeBlocks spectra_to_merge;
        Size idx_block(0);
        SignedSize block_size_count(rt_block_size + 1);
        Size idx_spectrum(0);
        for (typename MapType::const_iterator it1 = exp.begin(); it1 != exp.end(); ++it1)
        {
          if (Int(it1->getMSLevel()) == *it_mslevel)
          {
            // block full if it contains a maximum number of scans or if maximum rt length spanned
            if (++block_size_count >= rt_block_size ||
                exp[idx_spectrum].getRT() - exp[idx_block].getRT() > rt_max_length)
            {
              block_size_count = 0;
              idx_block = idx_spectrum;
            }
            else
            {
              spectra_to_merge[idx_block].push_back(idx_spectrum);
            }
          }

          ++idx_spectrum;
        }
        // check if last block had sacrifice spectra
        if (block_size_count == 0) //block just got initialized
        {
          spectra_to_merge[idx_block] = std::vector<Size>();
        }

        // merge spectra, remove all old MS spectra and add new consensus spectra
        mergeSpectra_(exp, spectra_to_merge, *it_mslevel);
      }

      exp.sortSpectra();

      return;
    }

    /// merges spectra with similar precursors (must have MS2 level)
    template <typename MapType>
    void mergeSpectraPrecursors(MapType& exp)
    {

      // convert spectra's precursors to clusterizable data
      Size data_size;
      std::vector<BinaryTreeNode> tree;
      Map<Size, Size> index_mapping;
      // local scope to save memory - we do not need the clustering stuff later
      {
        std::vector<BaseFeature> data;

        for (Size i = 0; i < exp.size(); ++i)
        {
          if (exp[i].getMSLevel() != 2)
          {
            continue;
          }

          // remember which index in distance data ==> experiment index
          index_mapping[data.size()] = i;

          // make cluster element
          BaseFeature bf;
          bf.setRT(exp[i].getRT());
          std::vector<Precursor> pcs = exp[i].getPrecursors();
          if (pcs.empty())
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Scan #") + String(i) + " does not contain any precursor information! Unable to cluster!");
          }
          if (pcs.size() > 1)
          {
            LOG_WARN << "More than one precursor found. Using first one!" << std::endl;
          }
          bf.setMZ(pcs[0].getMZ());
          data.push_back(bf);
        }
        data_size = data.size();

        SpectraDistance_ llc;
        llc.setParameters(param_.copy("precursor_method:", true));
        SingleLinkage sl;
        DistanceMatrix<float> dist; // will be filled
        ClusterHierarchical ch;

        //ch.setThreshold(0.99);
        // clustering ; threshold is implicitly at 1.0, i.e. distances of 1.0 (== similarity 0) will not be clustered
        ch.cluster<BaseFeature, SpectraDistance_>(data, llc, sl, tree, dist);
      }

      // extract the clusters
      ClusterAnalyzer ca;
      std::vector<std::vector<Size> > clusters;
      // count number of real tree nodes (not the -1 ones):
      Size node_count = 0;
      for (Size ii = 0; ii < tree.size(); ++ii)
      {
        if (tree[ii].distance >= 1)
        {
          tree[ii].distance = -1;  // manually set to disconnect, as SingleLinkage does not support it
        }
        if (tree[ii].distance != -1)
        {
          ++node_count;
        }
      }
      ca.cut(data_size - node_count, tree, clusters);

      //std::cerr << "Treesize: " << (tree.size()+1) << "   #clusters: " << clusters.size() << std::endl;
      //std::cerr << "tree:\n" << ca.newickTree(tree, true) << "\n";

      // convert to blocks
      MergeBlocks spectra_to_merge;

      for (Size i_outer = 0; i_outer < clusters.size(); ++i_outer)
      {
        if (clusters[i_outer].size() <= 1)
        {
          continue;
        }
        // init block with first cluster element
        Size cl_index0 = clusters[i_outer][0];
        spectra_to_merge[index_mapping[cl_index0]] = std::vector<Size>();
        // add all other elements
        for (Size i_inner = 1; i_inner < clusters[i_outer].size(); ++i_inner)
        {
          Size cl_index = clusters[i_outer][i_inner];
          spectra_to_merge[index_mapping[cl_index0]].push_back(index_mapping[cl_index]);
        }
      }

      // do it
      mergeSpectra_(exp, spectra_to_merge, 2);

      exp.sortSpectra();

      return;
    }

    /**
     * @brief average over neighbouring spectra
     *
     * @param exp   experimental data to be averaged
     * @param average_type    averaging type to be used ("gaussian" or "tophat")
     */
    template <typename MapType>
    void average(MapType& exp, String average_type)
    {
      // MS level to be averaged
      int ms_level = param_.getValue("average_gaussian:ms_level");
      if (average_type == "tophat")
      {
        ms_level = param_.getValue("average_tophat:ms_level");
      }
      
      // spectrum type (profile, centroid or automatic)
      String spectrum_type = param_.getValue("average_gaussian:spectrum_type");
      if (average_type == "tophat")
      {
        spectrum_type = param_.getValue("average_tophat:spectrum_type");
      }

      // parameters for Gaussian averaging
      double fwhm(param_.getValue("average_gaussian:rt_FWHM"));
      double factor = -4 * log(2.0) / (fwhm * fwhm); // numerical factor within Gaussian
      double cutoff(param_.getValue("average_gaussian:cutoff"));

      // parameters for Top-Hat averaging
      bool unit(param_.getValue("average_tophat:rt_unit") == "scans"); // true if RT unit is 'scans', false if RT unit is 'seconds'
      double range(param_.getValue("average_tophat:rt_range")); // range of spectra to be averaged over
      double range_seconds = range / 2; // max. +/- <range_seconds> seconds from master spectrum
      int range_scans = range;
      if ((range_scans % 2) == 0)
      {
        ++range_scans;
      }
      range_scans = (range_scans - 1) / 2; // max. +/- <range_scans> scans from master spectrum

      AverageBlocks spectra_to_average_over;

      // loop over RT
      int n(0); // spectrum index
      for (typename MapType::const_iterator it_rt = exp.begin(); it_rt != exp.end(); ++it_rt)
      {
        if (Int(it_rt->getMSLevel()) == ms_level)
        {
          int m; // spectrum index
          int steps;
          bool terminate_now;
          typename MapType::const_iterator it_rt_2;

          // go forward (start at next downstream spectrum; the current spectrum will be covered when looking backwards)
          steps = 0;
          m = n + 1;
          it_rt_2 = it_rt + 1;
          terminate_now = false;
          while (it_rt_2 != exp.end() && !terminate_now)
          {
            if (Int(it_rt_2->getMSLevel()) == ms_level)
            {
              double weight = 1;
              if (average_type == "gaussian")
              {
                weight = std::exp(factor * pow(it_rt_2->getRT() - it_rt->getRT(), 2));
              }
              std::pair<Size, double> p(m, weight);
              spectra_to_average_over[n].push_back(p);
              ++steps;
            }
            if (average_type == "gaussian")
            {
              // Gaussian
              terminate_now = std::exp(factor * pow(it_rt_2->getRT() - it_rt->getRT(), 2)) < cutoff;
            }
            else if (unit)
            {
              // Top-Hat with RT unit = scans
              terminate_now = (steps > range_scans);
            }
            else
            {
              // Top-Hat with RT unit = seconds
              terminate_now = (std::abs(it_rt_2->getRT() - it_rt->getRT()) > range_seconds);
            }
            ++m;
            ++it_rt_2;
          }

          // go backward
          steps = 0;
          m = n;
          it_rt_2 = it_rt;
          terminate_now = false;
          while (it_rt_2 != exp.begin() && !terminate_now)
          {
            if (Int(it_rt_2->getMSLevel()) == ms_level)
            {
              double weight = 1;
              if (average_type == "gaussian")
              {
                weight = std::exp(factor * pow(it_rt_2->getRT() - it_rt->getRT(), 2));
              }
              std::pair<Size, double> p(m, weight);
              spectra_to_average_over[n].push_back(p);
              ++steps;
            }
            if (average_type == "gaussian")
            {
              // Gaussian
              terminate_now = std::exp(factor * pow(it_rt_2->getRT() - it_rt->getRT(), 2)) < cutoff;
            }
            else if (unit)
            {
              // Top-Hat with RT unit = scans
              terminate_now = (steps > range_scans);
            }
            else
            {
              // Top-Hat with RT unit = seconds
              terminate_now = (std::abs(it_rt_2->getRT() - it_rt->getRT()) > range_seconds);
            }
            --m;
            --it_rt_2;
          }

        }
        ++n;
      }

      // normalize weights
      for (AverageBlocks::Iterator it = spectra_to_average_over.begin(); it != spectra_to_average_over.end(); ++it)
      {
        double sum(0.0);
        for (std::vector<std::pair<Size, double> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          sum += it2->second;
        }

        for (std::vector<std::pair<Size, double> >::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          (*it2).second /= sum;
        }
      }

      // determine type of spectral data (profile or centroided)
      SpectrumSettings::SpectrumType type;
      if (spectrum_type == "automatic")
      {
        Size idx = spectra_to_average_over.begin()->first; // index of first spectrum to be averaged
        type = exp[idx].getType();
        if (type == SpectrumSettings::UNKNOWN)
        {
          type = PeakTypeEstimator().estimateType(exp[idx].begin(), exp[idx].end());
        }
      }
      else if (spectrum_type == "profile")
      {
        type = SpectrumSettings::PROFILE;
      }
      else if (spectrum_type == "centroid")
      {
        type = SpectrumSettings::CENTROID;
      }

      // generate new spectra
      if (type == SpectrumSettings::CENTROID)
      {
        averageCentroidSpectra_(exp, spectra_to_average_over, ms_level);
      }
      else
      {
        averageProfileSpectra_(exp, spectra_to_average_over, ms_level);
      }

      exp.sortSpectra();

      return;
    }

    // @}

protected:

    /**
        @brief merges blocks of spectra of a certain level

        Merges spectra belonging to the same block, setting their MS level to @p ms_level.
        All old spectra of level @p ms_level are removed, and the new consensus spectra (one per block)
        are added.
        All spectra with other MS levels remain untouched.
        The resulting map is NOT sorted!

    */
    template <typename MapType>
    void mergeSpectra_(MapType& exp, const MergeBlocks& spectra_to_merge, const UInt ms_level)
    {
      double mz_binning_width(param_.getValue("mz_binning_width"));
      String mz_binning_unit(param_.getValue("mz_binning_width_unit"));

      // merge spectra
      MapType merged_spectra;

      Map<Size, Size> cluster_sizes;
      std::set<Size> merged_indices;

      // set up alignment
      SpectrumAlignment sas;
      Param p;
      p.setValue("tolerance", mz_binning_width);
      if (!(mz_binning_unit == "Da" || mz_binning_unit == "ppm"))
      {
        throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);  // sanity check
      }

      p.setValue("is_relative_tolerance", mz_binning_unit == "Da" ? "false" : "true");
      sas.setParameters(p);
      std::vector<std::pair<Size, Size> > alignment;

      Size count_peaks_aligned(0);
      Size count_peaks_overall(0);

      // each BLOCK
      for (auto it = spectra_to_merge.begin(); it != spectra_to_merge.end(); ++it)
      {
        ++cluster_sizes[it->second.size() + 1]; // for stats

        typename MapType::SpectrumType consensus_spec = exp[it->first];
        consensus_spec.setMSLevel(ms_level);

        //consensus_spec.unify(exp[it->first]); // append meta info
        merged_indices.insert(it->first);

        //typename MapType::SpectrumType all_peaks = exp[it->first];
        double rt_average = consensus_spec.getRT();
        double precursor_mz_average = 0.0;
        Size precursor_count(0);
        if (!consensus_spec.getPrecursors().empty())
        {
          precursor_mz_average = consensus_spec.getPrecursors()[0].getMZ();
          ++precursor_count;
        }

        count_peaks_overall += consensus_spec.size();

        // block elements
        for (auto sit = it->second.begin(); sit != it->second.end(); ++sit)
        {
          consensus_spec.unify(exp[*sit]); // append meta info
          merged_indices.insert(*sit);

          rt_average += exp[*sit].getRT();
          if (ms_level >= 2 && exp[*sit].getPrecursors().size() > 0)
          {
            precursor_mz_average += exp[*sit].getPrecursors()[0].getMZ();
            ++precursor_count;
          }

          // merge data points
          sas.getSpectrumAlignment(alignment, consensus_spec, exp[*sit]);
          //std::cerr << "alignment of " << it->first << " with " << *sit << " yielded " << alignment.size() << " common peaks!\n";
          count_peaks_aligned += alignment.size();
          count_peaks_overall += exp[*sit].size();

          Size align_index(0);
          Size spec_b_index(0);

          // sanity check for number of peaks
          Size spec_a = consensus_spec.size(), spec_b = exp[*sit].size(), align_size = alignment.size();
          for (auto pit = exp[*sit].begin(); pit != exp[*sit].end(); ++pit)
          {
            if (alignment.size() == 0 || alignment[align_index].second != spec_b_index)
              // ... add unaligned peak
            {
              consensus_spec.push_back(*pit);
            }
            // or add aligned peak height to ALL corresponding existing peaks
            else
            {
              Size counter(0);
              Size copy_of_align_index(align_index);

              while (alignment.size() > 0 && 
                     copy_of_align_index < alignment.size() && 
                     alignment[copy_of_align_index].second == spec_b_index)
              {
                ++copy_of_align_index;
                ++counter;
              } // Count the number of peaks in a which correspond to a single b peak.

              while (alignment.size() > 0 &&
                     align_index < alignment.size() &&  
                     alignment[align_index].second == spec_b_index)
              {
                consensus_spec[alignment[align_index].first].setIntensity(consensus_spec[alignment[align_index].first].getIntensity() +
                    (pit->getIntensity() / (double)counter)); // add the intensity divided by the number of peaks
                ++align_index; // this aligned peak was explained, wait for next aligned peak ...
                if (align_index == alignment.size())
                {
                  alignment.clear();  // end reached -> avoid going into this block again
                }
              }
              align_size = align_size + 1 - counter; //Decrease align_size by number of
            }
            ++spec_b_index;
          }
          consensus_spec.sortByPosition(); // sort, otherwise next alignment will fail
          if (spec_a + spec_b - align_size != consensus_spec.size())
          {
            LOG_WARN << "wrong number of features after merge. Expected: " << spec_a + spec_b - align_size << " got: " << consensus_spec.size() << "\n";
          }
        }
        rt_average /= it->second.size() + 1;
        consensus_spec.setRT(rt_average);

        if (ms_level >= 2)
        {
          if (precursor_count)
          {
            precursor_mz_average /= precursor_count;
          }
          std::vector<Precursor> pcs = consensus_spec.getPrecursors();
          //if (pcs.size()>1) LOG_WARN << "Removing excessive precursors - leaving only one per MS2 spectrum.\n";
          pcs.resize(1);
          pcs[0].setMZ(precursor_mz_average);
          consensus_spec.setPrecursors(pcs);
        }

        if (consensus_spec.empty())
        {
          continue;
        }
        else
        {
          merged_spectra.addSpectrum(consensus_spec);
        }
      }

      LOG_INFO << "Cluster sizes:\n";
      for (Map<Size, Size>::const_iterator it = cluster_sizes.begin(); it != cluster_sizes.end(); ++it)
      {
        LOG_INFO << "  size " << it->first << ": " << it->second << "x\n";
      }

      char buffer[200];
      sprintf(buffer, "%d/%d (%.2f %%) of blocked spectra", (int)count_peaks_aligned,
              (int)count_peaks_overall, float(count_peaks_aligned) / float(count_peaks_overall) * 100.);
      LOG_INFO << "Number of merged peaks: " << String(buffer) << "\n";

      // remove all spectra that were within a cluster
      typename MapType::SpectrumType empty_spec;
      MapType exp_tmp;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (merged_indices.count(i) == 0) // save unclustered ones
        {
          exp_tmp.addSpectrum(exp[i]);
          exp[i] = empty_spec;
        }
      }

      //typedef std::vector<typename MapType::SpectrumType> Base;
      //exp.Base::operator=(exp_tmp);
      exp.clear(false);
      exp.getSpectra().insert(exp.end(), exp_tmp.begin(), exp_tmp.end());

      // exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<typename MapType::SpectrumType>(ListUtils::create<int>(String(ms_level)), false)), exp.end());

      // ... and add consensus spectra
      exp.getSpectra().insert(exp.end(), merged_spectra.begin(), merged_spectra.end());

    }

    /**
     * @brief average spectra (profile mode)
     *
     * Averages spectra in profile mode of one MS level in an experiment. The
     * blocks of spectra to be combined and their relative weights have
     * previously been determined. The averaged spectra are generated in two
     * steps:
     * (1) The m/z of all spectra in a block are collected and sorted. m/z
     *     positions closer than mz_binning_width are removed.
     * (2) At these positions the weighted sum of all spline interpolations is
     *     calculated.
     *
     * The first step ensures roughly the same sampling rate as the one of the
     * original spectra. The exact m/z position is not crucial, since not the
     * original intensities but the spline-interpolated intensities are used.
     *
     * @param exp   experimental data to be averaged
     * @param spectra_to_average_over    mapping of spectral index to set of spectra to average over with corresponding weights
     * @param ms_level    MS level of spectra to be averaged
     */
    template <typename MapType>
    void averageProfileSpectra_(MapType& exp, const AverageBlocks& spectra_to_average_over, const UInt ms_level)
    {
      MapType exp_tmp; // temporary experiment for averaged spectra

      double mz_binning_width(param_.getValue("mz_binning_width"));
      String mz_binning_unit(param_.getValue("mz_binning_width_unit"));

      unsigned progress = 0;
      std::stringstream progress_message;
      progress_message << "averaging profile spectra of MS level " << ms_level;
      startProgress(0, spectra_to_average_over.size(), progress_message.str());

      // loop over blocks
      for (AverageBlocks::ConstIterator it = spectra_to_average_over.begin(); it != spectra_to_average_over.end(); ++it)
      {
        setProgress(++progress);

        // loop over spectra in blocks
        std::vector<double> mz_positions_all; // m/z positions from all spectra
        for (std::vector<std::pair<Size, double> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          // loop over m/z positions
          for (typename MapType::SpectrumType::ConstIterator it_mz = exp[it2->first].begin(); it_mz < exp[it2->first].end(); ++it_mz)
          {
            mz_positions_all.push_back(it_mz->getMZ());
          }
        }

        sort(mz_positions_all.begin(), mz_positions_all.end());

        std::vector<double> mz_positions; // positions at which the averaged spectrum should be evaluated
        std::vector<double> intensities;
        double last_mz = std::numeric_limits<double>::min(); // last m/z position pushed through from mz_position to mz_position_2
        double delta_mz(mz_binning_width); // for m/z unit Da
        for (std::vector<double>::iterator it_mz = mz_positions_all.begin(); it_mz < mz_positions_all.end(); ++it_mz)
        {
          if (mz_binning_unit == "ppm")
          {
            delta_mz = mz_binning_width * (*it_mz) / 1000000;
          }

          if (((*it_mz) - last_mz) > delta_mz)
          {
            mz_positions.push_back(*it_mz);
            intensities.push_back(0.0);
            last_mz = *it_mz;
          }
        }

        // loop over spectra in blocks
        for (std::vector<std::pair<Size, double> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          SplineSpectrum spline(exp[it2->first]);
          SplineSpectrum::Navigator nav = spline.getNavigator();

          // loop over m/z positions
          for (Size i = 0; i < mz_positions.size(); ++i)
          {
            if ((spline.getMzMin() < mz_positions[i]) && (mz_positions[i] < spline.getMzMax()))
            {
              intensities[i] += nav.eval(mz_positions[i]) * (it2->second); // spline-interpolated intensity * weight
            }
          }
        }

        // update spectrum
        typename MapType::SpectrumType average_spec = exp[it->first];
        average_spec.clear(false); // Precursors are part of the meta data, which are not deleted.
        //average_spec.setMSLevel(ms_level);

        // refill spectrum
        for (Size i = 0; i < mz_positions.size(); ++i)
        {
          typename MapType::PeakType peak;
          peak.setMZ(mz_positions[i]);
          peak.setIntensity(intensities[i]);
          average_spec.push_back(peak);
        }

        // store spectrum temporarily
        exp_tmp.addSpectrum(average_spec);
      }

      endProgress();

      // loop over blocks
      int n(0);
      //typename MapType::SpectrumType empty_spec;
      for (AverageBlocks::ConstIterator it = spectra_to_average_over.begin(); it != spectra_to_average_over.end(); ++it)
      {
        exp[it->first] = exp_tmp[n];
        //exp_tmp[n] = empty_spec;
        ++n;
      }

    }

    /**
     * @brief average spectra (centroid mode)
     *
     * Averages spectra in centroid mode of one MS level in an experiment. The
     * blocks of spectra to be combined and their relative weights have
     * previously determined. The averaged spectra are generated in two steps:
     * (1) The m/z of all spectra in a block are collected and sorted. Their
     *     corresponding intensities are weighted.
     * (2) m/z positions closer than mz_binning_width are combined to a single
     *     peak. The m/z are averaged and the corresponding intensities summed.
     *
     * @param exp   experimental data to be averaged
     * @param spectra_to_average_over    mapping of spectral index to set of spectra to average over with corresponding weights
     * @param ms_level    MS level of spectra to be averaged
     */
    template <typename MapType>
    void averageCentroidSpectra_(MapType& exp, const AverageBlocks& spectra_to_average_over, const UInt ms_level)
    {
      MapType exp_tmp; // temporary experiment for averaged spectra

      double mz_binning_width(param_.getValue("mz_binning_width"));
      String mz_binning_unit(param_.getValue("mz_binning_width_unit"));

      unsigned progress = 0;
      ProgressLogger logger;
      std::stringstream progress_message;
      progress_message << "averaging centroid spectra of MS level " << ms_level;
      logger.startProgress(0, spectra_to_average_over.size(), progress_message.str());

      // loop over blocks
      for (AverageBlocks::ConstIterator it = spectra_to_average_over.begin(); it != spectra_to_average_over.end(); ++it)
      {
        logger.setProgress(++progress);

        // collect peaks from all spectra
        // loop over spectra in blocks
        std::vector<std::pair<double, double> > mz_intensity_all; // m/z positions and peak intensities from all spectra
        for (std::vector<std::pair<Size, double> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
          // loop over m/z positions
          for (typename MapType::SpectrumType::ConstIterator it_mz = exp[it2->first].begin(); it_mz < exp[it2->first].end(); ++it_mz)
          {
            std::pair<double, double> mz_intensity(it_mz->getMZ(), (it_mz->getIntensity() * it2->second)); // m/z, intensity * weight
            mz_intensity_all.push_back(mz_intensity);
          }
        }

        sort(mz_intensity_all.begin(), mz_intensity_all.end(), SpectraMerger::compareByFirst);

        // generate new spectrum
        std::vector<double> mz_new;
        std::vector<double> intensity_new;
        double last_mz = std::numeric_limits<double>::min();
        double delta_mz = mz_binning_width;
        double sum_mz(0);
        double sum_intensity(0);
        Size count(0);
        for (std::vector<std::pair<double, double> >::const_iterator it_mz = mz_intensity_all.begin(); it_mz != mz_intensity_all.end(); ++it_mz)
        {
          if (mz_binning_unit == "ppm")
          {
            delta_mz = mz_binning_width * (it_mz->first) / 1000000;
          }

          if (((it_mz->first - last_mz) > delta_mz) && (count > 0))
          {
            mz_new.push_back(sum_mz / count);
            intensity_new.push_back(sum_intensity); // intensities already weighted

            sum_mz = 0;
            sum_intensity = 0;

            last_mz = it_mz->first;
            count = 0;
          }

          sum_mz += it_mz->first;
          sum_intensity += it_mz->second;
          ++count;
        }
        if (count > 0)
        {
          mz_new.push_back(sum_mz / count);
          intensity_new.push_back(sum_intensity); // intensities already weighted
        }

        // update spectrum
        typename MapType::SpectrumType average_spec = exp[it->first];
        average_spec.clear(false); // Precursors are part of the meta data, which are not deleted.
        //average_spec.setMSLevel(ms_level);

        // refill spectrum
        for (Size i = 0; i < mz_new.size(); ++i)
        {
          typename MapType::PeakType peak;
          peak.setMZ(mz_new[i]);
          peak.setIntensity(intensity_new[i]);
          average_spec.push_back(peak);
        }

        // store spectrum temporarily
        exp_tmp.addSpectrum(average_spec);

      }

      logger.endProgress();

      // loop over blocks
      int n(0);
      for (AverageBlocks::ConstIterator it = spectra_to_average_over.begin(); it != spectra_to_average_over.end(); ++it)
      {
        exp[it->first] = exp_tmp[n];
        ++n;
      }

    }

    /**
     * @brief comparator for sorting peaks (m/z, intensity)
     */
    bool static compareByFirst(std::pair<double, double> i, std::pair<double, double> j)
    {
      return i.first < j.first;
    }

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_SPECTRAMERGER_H
