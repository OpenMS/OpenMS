// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/SILACAnalyzer.h>

using namespace std;

namespace OpenMS 
{

  void SILACAnalyzer::filterData(MSExperiment<Peak1D> & exp, const PeakWidthEstimator::Result & peak_width, vector<vector<SILACPattern> > & data)
  {
    list<SILACFilter> filters;

    // create filters for all numbers of isotopes per peptide, charge states and mass shifts
    // iterate over all number for peaks per peptide (from max to min)
    for (UInt isotopes_per_peptide = isotopes_per_peptide_max; isotopes_per_peptide >= isotopes_per_peptide_min; isotopes_per_peptide--)
    {
      // iterate over all charge states (from max to min)
      for (UInt charge = charge_max; charge >= charge_min; charge--)
      {
        // iterate over all mass shifts
        for (UInt i = 0; i < massShifts.size(); i++)
        {
          // convert std::vector<DoubleReal> to set<DoubleReal> for SILACFilter
          std::vector<DoubleReal> massShifts_set = massShifts[i];

          //copy(massShifts[i].begin(), massShifts[i].end(), inserter(massShifts_set, massShifts_set.end()));
          filters.push_back(SILACFilter(massShifts_set, charge, model_deviation, isotopes_per_peptide, intensity_cutoff, intensity_correlation, allow_missing_peaks));
        }
      }
    }

    // create filtering
    SILACFiltering filtering(exp, peak_width, intensity_cutoff, out_debug);
    filtering.setLogType(getLogType());

    // register filters to the filtering
    for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
    {
      filtering.addFilter(*filter_it);
    }

    // perform filtering
    filtering.filterDataPoints();

    // retrieve filtered data points
    for (SILACFiltering::Filters::iterator filter_it = filtering.filters_.begin(); filter_it != filtering.filters_.end(); ++filter_it)
    {
      data.push_back(filter_it->getElements());
    }


    //--------------------------------------------------
    // combine DataPoints to improve the clustering
    //--------------------------------------------------

    // DataPoints that originate from filters with same charge state and mass shift(s)
    // and whose filters only differ in number of isotopes per peptide are combined
    // to get one cluster for peptides whose elution profile varies in number of isotopes per peptide

    // perform combination only if the user specified a peaks_per_peptide range > 1
    if (isotopes_per_peptide_min != isotopes_per_peptide_max)
    {
      // erase empty filter results from "data"
      std::vector<std::vector<SILACPattern> > data_temp;

      for (std::vector<std::vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        if (data_it->size() != 0)
        {
          data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
        }
      }

      data.swap(data_temp);     // data = data_temp
      data_temp.clear();      // clear "data_temp"

      if (data.size() >= 2)
      {
        Int temp = 0;
        // combine corresponding DataPoints
        std::vector<std::vector<SILACPattern> >::iterator data_it_1 = data.begin();      // first iterator over "data" to get first DataPoint for combining
        std::vector<std::vector<SILACPattern> >::iterator data_it_2 = data_it_1 + 1;     // second iterator over "data" to get second DataPoint for combining
        std::vector<std::vector<SILACPattern> >::iterator data_it_end = data.end() - 1;      // pointer to second last elemnt of "data"
        std::vector<SILACPattern>::iterator it_1;     // first inner iterator over elements of first DataPoint
        std::vector<SILACPattern>::iterator it_2;     // second inner iterator over elements of second DataPoint

        while (data_it_1 < data_it_end)      // check for combining as long as first DataPoint is not second last elment of "data"
        {
          while (data_it_1->size() == 0 && data_it_1 < data_it_end)
          {
            ++data_it_1;      // get next first DataPoint
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          if (data_it_1 == data_it_end && data_it_2 == data.end())     // if first iterator points to last element of "data" and second iterator points to end of "data"
          {
            break;      // stop combining
          }

          while (data_it_2 < data.end() && data_it_2->size() == 0)      // as long as current second DataPoint is empty and second iterator does not point to end of "data"
          {
            ++data_it_2;      // get next second DataPoint
          }

          if (data_it_2 == data.end())      // if second iterator points to end of "data"
          {
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          it_1 = data_it_1->begin();      // set first inner iterator to first element of first DataPoint
          it_2 = data_it_2->begin();      // set second inner iterator to first element of second DataPoint

          // check if DataPoints are not empty
          if (data_it_1->size() != 0 && data_it_2->size() != 0)
          {
            // check if DataPoints have the same charge state and mass shifts
            if (it_1->charge != it_2->charge || it_1->mass_shifts != it_2->mass_shifts)
            {
              if (data_it_2 < data_it_end)     // if DataPpoints differ and second DataPoint is not second last element of "data"
              {
                temp++;
                ++data_it_2;      // get next second DataPoint
                if (temp > 50000)
                {
                  ++data_it_1;
                  temp = 0;
                }
              }
              else if (data_it_2 == data_it_end && data_it_1 < data.end() - 2)     // if DataPpoints differ and second DataPoint is second last element of "data" and first DataPoint is not third last element of "data"
              {
                ++data_it_1;      // get next first DataPoint
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }
              else
              {
                ++data_it_1;      // get next first DataPoint
              }
            }
            else
            {
              // perform combining
              (*data_it_1).insert(data_it_1->end(), data_it_2->begin(), data_it_2->end());      // append second DataPoint to first DataPoint
              (*data_it_2).clear();     // clear second Datapoint to keep iterators valid and to keep size of "data"

              if (data_it_2 < data_it_end)     // if second DataPoint is not second last element of "data"
              {
                ++data_it_2;      // get next second DataPoint
              }
              else
              {
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }
            }
          }
          else
          {
            ++data_it_1;      // get next first DataPoint
          }
        }

        // erase empty DataPoints from "data"
        std::vector<std::vector<SILACPattern> > data_temp;

        for (std::vector<std::vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
        {
          if (data_it->size() != 0)
          {
            data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
          }
        }

        data.swap(data_temp);     // data = data_temp
        data_temp.clear();      // clear "data_temp"
      }
    }


  }

  void SILACAnalyzer::clusterData(const MSExperiment<> & exp, const PeakWidthEstimator::Result & peak_width, vector<Clustering *> & cluster_data, vector<vector<SILACPattern> > & data)
  {
    typedef Clustering::PointCoordinate PointCoordinate;

    ProgressLogger progresslogger;
    progresslogger.setLogType(getLogType());
    progresslogger.startProgress(0, data.size(), "clustering data");

    // Use peak half width @1000 Th for mz threshold
    DoubleReal mz_threshold = peak_width(1000);

    // Use double median of spectrum spacing for max rt spacing
    DoubleReal rt_max_spacing = 0;
    {
      // Calculate distance between each spectrum; this needs sorted spectra
      std::vector<Real> space;
      MSExperiment<>::const_iterator it1 = exp.begin();
      MSExperiment<>::const_iterator it2 = exp.begin(); ++it2;
      for (; it2 != exp.end(); ++it1, ++it2)
      {
        DoubleReal s = it2->getRT() - it1->getRT();
        space.push_back(s);
      }

      sort(space.begin(), space.end());

      // Calculate median by extracting the middle element (okay, the upper median)
      // Set max spacing to five times the median spectrum spacing
      // The five is an empirical value
      if (space.size())
        rt_max_spacing = space[space.size() / 2 + 1] * 5;
    }

    UInt data_id = 0;

    for (vector<vector<SILACPattern> >::iterator data_it = data.begin();
         data_it != data.end();
         ++data_it, ++data_id)
    {
      const PointCoordinate max_delta(rt_threshold, mz_threshold);
      Clustering * clustering = new Clustering(max_delta, rt_min, rt_max_spacing);

      for (vector<SILACPattern>::iterator it = data_it->begin(); it != data_it->end(); ++it)
      {
        const PointCoordinate key(it->rt, it->mz);
        SILACPattern & p = *it;
        clustering->insertPoint(key, &p);
      }

      clustering->cluster();

      cluster_data.push_back(clustering);

      progresslogger.setProgress(data_id);
    }

    progresslogger.endProgress();
  }

  PeakWidthEstimator::Result SILACAnalyzer::estimatePeakWidth(const MSExperiment<Peak1D> & exp)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(getLogType());
    progresslogger.startProgress(0, 1, "estimate peak width");

    PeakWidthEstimator::Result ret = PeakWidthEstimator::estimateFWHM(exp);

    progresslogger.endProgress();
    std::cout << "Estimated peak width: e ^ (" << ret.c0 << " + " << ret.c1 << " * log mz)" << std::endl;
    return ret;
  }

  void SILACAnalyzer::generateClusterConsensusByCluster(ConsensusMap & out, const Clustering & clustering) const
  {
    // iterate over clusters
    for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it)
    {
      ConsensusFeature consensus;

      // determine the number of peptides
      // for that we look at the first point in the first pattern
      const SILACPattern & firstPattern = *(cluster_it->second.begin())->second;
      const SILACPoint & firstPoint = *(firstPattern.points.begin());
      UInt numberPeptides = firstPoint.intensities.size();
      UInt charge = firstPoint.charge;

      // sums for each peptide of the pair (triplet, singlet, ...)
      std::vector<DoubleReal> sumMzIntensities(numberPeptides, 0);    // sum m/z * intensity (for intensity-weighted m/z average)
      std::vector<DoubleReal> sumRtIntensities(numberPeptides, 0);    // sum rt * intensity (for intensity-weighted rt average)
      std::vector<DoubleReal> sumIntensities(numberPeptides, 0);    // sum intensity (for 'feature volume' = peptide intensity)
      std::vector<DoubleReal> sumIntensitiesMonoisotopic(numberPeptides, 0);    // sum intensity of monoisotopic mass trace (for normalisation of intensity-weighted m/z average)
      std::vector<DoubleReal> maxIntensityXIC(numberPeptides, 0);    // tracks maximum of sumIntensitiesXIC
      std::vector<DoubleReal> RtAtMaxIntensityXIC(numberPeptides, 0);    // tracks rt at maximum of sumIntensitiesXIC

      // iterate over SILAC patterns in each cluster
      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
           pattern_it != cluster_it->second.end();
           ++pattern_it)
      {
        const SILACPattern & pattern = *pattern_it->second;
        std::vector<DoubleReal> sumIntensitiesXIC(numberPeptides, 0);    // sums intensities at fixed rt (for XIC)

        // iterate over SILAC points in each SILAC pattern
        for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
             point_it != pattern.points.end();
             ++point_it)
        {
          const SILACPoint & point = *point_it;

          // iterate over peptides in doublet (or triplet, ...)
          UInt peptide = 0;    // peptide for which intensity, retention time and m/z are to be calculated
          for (std::vector<std::vector<DoubleReal> >::const_iterator peptide_it = point.intensities.begin();
               peptide_it != point.intensities.end();
               ++peptide_it)
          {
            // iterate over isotopes in peptide
            UInt isotope = 0;
            for (std::vector<DoubleReal>::const_iterator isotope_it = peptide_it->begin();
                 isotope_it != peptide_it->end();
                 ++isotope_it)
            {
              sumIntensities[peptide] += point.intensities[peptide][isotope];
              sumIntensitiesXIC[peptide] += point.intensities[peptide][isotope];
              sumRtIntensities[peptide] += point.rt * point.intensities[peptide][isotope];
              if (isotope == 0)
              {
                sumMzIntensities[peptide] += pattern.mz_positions[peptide][isotope] * point.intensities[peptide][isotope];
                sumIntensitiesMonoisotopic[peptide] += point.intensities[peptide][isotope];
              }
              //sumMzIntensities[peptide] += (pattern.mz_positions[peptide][isotope] - (isotope * 1.003355 / point.charge)) * point.intensities[peptide][isotope];
              ++isotope;
            }
            ++peptide;
          }

        }

        // check for each peptide if its XIC intensity has been raised
        for (UInt peptide = 0; peptide < numberPeptides; ++peptide)
        {
          if (sumIntensitiesXIC[peptide] > maxIntensityXIC[peptide])
          {
            maxIntensityXIC[peptide] = sumIntensitiesXIC[peptide];
            RtAtMaxIntensityXIC[peptide] = pattern.rt;
          }
        }
      }
      /*cout << "light m/z: " << sumMzIntensities[0]/sumIntensitiesMonoisotopic[0] << '\n';
       cout << "light rt (intensity averaged): " << sumRtIntensities[0]/sumIntensities[0] << '\n';
       cout << "light rt (at max XIC): " << RtAtMaxIntensityXIC[0] << '\n';
       cout << "heavy m/z: " << sumMzIntensities[1]/sumIntensitiesMonoisotopic[1] << '\n';
       cout << "heavy rt (intensity averaged): " << sumRtIntensities[1]/sumIntensities[1] << '\n';
       cout << "heavy rt (at max XIC): " << RtAtMaxIntensityXIC[1] << '\n' << '\n';*/

      // consensus feature has coordinates of the light peptide
      consensus.setMZ(sumMzIntensities[0] / sumIntensitiesMonoisotopic[0]);    // intensity-average only over the mono-isotopic peak
      consensus.setRT(sumRtIntensities[0] / sumIntensities[0]);    // intensity-average over the entire peptide, i.e. all peptides
      //consensus.setRT(RtAtMaxIntensityXIC[0]);
      consensus.setIntensity(sumIntensities[0]);
      consensus.setCharge(charge);
      consensus.setQuality(std::floor(firstPattern.mass_shifts[1] * charge));    // set Quality to the first mass shift (allows later to filter in consensXML)

      // attach features to consensus
      for (UInt peptide = 0; peptide < numberPeptides; ++peptide)
      {
        FeatureHandle feature;

        feature.setMZ(sumMzIntensities[peptide] / sumIntensitiesMonoisotopic[peptide]);
        feature.setRT(sumRtIntensities[peptide] / sumIntensities[peptide]);
        //feature.setRT(RtAtMaxIntensityXIC[peptide]);
        feature.setIntensity(sumIntensities[peptide]);
        feature.setCharge(charge);
        feature.setMapIndex(peptide);
        out.getFileDescriptions()[peptide].size++;

        consensus.insert(feature);
      }

      // add consensus to consensus map
      out.push_back(consensus);
    }
  }

  void SILACAnalyzer::generateClusterConsensusByPattern(ConsensusMap & out, const Clustering & clustering, UInt & cluster_id) const
  {
    for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it, ++cluster_id)
    {
      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin(); pattern_it != cluster_it->second.end(); ++pattern_it)
      {
        ConsensusFeature consensus = generateSingleConsensusByPattern(*pattern_it->second);

        consensus.setMetaValue("color", selectColor(cluster_id));
        consensus.setMetaValue("Cluster ID", cluster_id);

        out.getFileDescriptions()[0].size++;

        out.push_back(consensus);
      }
    }
  }

  void SILACAnalyzer::generateClusterDebug(std::ostream & out, const Clustering & clustering, UInt & cluster_id) const
  {
    for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin();
         cluster_it != clustering.grid.end();
         ++cluster_it, ++cluster_id)
    {
      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
           pattern_it != cluster_it->second.end();
           ++pattern_it)
      {
        const SILACPattern & pattern = *pattern_it->second;

        std::ostringstream preamble;

        preamble
        << std::fixed << std::setprecision(8)
        << cluster_id << ','
        << pattern.rt << ','
        << pattern.mz << ','
        << pattern.charge << ',';

        for (std::vector<DoubleReal>::const_iterator shift_it = pattern.mass_shifts.begin();
             shift_it != pattern.mass_shifts.end();
             ++shift_it)
        {
          preamble
          << *++shift_it * pattern.charge << ',';
        }

        for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = pattern.intensities.begin();
             shift_inten_it != pattern.intensities.end();
             ++shift_inten_it)
        {
          UInt peak_inten_id = 0;
          for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
               peak_inten_it != shift_inten_it->end();
               ++peak_inten_it, ++peak_inten_id)
          {
            preamble
            << *peak_inten_it << ',';
          }
          for (; peak_inten_id < isotopes_per_peptide_max; ++peak_inten_id)
          {
            preamble
            << "NA,";
          }
        }

        for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
             point_it != pattern.points.end();
             ++point_it)
        {
          const SILACPoint & point = *point_it;

          out
          << preamble.str()
          << point.mz;

          // write INT_RAW_...
          for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = point.intensities.begin();
               shift_inten_it != point.intensities.end();
               ++shift_inten_it)
          {
            UInt peak_inten_id = 0;
            for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
                 peak_inten_it != shift_inten_it->end();
                 ++peak_inten_it, ++peak_inten_id)
            {
              out << ','
                  << *peak_inten_it;
            }
            for (; peak_inten_id < isotopes_per_peptide_max; ++peak_inten_id)
            {
              out
              << ",NA";
            }
          }

          // write MZ_RAW_...
          for (std::vector<std::vector<DoubleReal> >::const_iterator shift_mz_it = point.mz_positions.begin();
               shift_mz_it != point.mz_positions.end();
               ++shift_mz_it)
          {
            UInt peak_mz_id = 0;
            for (std::vector<DoubleReal>::const_iterator peak_mz_it = shift_mz_it->begin();
                 peak_mz_it != shift_mz_it->end();
                 ++peak_mz_it, ++peak_mz_id)
            {
              out << ','
                  << *peak_mz_it;
            }
            for (; peak_mz_id < isotopes_per_peptide_max; ++peak_mz_id)
            {
              out
              << ",NA";
            }
          }

          out << '\n';
        }
      }
    }
  }

  void SILACAnalyzer::generateFilterConsensusByPattern(ConsensusMap & out, const std::vector<SILACPattern> & pattern) const
  {
    for (std::vector<SILACPattern>::const_iterator pattern_it = pattern.begin(); pattern_it != pattern.end(); ++pattern_it)
    {
      out.push_back(generateSingleConsensusByPattern(*pattern_it));
    }
  }

  ConsensusFeature SILACAnalyzer::generateSingleConsensusByPattern(const SILACPattern & pattern) const
  {
    // XXX: get from experiment
    Int charge = pattern.charge;

    ConsensusFeature consensus;
    consensus.setRT(pattern.rt);
    consensus.setMZ(pattern.mz);
    consensus.setIntensity(pattern.intensities[0][0]);
    consensus.setCharge(charge);

    consensus.setMetaValue("Peaks per peptide", pattern.isotopes_per_peptide);

    // Output mass shifts
    {
      std::ostringstream out;
      out << std::fixed << std::setprecision(4);
      for (vector<DoubleReal>::const_iterator shift_it = pattern.mass_shifts.begin() + 1; shift_it != pattern.mass_shifts.end(); ++shift_it)
      {
        out << *shift_it * charge << ';';
      }
      // Remove the last delimiter
      std::string outs = out.str(); outs.erase(outs.end() - 1);
      consensus.setQuality(std::floor(pattern.mass_shifts.at(1) * charge));
      consensus.setMetaValue("SILAC", outs);
    }

    // Output all intensities per peptide as list
    {
      std::ostringstream out;
      for (vector<vector<DoubleReal> >::const_iterator inten_it = pattern.intensities.begin(); inten_it != pattern.intensities.end(); ++inten_it)
      {
        std::ostringstream out2;
        out2 << std::fixed << std::setprecision(4);
        for (vector<DoubleReal>::const_iterator inten2_it = inten_it->begin(); inten2_it != inten_it->end(); ++inten2_it)
        {
          out2 << *inten2_it << ',';
        }
        // Remove the last delimiter
        std::string out2s = out2.str(); out2s.erase(out2s.end() - 1);
        out << out2s << ';';
      }
      // Remove the last delimiter
      std::string outs = out.str(); outs.erase(outs.end() - 1);
      consensus.setMetaValue("Intensities", outs);
    }

    UInt point_id = 0;
    for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
         point_it != pattern.points.end();
         ++point_it, ++point_id)
    {
      FeatureHandle point;
      point.setRT(point_it->rt);
      point.setMZ(point_it->mz);
      point.setUniqueId(point_id);

      consensus.insert(point);
    }

    return consensus;
  }

  void SILACAnalyzer::generateClusterFeatureByCluster(FeatureMap<> & out, const Clustering & clustering) const
  {
    for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it)
    {
      // RT value as weighted RT position of all peaks
      DoubleReal global_rt = 0;
      // Total intensity
      DoubleReal global_intensity = 0;

      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
           pattern_it != cluster_it->second.end();
           ++pattern_it)
      {
        SILACPattern & pattern = *pattern_it->second;

        for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = pattern.intensities.begin();
             shift_inten_it != pattern.intensities.end();
             ++shift_inten_it)
        {
          for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
               peak_inten_it != shift_inten_it->end();
               ++peak_inten_it)
          {
            DoubleReal intensity = *peak_inten_it;

            // Add to RT value and global intensity
            global_rt += intensity * pattern.rt;
            global_intensity += intensity;
          }
        }
      }

      // Calculate global RT value
      global_rt /= global_intensity;

      SILACPattern & pattern_first = *cluster_it->second.begin()->second;

      for (UInt shift_id = 0; shift_id < pattern_first.mass_shifts.size(); ++shift_id)
      {
        // XXX: Feature detection produces a stray 0 mass shift
        if (shift_id > 0 && pattern_first.mass_shifts[shift_id] == 0)
          continue;

        Feature feature;

        // MZ value as weighted MZ position of monoisotopic peaks of given mass shift
        DoubleReal shift_mz = 0;
        // Total intensity
        DoubleReal shift_intensity = 0;
        // Total intensity of monoisotopic peak
        DoubleReal shift_intensity0 = 0;

        // Bounding box per peak
        std::map<UInt, DBoundingBox<2> > bboxs;

        for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
             pattern_it != cluster_it->second.end();
             ++pattern_it)
        {
          SILACPattern & pattern = *pattern_it->second;

          const std::vector<DoubleReal> & intensities = pattern.intensities[shift_id];
          DoubleReal mz = pattern.mz + pattern.mass_shifts[shift_id];
          DoubleReal intensity0 = intensities[0];

          // Add to MZ value and shift intensity of monoisotopic peak
          shift_mz += intensity0 * mz;
          shift_intensity0 += intensity0;

          // Iterator over every peak
          UInt peak_id = 0;
          std::vector<DoubleReal>::const_iterator peak_inten_it = intensities.begin();
          DoubleReal peak_mz = mz;
          for (;
               peak_inten_it != intensities.end();
               ++peak_id, ++peak_inten_it, peak_mz += 1. / pattern.charge)
          {
            shift_intensity += *peak_inten_it;
            bboxs[peak_id].enlarge(pattern.rt, peak_mz);
          }
        }

        // Add each bbox as convex hulls to the cluster
        for (std::map<UInt, DBoundingBox<2> >::const_iterator bboxs_it = bboxs.begin();
             bboxs_it != bboxs.end();
             ++bboxs_it)
        {
          ConvexHull2D hull;
          hull.addPoint(bboxs_it->second.min_);
          hull.addPoint(bboxs_it->second.max_);
          feature.getConvexHulls().push_back(hull);
        }

        // XXX: Real quality?
        feature.setOverallQuality(1);
        feature.setCharge(pattern_first.charge);

        // Calculate MZ value
        shift_mz /= shift_intensity0;

        feature.setRT(global_rt);
        feature.setMZ(shift_mz);
        feature.setIntensity(shift_intensity);

        out.push_back(feature);
      }
    }
  }

  void SILACAnalyzer::readFilterConsensusByPattern(ConsensusMap & in, vector<vector<SILACPattern> > & data)
  {
    std::map<std::pair<Int, Int>, std::vector<SILACPattern> > layers;

    for (ConsensusMap::const_iterator pattern_it = in.begin(); pattern_it != in.end(); ++pattern_it)
    {
      SILACPattern pattern;
      pattern.rt = pattern_it->getRT();
      pattern.mz = pattern_it->getMZ();
      pattern.charge = pattern_it->getCharge();
      pattern.quality = pattern_it->getQuality();

      pattern.isotopes_per_peptide = pattern_it->getMetaValue("Peaks per peptide");

      StringList text = StringList::create(pattern_it->getMetaValue("Mass shifts [Da]"), ';');
      pattern.mass_shifts.push_back(0);
      for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
      {
        pattern.mass_shifts.push_back(text_it->toDouble() / pattern.charge);
      }

      text = StringList::create(pattern_it->getMetaValue("Intensities"), ';');
      for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
      {
        StringList text2 = StringList::create(*text_it, ',');
        vector<DoubleReal> inten;
        for (StringList::const_iterator text2_it = text2.begin(); text2_it != text2.end(); ++text2_it)
        {
          inten.push_back(text2_it->toDouble());
        }
        pattern.intensities.push_back(inten);
      }

      for (ConsensusFeature::const_iterator point_it = pattern_it->begin(); point_it != pattern_it->end(); ++point_it)
      {
        SILACPoint point;
        point.rt = point_it->getRT();
        point.mz = point_it->getMZ();

        pattern.points.push_back(point);
      }

      layers[std::make_pair(Int(pattern.mass_shifts.at(1)), pattern.charge)].push_back(pattern);
    }

    for (std::map<std::pair<Int, Int>, std::vector<SILACPattern> >::iterator it = layers.begin(); it != layers.end(); ++it)
    {
      data.push_back(it->second);
    }
  }

  const String & SILACAnalyzer::selectColor(UInt nr)
  {
    // 15 HTML colors
    static const String colors[] =
    {
      "#00FFFF", "#000000", "#0000FF", "#FF00FF", "#008000",
      "#808080", "#00FF00", "#800000", "#000080", "#808000",
      "#800080", "#FF0000", "#C0C0C0", "#008080", "#FFFF00",
    };
    const Int colors_len = 15;

    return colors[nr % colors_len];
  }

}

