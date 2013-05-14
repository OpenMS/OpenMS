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

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACANALYZER_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACANALYZER_H

//OpenMS includes
#include <OpenMS/config.h>
// #include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/METADATA/MSQuantifications.h>

#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace std;

namespace OpenMS
{
  /**
   * @brief Algorithm to use for SILACAnalysis
   *
   * Please initialize before usage using the initialize function.
   * Next, one can estimate the peak width before filtering the data. 
   *
   * In the final step, clusterData will generate the output data.
   *
   * @see SILACFiltering
   */
  class OPENMS_DLLAPI SILACAnalyzer :
    public ProgressLogger
  {
  private:

    // input and output files
    String in;
    String out;
    String out_clusters;
    String out_features;
    String out_mzq;

    String out_filters;
    String in_filters;
    String out_debug;

    // section "sample"
    String selected_labels;
    UInt charge_min;
    UInt charge_max;
    Int missed_cleavages;
    UInt isotopes_per_peptide_min;
    UInt isotopes_per_peptide_max;

    // section "algorithm"
    DoubleReal rt_threshold;
    DoubleReal rt_min;
    DoubleReal intensity_cutoff;
    DoubleReal intensity_correlation;
    DoubleReal model_deviation;
    bool allow_missing_peaks;

    // section "labels"
    std::vector<std::vector<String> > SILAClabels;    // list of SILAC labels, e.g. selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"
    std::vector<std::vector<DoubleReal> > massShifts; // list of mass shifts

    typedef SILACClustering Clustering;

    MSQuantifications msq;

  public:
    SILACAnalyzer() :
      allow_missing_peaks(true)
    {
    }

    /**
     * @brief Initializes the algorithm with parameters
     *
     * @param selected_labels_ Labels used for labelling the sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. The used labels must be described in the label_identifiers parameter
     * @param label_identifiers a map of labels and corresponding mass shits e.g. "Arg6" => 6.0201290268
     
     */
    void initialize(
      // section "sample"
      String selected_labels_,
      UInt charge_min_,
      UInt charge_max_,
      Int missed_cleavages_,
      UInt isotopes_per_peptide_min_,
      UInt isotopes_per_peptide_max_, 
      // section "algorithm"
      DoubleReal rt_threshold_,
      DoubleReal rt_min_,
      DoubleReal intensity_cutoff_,
      DoubleReal intensity_correlation_,
      DoubleReal model_deviation_,
      bool allow_missing_peaks_,
      // labels part
      map<String, DoubleReal> label_identifiers)
    {
      selected_labels          = selected_labels_;
      charge_min               = charge_min_;
      charge_max               = charge_max_;
      missed_cleavages         = missed_cleavages_;
      isotopes_per_peptide_min = isotopes_per_peptide_min_;
      isotopes_per_peptide_max = isotopes_per_peptide_max_;

      rt_threshold           = rt_threshold_;
      rt_min                 = rt_min_;
      intensity_cutoff       = intensity_cutoff_;
      intensity_correlation  = intensity_correlation_;
      model_deviation        = model_deviation_;
      allow_missing_peaks    = allow_missing_peaks_;

      calculateLabelsAndMassShifts(label_identifiers);
    }

    /**
     * @brief Calculate the internal massShift and label datastructures from a
     * map of identifiers and mass shifts (e.g. "Arg6" => 6.0201290268).
     *
     * Note that this part is necessary for the algorithm to work and this
     * function is called from the initialize section.  The algorithm has to be
     * initialized _first_ with the list of actually used labels
     * (selected_labels) using the initialize_sample call.
     */
    void calculateLabelsAndMassShifts(map<String, DoubleReal> label_identifiers);

    void run_all(MSExperiment<Peak1D> & exp, ConsensusMap & out_map)
    {
      PeakWidthEstimator::Result peak_width;
      vector<vector<SILACPattern> > data;
      MSQuantifications msq;
      vector<Clustering *> cluster_data;

      peak_width = estimatePeakWidth(exp);
      filterData(exp, peak_width, data); 
      clusterData(exp, peak_width, cluster_data, data);

      // write output to consensus map
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterConsensusByCluster(out_map, **it);
      }
    }

    /**
     * @brief Peak width estimation
     */
    PeakWidthEstimator::Result estimatePeakWidth(const MSExperiment<Peak1D> & exp);

    /**
     * @brief Filtering
     */
    void filterData(MSExperiment<Peak1D> & exp, const PeakWidthEstimator::Result & peak_width, vector<vector<SILACPattern> > & data);

    /**
     * @brief Clustering
     */
    void clusterData(const MSExperiment<> &, const PeakWidthEstimator::Result &, vector<Clustering *> &, vector<vector<SILACPattern> > & data);


    /**
     * @brief Returns the list of SILAC labels, e.g.
     * selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"
     */
    std::vector<std::vector<String> > & getSILAClabels()
    {
      return SILAClabels;
    }

    /**
     * @brief Returns the list of mass shifts calcualted
     */
    std::vector<std::vector<DoubleReal> > & getMassShifts()
    {
      return massShifts;
    }

  public:

    /**
     * @brief Generate ConsensusMap from clustering result
     */
    void generateClusterConsensusByCluster(ConsensusMap &, const Clustering &) const;

    /**
     * @brief Generate ConsensusMap from clustering result, one consensus per pattern
     */
    void generateClusterConsensusByPattern(ConsensusMap &, const Clustering &, UInt & cluster_id) const;

    /**
     * @brief Generate debug output from clustering result
     */
    void generateClusterDebug(std::ostream & out, const Clustering & clustering, UInt & cluster_id) const;

  public:
    /**
     * @brief Generate ConsensusMap from filter result
     */
    void generateFilterConsensusByPattern(ConsensusMap &, const std::vector<SILACPattern> &) const;

  private:

    /**
     * @brief Generate a consensus entry from a pattern
     */
    ConsensusFeature generateSingleConsensusByPattern(const SILACPattern &) const;


  public:

    /**
     * @brief Generate FeatureMap from clustering result
     */
    void generateClusterFeatureByCluster(FeatureMap<> &, const Clustering &) const;

    /**
     * @brief Read filter result from ConsensusMap
     */
    void readFilterConsensusByPattern(ConsensusMap &, vector<vector<SILACPattern> > &);

    static const String & selectColor(UInt nr);

    /**
     * @brief Read consensusXML from file to ConsensusMap
     */
    void readConsensus(const String & filename, ConsensusMap & in) const
    {
      ConsensusXMLFile c_file;
      c_file.load(filename, in);
    }

    /**
     * @brief Write consensusXML from ConsensusMap to file
     */
    void writeConsensus(const String & filename, ConsensusMap & out) const
    {
      out.sortByPosition();
      out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      out.setExperimentType("silac");

      ConsensusXMLFile c_file;
      c_file.store(filename, out);
    }

    /**
     * @brief Write MzQuantML from ConsensusMap to file
     */
    void writeMzQuantML(const String & filename, MSQuantifications & msq) const
    {
      //~ TODO apply above to ConsensusMap befor putting into Msq
      //~ out.sortByPosition();
      //~ out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
      //~ out.setExperimentType("SILAC");

      MzQuantMLFile file;
      file.store(filename, msq);
    }

    /**
     * @brief Write featureXML from FeatureMap to file
     */
    void writeFeatures(const String & filename, FeatureMap<> & out) const
    {
      out.sortByPosition();
      out.applyMemberFunction(&UniqueIdInterface::setUniqueId);

      FeatureXMLFile f_file;
      f_file.store(filename, out);
    }

  };
}

#endif /* SILACANALYZER_H_ */
