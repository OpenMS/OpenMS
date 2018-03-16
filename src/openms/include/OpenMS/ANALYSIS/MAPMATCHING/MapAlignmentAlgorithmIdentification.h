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
// $Maintainer: Hendrik Weisser $
// $Authors: Eva Lange, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <cmath> // for "abs"
#include <limits> // for "max"
#include <map>

namespace OpenMS
{
  /**
    @brief A map alignment algorithm based on peptide identifications from MS2 spectra.

    PeptideIdentification instances are grouped by sequence of the respective best-scoring PeptideHit and retention time data is collected (PeptideIdentification::getRT()).
    ID groups with the same sequence in different maps represent points of correspondence between the maps and form the basis of the alignment.

    Each map is aligned to a reference retention time scale.
    This time scale can either come from a reference file (@p reference parameter) or be computed as a consensus of the input maps (median retention times over all maps of the ID groups).
    The maps are then aligned to this scale as follows:\n
    The median retention time of each ID group in a map is mapped to the reference retention time of this group.
    Cubic spline smoothing is used to convert this mapping to a smooth function.
    Retention times in the map are transformed to the consensus scale by applying this function.

    @htmlinclude OpenMS_MapAlignmentAlgorithmIdentification.parameters

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI MapAlignmentAlgorithmIdentification :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MapAlignmentAlgorithmIdentification();

    /// Destructor
    ~MapAlignmentAlgorithmIdentification() override;

    // Set a reference for the alignment
    template <typename DataType> void setReference(DataType& data)
    {
      reference_.clear();
      if (data.empty()) return; // empty input resets the reference
      SeqToList rt_data;
      bool sorted = getRetentionTimes_(data, rt_data);
      computeMedians_(rt_data, reference_, sorted);
      if (reference_.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not extract retention time information from the reference file");
      }
    }

    /**
      @brief Align feature maps, consensus maps, peak maps, or peptide identifications.

      @param data Vector of input data (FeatureMap, ConsensusMap, PeakMap or @p vector<PeptideIdentification>) that should be aligned.
      @param transformations Vector of RT transformations that will be computed.
      @param reference_index Index in @p data of the reference to align to, if any
    */
    template <typename DataType>
    void align(std::vector<DataType>& data,
               std::vector<TransformationDescription>& transformations,
               Int reference_index = -1)
    {
      checkParameters_(data.size());
      startProgress(0, 3, "aligning maps");

      reference_index_ = reference_index;
      // is reference one of the input files?
      bool use_internal_reference = (reference_index >= 0);
      if (use_internal_reference)
      {
        if (reference_index >= data.size())
        {
          throw Exception::IndexOverflow(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, reference_index,
                                         data.size());
        }
        setReference(data[reference_index]);
      }

      // one set of RT data for each input map, except reference (if any):
      std::vector<SeqToList> rt_data(data.size() - use_internal_reference);
      bool all_sorted = true;
      for (Size i = 0, j = 0; i < data.size(); ++i)
      {
        if ((reference_index >= 0) && (i == Size(reference_index)))
        {
          continue; // skip reference map, if any
        }
        all_sorted &= getRetentionTimes_(data[i], rt_data[j++]);
      }
      setProgress(1);

      computeTransformations_(rt_data, transformations, all_sorted);
      setProgress(2);

      setProgress(3);
      endProgress();
    }

protected:

    /// Type to store retention times given for individual peptide sequences
    typedef std::map<String, DoubleList> SeqToList;

    /// Type to store one representative retention time per peptide sequence
    typedef std::map<String, double> SeqToValue;

    /// Index of input file to use as reference (if any)
    Int reference_index_;

    /// Reference retention times (per peptide sequence)
    SeqToValue reference_;

    /// Minimum number of runs a peptide must occur in
    Size min_run_occur_;

    /**
      @brief Compute the median retention time for each peptide sequence

      @param rt_data Lists of RT values for diff. peptide sequences (input, will be sorted)
      @param medians Median RT values for the peptide sequences (output)
      @param sorted Are RT lists already sorted?

      @throw Exception::IllegalArgument if the input list is empty
    */
    void computeMedians_(SeqToList& rt_data, SeqToValue& medians,
                         bool sorted = false);

    /**
      @brief Collect retention time data ("RT" MetaInfo) from peptide IDs

      @param peptides Input peptide IDs (lists of peptide hits will be sorted)
      @param rt_data Lists of RT values for diff. peptide sequences (output)

      @return Are the RTs already sorted? (Here: false)
    */
    bool getRetentionTimes_(std::vector<PeptideIdentification>& peptides,
                            SeqToList& rt_data);

    /**
      @brief Collect retention time data ("RT" MetaInfo) from peptide IDs annotated to spectra

      @param experiment Input map for RT data
      @param rt_data Lists of RT values for diff. peptide sequences (output)

      @return Are the RTs already sorted? (Here: false)
    */
    bool getRetentionTimes_(PeakMap& experiment, SeqToList& rt_data);

    /**
      @brief Collect retention time data ("RT" MetaInfo) from peptide IDs contained in feature maps or consensus maps

      The following global flags (mutually exclusive) influence the processing:\n
      Depending on @p use_unassigned_peptides, unassigned peptide IDs are used in addition to IDs annotated to features.\n
      Depending on @p use_feature_rt, feature retention times are used instead of peptide retention times.

      @param features Input features for RT data
      @param rt_data Lists of RT values for diff. peptide sequences (output)

      @return Are the RTs already sorted? (Here: true)
    */
    template <typename MapType>
    bool getRetentionTimes_(MapType& features, SeqToList& rt_data)
    {
      bool use_feature_rt = param_.getValue("use_feature_rt").toBool();
      for (typename MapType::Iterator feat_it = features.begin();
           feat_it != features.end(); ++feat_it)
      {
        if (use_feature_rt)
        {
          // find the peptide ID closest in RT to the feature centroid:
          String sequence;
          double rt_distance = std::numeric_limits<double>::max();
          bool any_hit = false;
          for (std::vector<PeptideIdentification>::iterator pep_it =
                 feat_it->getPeptideIdentifications().begin(); pep_it !=
                 feat_it->getPeptideIdentifications().end(); ++pep_it)
          {
            if (!pep_it->getHits().empty())
            {
              any_hit = true;
              double current_distance = fabs(pep_it->getRT() -
                                             feat_it->getRT());
              if (current_distance < rt_distance)
              {
                pep_it->sort();
                sequence = pep_it->getHits()[0].getSequence().toString();
                rt_distance = current_distance;
              }
            }
          }

          if (any_hit) rt_data[sequence].push_back(feat_it->getRT());
        }
        else
        {
          getRetentionTimes_(feat_it->getPeptideIdentifications(), rt_data);
        }
      }

      if (!use_feature_rt &&
          param_.getValue("use_unassigned_peptides").toBool())
      {
        getRetentionTimes_(features.getUnassignedPeptideIdentifications(),
                           rt_data);
      }

      // remove duplicates (can occur if a peptide ID was assigned to several
      // features due to overlap or annotation tolerance):
      for (SeqToList::iterator rt_it = rt_data.begin(); rt_it != rt_data.end();
           ++rt_it)
      {
        DoubleList& rt_values = rt_it->second;
        sort(rt_values.begin(), rt_values.end());
        DoubleList::iterator it = unique(rt_values.begin(), rt_values.end());
        rt_values.resize(it - rt_values.begin());
      }
      return true; // RTs were already sorted for duplicate detection
    }

    /**
      @brief Compute retention time transformations from RT data grouped by peptide sequence

      @param rt_data Lists of RT values for diff. peptide sequences, per dataset (input, will be sorted)
      @param transforms Resulting transformations, per dataset (output)
      @param sorted Are RT lists already sorted?
    */
    void computeTransformations_(std::vector<SeqToList>& rt_data,
                                 std::vector<TransformationDescription>&
                                 transforms, bool sorted = false);

    /**
      @brief Check that parameter values are valid

      Currently only 'min_run_occur' is checked.

      @param runs Number of runs (input files) to be aligned
    */
    void checkParameters_(const Size runs);

    /**
      @brief Get reference retention times

      If a reference file is supplied via the @p reference parameter, extract retention time
      information and store it in #reference_.
    */
    void getReference_();

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentAlgorithmIdentification(const MapAlignmentAlgorithmIdentification&);

    ///Assignment operator intentionally not implemented -> private
    MapAlignmentAlgorithmIdentification& operator=(const MapAlignmentAlgorithmIdentification&);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H
