// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <cmath>
#include <limits>
#include <map>

namespace OpenMS
{
  class OPENMS_DLLAPI MapAlignmentAlgorithmIdentification :
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    MapAlignmentAlgorithmIdentification();
    ~MapAlignmentAlgorithmIdentification() override;

    template <typename DataType>
    void setReference(const DataType& data)
    {
      reference_.clear();
      if (data.empty()) return;
      SeqToList rt_data;
      use_feature_rt_ = param_.getValue("use_feature_rt").toBool();
      score_cutoff_ = param_.getValue("score_cutoff").toBool();
      score_type_ = (std::string)param_.getValue("score_type");
      bool sorted = getRetentionTimes_(data, rt_data); // CHANGED: Ensures input is const-read
      computeMedians_(rt_data, reference_, sorted);

      if (reference_.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not extract retention time information from the reference file");
      }
    }

    template <typename DataType>
    void align(const std::vector<DataType>& data, // CHANGED: now `const`
               std::vector<TransformationDescription>& transformations,
               Int reference_index = -1)
    {
      checkParameters_(data.size());
      startProgress(0, 3, "aligning maps");

      reference_index_ = reference_index;
      bool use_internal_reference = (reference_index >= 0);
      if (use_internal_reference)
      {
        if (reference_index >= Int(data.size()))
        {
          throw Exception::IndexOverflow(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION,
                                         reference_index, data.size());
        }
        setReference(data[reference_index]); // CHANGED: safe as setReference also accepts const
      }

      std::vector<SeqToList> rt_data(data.size() - use_internal_reference);
      bool all_sorted = true;
      for (Size i = 0, j = 0; i < data.size(); ++i)
      {
        if ((reference_index >= 0) && (i == Size(reference_index)))
        {
          continue;
        }
        all_sorted &= getRetentionTimes_(data[i], rt_data[j++]); // CHANGED: const-read
      }

      setProgress(1);
      computeTransformations_(rt_data, transformations, all_sorted);
      setProgress(2);
      setProgress(3);
      endProgress();
    }

  protected:
    typedef std::map<String, DoubleList> SeqToList;
    typedef std::map<String, double> SeqToValue;

    Int reference_index_;
    SeqToValue reference_;
    Size min_run_occur_;
    bool use_feature_rt_{};
    bool use_adducts_{};
    double min_score_;
    bool score_cutoff_{};
    String score_type_;
    bool (*better_) (double, double) = [](double, double) { return true; };

    void computeMedians_(SeqToList& rt_data, SeqToValue& medians, bool sorted = false);

    bool getRetentionTimes_(const std::vector<PeptideIdentification>& peptides,
                            SeqToList& rt_data)
    {
      for (const auto& pep_id : peptides)
      {
        if (!pep_id.getHits().empty())
        {
          auto hits = pep_id.getHits();
          std::sort(hits.begin(), hits.end()); // CHANGED: sort local copy, not input
          if (better_(hits[0].getScore(), min_score_))
          {
            String seq = hits[0].getSequence().toString();
            rt_data[seq].push_back(pep_id.getRT());
          }
        }
      }
      return false;
    }

    bool getRetentionTimes_(IdentificationData& id_data, SeqToList& rt_data);

    bool getRetentionTimes_(PeakMap& experiment, SeqToList& rt_data);

    template <typename MapType>
    bool getRetentionTimes_(const MapType& features, SeqToList& rt_data)
    {
      if (!score_cutoff_)
      {
        better_ = [](double, double) { return true; };
      }
      else if (!features.empty() &&
               !features[0].getPeptideIdentifications().empty() &&
               features[0].getPeptideIdentifications()[0].isHigherScoreBetter())
      {
        better_ = [](double a, double b) { return a >= b; };
      }
      else
      {
        better_ = [](double a, double b) { return a <= b; };
      }

      for (auto feat_it = features.begin(); feat_it != features.end(); ++feat_it)
      {
        if (use_feature_rt_)
        {
          String sequence;
          double rt_distance = std::numeric_limits<double>::max();
          bool any_hit = false;
          for (const auto& pep_id : feat_it->getPeptideIdentifications())
          {
            if (!pep_id.getHits().empty())
            {
              any_hit = true;
              double current_distance = std::fabs(pep_id.getRT() - feat_it->getRT());
              if (current_distance < rt_distance)
              {
                auto hits = pep_id.getHits();
                std::sort(hits.begin(), hits.end()); // CHANGED: sort copy
                if (better_(hits[0].getScore(), min_score_))
                {
                  sequence = hits[0].getSequence().toString();
                  rt_distance = current_distance;
                }
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

      if (!use_feature_rt_ &&
          param_.getValue("use_unassigned_peptides").toBool())
      {
        getRetentionTimes_(features.getUnassignedPeptideIdentifications(), rt_data);
      }

      for (auto& rt_it : rt_data)
      {
        DoubleList& rt_values = rt_it.second;
        std::sort(rt_values.begin(), rt_values.end()); // CHANGED: only sorting local copies
        auto it = std::unique(rt_values.begin(), rt_values.end());
        rt_values.resize(std::distance(rt_values.begin(), it));
      }
      return true;
    }

    void computeTransformations_(std::vector<SeqToList>& rt_data,
                                 std::vector<TransformationDescription>& transforms,
                                 bool sorted = false);

    void checkParameters_(const Size runs);
    void getReference_();
    IdentificationData::ScoreTypeRef handleIdDataScoreType_(const IdentificationData& id_data);

  private:
    MapAlignmentAlgorithmIdentification(const MapAlignmentAlgorithmIdentification&);
    MapAlignmentAlgorithmIdentification& operator=(const MapAlignmentAlgorithmIdentification&);
  };

} // namespace OpenMS
