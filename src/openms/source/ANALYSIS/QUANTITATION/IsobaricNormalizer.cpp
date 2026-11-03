// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <map>

namespace OpenMS
{

  IsobaricNormalizer::IsobaricNormalizer(const IsobaricQuantitationMethod* const quant_method) :
    quant_meth_(quant_method)
  {
    reference_channel_name_ = quant_meth_->getChannelInformation()[quant_meth_->getReferenceChannel()].name;
  }

  IsobaricNormalizer::IsobaricNormalizer(const IsobaricNormalizer& other) :
    quant_meth_(other.quant_meth_),
    reference_channel_name_(other.reference_channel_name_)
  {
  }

  IsobaricNormalizer& IsobaricNormalizer::operator=(const IsobaricNormalizer& rhs)
  {
    if (this == &rhs)
      return *this;

    quant_meth_ = rhs.quant_meth_;
    reference_channel_name_ = rhs.reference_channel_name_;

    return *this;
  }

  ConsensusFeature::HandleSetType::iterator IsobaricNormalizer::findReferenceChannel_(ConsensusFeature& cf, const ConsensusMap& consensus_map) const
  {
    for (ConsensusFeature::HandleSetType::iterator it_elements = cf.begin();
         it_elements != cf.end();
         ++it_elements)
    {
      if (consensus_map.getColumnHeaders().find(it_elements->getMapIndex())->second.getMetaValue("channel_name") == reference_channel_name_)
      {
        return it_elements;
      }
    }

    return cf.end();
  }

  void IsobaricNormalizer::buildVectorIndex_(const ConsensusMap& consensus_map)
  {
    // clear old values
    ref_vec_idx_ = 0;
    map_to_vec_index_.clear();

    Size index = 0;
    for (ConsensusMap::ColumnHeaders::const_iterator file_it = consensus_map.getColumnHeaders().begin();
         file_it != consensus_map.getColumnHeaders().end();
         ++file_it)
    {
      // TODO I think this does not work if you are normalizing a consensus map including multiple files
      //  which all have that reference channel name.
      if (file_it->second.getMetaValue("channel_name") == reference_channel_name_)
      {
        // the vector index, not the map index: they differ if the map indices do not start at 0
        // (e.g. the channels of a later file of a multi-file experiment)
        ref_vec_idx_ = index;
      }
      map_to_vec_index_[file_it->first] = index;
      ++index;
    }
  }

  void IsobaricNormalizer::collectRatios_(const ConsensusFeature& cf, const Peak2D::IntensityType& ref_intensity)
  {
    for (ConsensusFeature::HandleSetType::const_iterator it_elements = cf.begin();
         it_elements != cf.end();
         ++it_elements)
    {
      if (ref_intensity == 0) //avoid nan's and inf's
      {
        if (it_elements->getIntensity() == 0) // 0/0 will give 'nan'
        {
          //so leave it out completely (there is no information to be gained)
        }
        else // x/0 is 'inf' but std::sort() has problems with that
        {
          peptide_ratios_[map_to_vec_index_[it_elements->getMapIndex()]].push_back(std::numeric_limits<Peak2D::IntensityType>::max());
        }
      }
      else // everything seems fine
      {
        peptide_ratios_[map_to_vec_index_[it_elements->getMapIndex()]].push_back(it_elements->getIntensity() / ref_intensity);
      }

      // control
      peptide_intensities_[map_to_vec_index_[it_elements->getMapIndex()]].push_back(it_elements->getIntensity());
    }

  }

  void IsobaricNormalizer::computeNormalizationFactors_(std::vector<Peak2D::IntensityType>& normalization_factors)
  {
    // ensure that the ref_(ratios|intensities) are sorted
    std::sort(peptide_ratios_[ref_vec_idx_].begin(), peptide_ratios_[ref_vec_idx_].end());
    std::sort(peptide_intensities_[ref_vec_idx_].begin(), peptide_intensities_[ref_vec_idx_].end());

    // reporting
    Peak2D::IntensityType max_deviation_from_control = 0;

    // find MEDIAN of ratios for each channel (store as 0th element in sorted vector)
    for (std::map<Size, Size>::const_iterator it_map = map_to_vec_index_.begin(); it_map != map_to_vec_index_.end(); ++it_map)
    {
      // this is solely for readability reasons, the compiler should optimize this anyway
      const Size vec_idx = it_map->second;

      // no ratio at all (e.g. the reference channel is empty in every feature): nothing to normalize by
      if (peptide_ratios_[vec_idx].empty())
      {
        OPENMS_LOG_WARN << "IsobaricNormalizer:  map-id " << (it_map->first) << " has no intensity ratio to the reference channel and is not normalized." << std::endl;
        normalization_factors[vec_idx] = 1.0;
        continue;
      }

      // sort vector (partial_sort might improve performance here)
      std::sort(peptide_ratios_[vec_idx].begin(), peptide_ratios_[vec_idx].end());

      // save median as first element
      normalization_factors[vec_idx] = peptide_ratios_[vec_idx][peptide_ratios_[vec_idx].size() / 2];

      // A median ratio of 0 (the channel has no signal in most features) or of the stand-in for x/0 (the
      // reference channel has no signal in most features, see collectRatios_) is nothing to normalize by
      // either: dividing by it would turn the intensities of the channel into NaN, infinity or (almost) 0.
      if (normalization_factors[vec_idx] == 0 || normalization_factors[vec_idx] == std::numeric_limits<Peak2D::IntensityType>::max())
      {
        OPENMS_LOG_WARN << "IsobaricNormalizer:  map-id " << (it_map->first) << " has no usable median intensity ratio to the reference channel"
                        << " (the channel or the reference channel lacks signal in most features) and is not normalized." << std::endl;
        normalization_factors[vec_idx] = 1.0;
        continue;
      }

      // sort control (intensities)
      std::sort(peptide_intensities_[vec_idx].begin(), peptide_intensities_[vec_idx].end());

      // find MEDIAN of control-method (intensities) for each channel
      peptide_intensities_[vec_idx][0] = peptide_intensities_[vec_idx][peptide_intensities_[vec_idx].size() / 2] /
                                         peptide_intensities_[ref_vec_idx_][peptide_intensities_[ref_vec_idx_].size() / 2];

      OPENMS_LOG_INFO << "IsobaricNormalizer:  map-id " << (it_map->first) << " has factor " << (normalization_factors[vec_idx]) << " (control: " << (peptide_intensities_[vec_idx][0]) << ")" << std::endl;

      Peak2D::IntensityType dev = (peptide_ratios_[vec_idx][0] - peptide_intensities_[vec_idx][0]) / normalization_factors[vec_idx];
      if (fabs(max_deviation_from_control) < fabs(dev))
      {
        max_deviation_from_control = dev;
      }
    }

    OPENMS_LOG_INFO << "IsobaricNormalizer: max ratio deviation of alternative method is " << (max_deviation_from_control * 100) << "%\n";
  }

  void IsobaricNormalizer::normalize(ConsensusMap& consensus_map)
  {
    // determine reference channel as vector index
    buildVectorIndex_(consensus_map);

    // build mapping of map_index to ratio_array_index
    peptide_ratios_.resize(quant_meth_->getNumberOfChannels());
    peptide_intensities_.resize(quant_meth_->getNumberOfChannels());

    //build up ratios for each peptide of non-reference channels
    ConsensusFeature::HandleSetType::iterator ref_it;

    for (ConsensusMap::Iterator cm_it = consensus_map.begin(); cm_it != consensus_map.end(); ++cm_it)
    {
      // find reference index (this is inefficient to do every time,
      // but the most robust against anyone who tries to change the internals of ConsensusFeature):
      ref_it = findReferenceChannel_(*cm_it, consensus_map);

      // reference channel not found in this ConsensusFeature
      if (ref_it == cm_it->end())
      {
        OPENMS_LOG_WARN << "IsobaricNormalizer::normalize() WARNING: ConsensusFeature "
                 << (cm_it - consensus_map.begin())
                 << " does not have a reference channel! Skipping"
                 << std::endl;
        continue;
      }

      collectRatios_(*cm_it, ref_it->getIntensity());
    } // ! collect ratios

    // vector to store the channel wise normalization factors
    std::vector<Peak2D::IntensityType> normalization_factors;
    normalization_factors.resize(quant_meth_->getNumberOfChannels());

    // compute the normalization factors based on the medians of the compute ratios
    computeNormalizationFactors_(normalization_factors);

    // free memory
    peptide_intensities_.clear();
    peptide_ratios_.clear();

    // adjust intensity ratios
    for (size_t i = 0; i < consensus_map.size(); ++i)
    {
      // find reference index (this is inefficient to do every time,
      // but the most robust against anyone who tries to change the
      // internals of ConsensusFeature):
      ref_it = findReferenceChannel_(consensus_map[i], consensus_map);

      // reference channel not found in this ConsensusFeature
      if (ref_it == consensus_map[i].end())
      {
        continue;
      }

      // now adjust the ratios, in place: the intensity does not take part in the order of the handles, and
      // copying the whole feature (with its peptide identifications) would be expensive
      for (ConsensusFeature::HandleSetType::iterator it_elements = consensus_map[i].begin();
           it_elements != consensus_map[i].end();
           ++it_elements)
      {
        if (it_elements == ref_it)
        {
          it_elements->asMutable().setIntensity(1);
        }
        else // divide current intensity by normalization factor (which was stored at position 0)
        {
          it_elements->asMutable().setIntensity(it_elements->getIntensity() / normalization_factors[map_to_vec_index_[it_elements->getMapIndex()]]);
        }
      }
    } // ! adjust ratios
  }

} // namespace
