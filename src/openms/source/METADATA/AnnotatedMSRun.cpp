// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt, Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/METADATA/AnnotatedMSRun.h>

#include <algorithm>


namespace OpenMS
{
  bool AnnotatedMSRun::operator==(const AnnotatedMSRun& rhs) const
  {
    if (!(data == rhs.data && id_data_ == rhs.id_data_)) return false;
    for (Size index = 0; index < std::max(spectrum_queries_.size(), rhs.spectrum_queries_.size()); ++index)
    {
      if (getSpectrumIDQueries(index) != rhs.getSpectrumIDQueries(index)) return false;
    }
    return true;
  }

  MSExperiment& AnnotatedMSRun::getMSExperiment()
  {
    return data;
  }

  const MSExperiment& AnnotatedMSRun::getMSExperiment() const
  {
    return data;
  }

  void AnnotatedMSRun::setMSExperiment(MSExperiment&& experiment)
  {
    data = std::move(experiment);
  }

  void AnnotatedMSRun::setMSExperiment(const MSExperiment& experiment)
  {
    data = experiment;
  }
  
  const std::set<IdentificationData::QueryReference>& AnnotatedMSRun::getSpectrumIDQueries(Size index) const
  {
    static const std::set<IdentificationData::QueryReference> none;
    return index < spectrum_queries_.size() ? spectrum_queries_[index] : none;
  }

  void AnnotatedMSRun::addSpectrumIDQuery(Size index, const IdentificationData::QueryReference& query)
  {
    if (index >= data.getSpectra().size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, data.getSpectra().size());
    }
    if (index >= spectrum_queries_.size()) spectrum_queries_.resize(index + 1);
    spectrum_queries_[index].insert(query);
  }

  void AnnotatedMSRun::clearSpectrumIDQueries()
  {
    spectrum_queries_.clear();
  }

  std::vector<IdentificationData::QueryMatches> AnnotatedMSRun::getSpectrumIdentifications(Size index) const
  {
    auto identifications = id_data_.resolveLinks(getSpectrumIDQueries(index), {});
    for (auto& identification : identifications)
    {
      for (const auto& match : identification.query->getMatches())
      {
        identification.matches.push_back(&match);
      }
    }
    return identifications;
  }

  bool AnnotatedMSRun::hasSpectrumIdentifications() const
  {
    return std::any_of(spectrum_queries_.begin(), spectrum_queries_.end(), [](const auto& queries) { return ! queries.empty(); });
  }

  void AnnotatedMSRun::clearIdentifications()
  {
    id_data_.clear();
    spectrum_queries_.clear();
  }

  Size AnnotatedMSRun::eraseMatches(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&,
                                                             const IdentificationData::Match&)>& remove)
  {
    return id_data_.eraseMatches(remove).size();
  }

  Size AnnotatedMSRun::eraseIdentifications(const std::function<bool(const IdentificationData::Run&, const IdentificationData::Identification&)>& remove)
  {
    const auto erased = id_data_.eraseIdentifications(remove).first;
    for (auto& queries : spectrum_queries_)
    {
      std::erase_if(queries, [&](const IdentificationData::QueryReference& query) { return erased.contains(query); });
    }
    return erased.size();
  }

}

