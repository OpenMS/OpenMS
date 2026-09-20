// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/AppliedProcessingStep.h>
#include <OpenMS/METADATA/ID/IdentifiedCompound.h>
#include <OpenMS/METADATA/ID/IdentifiedSequence.h>
#include <OpenMS/METADATA/ID/InputFile.h>
#include <OpenMS/METADATA/ID/Observation.h>
#include <OpenMS/METADATA/ID/ObservationMatch.h>
#include <OpenMS/METADATA/ID/ObservationMatchGroup.h>
#include <OpenMS/METADATA/ID/ParentGroup.h>
#include <OpenMS/METADATA/ID/ParentSequence.h>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index_container.hpp>

namespace OpenMS::IdentificationDataInternal
{
template<typename Value>
struct ContainerTraits;
template<>
struct ContainerTraits<InputFile>
{
  typedef boost::multi_index_container<
    InputFile,
    boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::member<InputFile, std::string, &InputFile::name>>>>
    Type;
};

template<>
struct ContainerTraits<Observation>
{
  typedef boost::multi_index_container<Observation,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::composite_key<
                                         Observation,
                                         boost::multi_index::member<Observation, InputFileRef, &Observation::input_file>,
                                         boost::multi_index::member<Observation, std::string, &Observation::data_id>>>>>
    Type;
};

template<>
struct ContainerTraits<ParentSequence>
{
  typedef boost::multi_index_container<ParentSequence,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
                                         boost::multi_index::member<ParentSequence, std::string, &ParentSequence::accession>>>>
    Type;
};

template<>
struct ContainerTraits<IdentifiedCompound>
{
  typedef boost::multi_index_container<IdentifiedCompound,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
                                         boost::multi_index::member<IdentifiedCompound, std::string, &IdentifiedCompound::identifier>>>>
    Type;
};

template<>
struct ContainerTraits<IdentifiedPeptide>
{
  typedef boost::multi_index_container<IdentifiedPeptide,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
                                         boost::multi_index::member<IdentifiedPeptide, AASequence, &IdentifiedPeptide::sequence>>>>
    Type;
};

template<>
struct ContainerTraits<IdentifiedOligo>
{
  typedef boost::multi_index_container<IdentifiedOligo,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
                                         boost::multi_index::member<IdentifiedOligo, NASequence, &IdentifiedOligo::sequence>>>>
    Type;
};

template<>
struct ContainerTraits<ObservationMatch>
{
  typedef boost::multi_index_container<ObservationMatch,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::composite_key<
                                         ObservationMatch,
                                         boost::multi_index::member<ObservationMatch, ObservationRef, &ObservationMatch::observation_ref>,
                                         boost::multi_index::member<ObservationMatch, IdentifiedMolecule, &ObservationMatch::identified_molecule_var>,
                                         boost::multi_index::member<ObservationMatch, AdductOpt, &ObservationMatch::adduct_opt>>>>>
    Type;
};

template<>
struct ContainerTraits<ParentGroup>
{
  typedef boost::multi_index_container<ParentGroup,
                                       boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
                                         boost::multi_index::member<ParentGroup, std::set<ParentSequenceRef>, &ParentGroup::parent_refs>>>>
    Type;
};

template<>
struct ContainerTraits<ObservationMatchGroup>
{
  typedef boost::multi_index_container<
    ObservationMatchGroup,
    boost::multi_index::indexed_by<boost::multi_index::ordered_unique<
      boost::multi_index::member<ObservationMatchGroup, std::set<ObservationMatchRef>, &ObservationMatchGroup::observation_match_refs>>>>
    Type;
};

template<>
struct ContainerTraits<AppliedProcessingStep>
{
  typedef boost::multi_index_container<
    AppliedProcessingStep,
    boost::multi_index::indexed_by<
      boost::multi_index::sequenced<>,
      boost::multi_index::ordered_unique<
        boost::multi_index::member<AppliedProcessingStep, std::optional<ProcessingStepRef>, &AppliedProcessingStep::processing_step_opt>>>>
    Type;
};

template<typename Value, typename Key, typename Prefix>
struct IDDataContainer<Value, Key, Prefix>::Impl
{
  typename ContainerTraits<Value>::Type records;
  auto& ordered()
  {
    if constexpr (std::is_same_v<Value, AppliedProcessingStep>) return records.template get<1>();
    else
      return records;
  }
  const auto& ordered() const
  {
    if constexpr (std::is_same_v<Value, AppliedProcessingStep>) return records.template get<1>();
    else
      return records;
  }
};

template<typename V, typename K, typename P>
IDDataContainer<V, K, P>::IDDataContainer() = default;
template<typename V, typename K, typename P>
IDDataContainer<V, K, P>::IDDataContainer(const IDDataContainer& other): impl_(other.impl_ ? std::make_unique<Impl>(*other.impl_) : nullptr)
{
}
template<typename V, typename K, typename P>
IDDataContainer<V, K, P>::IDDataContainer(IDDataContainer&& other) noexcept = default;
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::operator=(const IDDataContainer& other) -> IDDataContainer&
{
  if (this != &other)
  {
    IDDataContainer copy(other);
    swap(copy);
  }
  return *this;
}
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::operator=(IDDataContainer&& other) noexcept -> IDDataContainer& = default;
template<typename V, typename K, typename P>
IDDataContainer<V, K, P>::~IDDataContainer() = default;

template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::begin() const -> iterator
{ return iterator(impl_.get(), empty() ? nullptr : std::addressof(*impl_->records.begin())); }
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::orderedBegin_() const -> ordered_iterator
{ return ordered_iterator(impl_.get(), empty() ? nullptr : std::addressof(*impl_->ordered().begin())); }
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::size() const -> size_type
{ return impl_ ? impl_->records.size() : 0; }
template<typename V, typename K, typename P>
void IDDataContainer<V, K, P>::clear()
{
  if (impl_) impl_->records.clear();
}
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::insert(const V& value) -> std::pair<iterator, bool>
{
  if (! impl_) impl_ = std::make_unique<Impl>();
  auto result = [&]() {
    if constexpr (std::is_same_v<V, AppliedProcessingStep>) return impl_->records.push_back(value);
    else
      return impl_->records.insert(value);
  }();
  return {iterator(impl_.get(), std::addressof(*result.first)), result.second};
}
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::erase(iterator position) -> iterator
{
  auto& records = impl_->records;
  auto next = records.erase(records.iterator_to(*position));
  return iterator(impl_.get(), next == records.end() ? nullptr : std::addressof(*next));
}
template<typename V, typename K, typename P>
bool IDDataContainer<V, K, P>::modify(iterator position, const std::function<void(V&)>& modifier)
{ return impl_->records.modify(impl_->records.iterator_to(*position), modifier); }
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::find(const K& key) const -> iterator
{
  if (! impl_) return iterator(nullptr, nullptr);
  const auto& records = impl_->ordered();
  auto result = records.find(key);
  return iterator(impl_.get(), result == records.end() ? nullptr : std::addressof(*result));
}
template<typename V, typename K, typename P>
auto IDDataContainer<V, K, P>::equal_range(const P& key) const -> std::pair<iterator, iterator>
{
  if (! impl_) return {end(), end()};
  if constexpr (std::is_same_v<V, AppliedProcessingStep>)
  {
    // This container iterates in insertion order; its unique-key range must
    // also end at the next element in that order, not in the ordered index.
    auto first = find(key);
    auto last = first;
    if (last != end()) ++last;
    return {first, last};
  }
  const auto& records = impl_->ordered();
  auto range = records.equal_range(key);
  auto wrap = [&](const auto& it) { return iterator(impl_.get(), it == records.end() ? nullptr : std::addressof(*it)); };
  return {wrap(range.first), wrap(range.second)};
}
template<typename V, typename K, typename P>
const V* IDDataContainer<V, K, P>::advance_(const Impl* owner, const V* value, bool forward, bool ordered)
{
  auto advance = [&](const auto& records) -> const V* {
    auto it = value ? records.iterator_to(*value) : records.end();
    if (forward) ++it;
    else
      --it;
    return it == records.end() ? nullptr : std::addressof(*it);
  };
  if (ordered) return advance(owner->ordered());
  return advance(owner->records);
}

template class IDDataContainer<InputFile, std::string, std::string>;
template class IDDataContainer<Observation, std::tuple<InputFileRef, std::string>, InputFileRef>;
template class IDDataContainer<ParentSequence, std::string, std::string>;
template class IDDataContainer<IdentifiedCompound, std::string, std::string>;
template class IDDataContainer<IdentifiedPeptide, AASequence, AASequence>;
template class IDDataContainer<IdentifiedOligo, NASequence, NASequence>;
template class IDDataContainer<ObservationMatch, std::tuple<ObservationRef, IdentifiedMolecule, AdductOpt>, ObservationRef>;
template class IDDataContainer<ParentGroup, std::set<ParentSequenceRef>, std::set<ParentSequenceRef>>;
template class IDDataContainer<ObservationMatchGroup, std::set<ObservationMatchRef>, std::set<ObservationMatchRef>>;
template class IDDataContainer<AppliedProcessingStep, std::optional<ProcessingStepRef>, std::optional<ProcessingStepRef>>;
} // namespace OpenMS::IdentificationDataInternal
