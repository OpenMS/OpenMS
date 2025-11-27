// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <vector>

namespace OpenMS
{
  class BaseFeature;

  /**
    @brief Representation of a Peak2D, RichPeak2D or Feature .

    The position and the intensity of the referenced feature are stored in the base class Peak2D.
    The original datapoint is referenced by the map and unique id.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI FeatureHandle :
    public Peak2D,
    public UniqueIdInterface
  {

public:
    class FeatureHandleMutable_;

    ///@name Type definitions
    //@{
    /// Charge type
    typedef Int ChargeType;
    /// Feature width type
    typedef float WidthType;
    //@}

    ///@name Constructors and destructor
    //@{
    /// Default constructor
    FeatureHandle();

    /// Constructor with map index, element index and position
    FeatureHandle(UInt64 map_index, const Peak2D& point, UInt64 element_index);

    /// Constructor from map index and basic feature
    FeatureHandle(UInt64 map_index, const BaseFeature& feature);

    /// Copy constructor
    FeatureHandle(const FeatureHandle& rhs);

    /// Assignment operator
    FeatureHandle& operator=(const FeatureHandle& rhs);

    /// Destructor
    ~FeatureHandle() override;

    /**
      @brief Override (most of all) constness.

      We provide this such that you can modify instances FeatureHandle which are
      stored within a ConsensusFeature.  Note that std::set does not provide
      non-const iterators, because these could be used to change the relative
      ordering of the elements, and iterators are (by design/concept) unaware of
      their containers.  Since ConsensusFeature uses the ordering by IndexLess
      (which see), you <i>must not</i> modify the map index of element index if
      there is more than one FeatureHandle stored in a ConsensusFeature.
      Consequently, we "disable" setMapIndex() or setElementIndex() even within
      FeatureHandleMutable_.  On the other hand, it is perfectly safe to apply
      FeatureHandle::setRT(), FeatureHandle::setMZ(),
      FeatureHandle::setIntensity(), FeatureHandle::setCharge(), etc..
    */
    FeatureHandleMutable_& asMutable() const;
    //@}

    ///@name Accessors
    //@{
    /// Returns the map index
    UInt64 getMapIndex() const;

    /// Set the map index
    void setMapIndex(UInt64 i);

    /// Sets the charge
    void setCharge(ChargeType charge);

    /// Returns the charge
    ChargeType getCharge() const;

    /// Sets the width (FWHM)
    void setWidth(WidthType width);

    /// Returns the width (FWHM)
    WidthType getWidth() const;

    //@}

    /// Equality operator
    bool operator==(const FeatureHandle& i) const;

    /// Equality operator
    bool operator!=(const FeatureHandle& i) const;

    /// Comparator by map and unique id
    struct IndexLess
    {
      bool operator()(FeatureHandle const& left, FeatureHandle const& right) const;
    };

protected:

    /// Index of the element's container
    UInt64 map_index_;
    /// Charge of the feature
    Int charge_;
    /// Width of the feature (FWHM)
    float width_;
  };

  /**
    @brief Helper class returned by FeatureHandle::asMutable(), which see.

    Note that the mutators for unique id and map index are declared private.
    This is done because these are used by IndexLess comparator.  This way it is
    a bit harder to use FeatureHandle::asMutable() for illegal purposes ;-)
   */
  class OPENMS_DLLAPI FeatureHandle::FeatureHandleMutable_ :
    public FeatureHandle
  {
private:
    using FeatureHandle::setUniqueId;
    using FeatureHandle::setMapIndex;
    FeatureHandleMutable_();
    FeatureHandleMutable_(const FeatureHandleMutable_&);
  };

  inline FeatureHandle::FeatureHandleMutable_& FeatureHandle::asMutable() const
  {
    // the const cast is to remove constness, but note that FeatureHandleMutable_ lacks some mutators
    // TODO use const_cast
    return static_cast<FeatureHandleMutable_&>(const_cast<FeatureHandle&>(*this));
  }

  inline bool FeatureHandle::IndexLess::operator()(FeatureHandle const& left, FeatureHandle const& right) const
  {
    // if map indices are equal, use unique ids
    if (left.map_index_ == right.map_index_)
    {
      return left.getUniqueId() < right.getUniqueId();
    }
    //else use map indices
    return left.map_index_ < right.map_index_;
  }

  ///Print the contents of a FeatureHandle to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const FeatureHandle& cons);
} // namespace OpenMS

