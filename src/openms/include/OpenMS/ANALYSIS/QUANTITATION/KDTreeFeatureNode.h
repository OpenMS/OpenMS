// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

namespace OpenMS
{

class KDTreeFeatureMaps;

/**
  @brief Lightweight kd-tree entry: a back-pointer to a @ref KDTreeFeatureMaps and the index of one feature in it.

  Used as the value type of the 2D kd-tree maintained by @ref KDTreeFeatureMaps. The node itself
  carries no feature data; instead, operator[]() forwards coordinate lookups (RT for axis 0, m/z
  for axis 1) to the parent @ref KDTreeFeatureMaps using the stored index. This keeps the kd-tree
  payload small and avoids duplicating per-feature attributes.

  Copy and assignment intentionally only copy the back-pointer; the wrapped data object is shared.
*/
class OPENMS_DLLAPI KDTreeFeatureNode
{

public:

  /**
    @brief Create a node referencing one feature inside @p data.

    @param[in] data Back-pointer to the @ref KDTreeFeatureMaps holding the feature; must outlive this node.
    @param[in] idx  Global index of the referenced feature inside @p data.
  */
  KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx);

  /// Copy constructor - copies the back-pointer; the wrapped data object is shared
  KDTreeFeatureNode(const KDTreeFeatureNode& rhs);

  /// Assignment operator - copies the back-pointer; the wrapped data object is shared
  KDTreeFeatureNode& operator=(KDTreeFeatureNode const& rhs);

  /// Destructor (does not delete the wrapped data)
  virtual ~KDTreeFeatureNode();

  /// libkdtree++ requires this typedef on the value type
  typedef double value_type;

  /// 2D coordinate accessor used by libkdtree++: @c [0] returns RT, @c [1] returns m/z (both forwarded to the parent @ref KDTreeFeatureMaps)
  value_type operator[](Size i) const;

  /// Global index of the referenced feature inside the parent @ref KDTreeFeatureMaps
  Size getIndex() const;

protected:

  /// Pointer to the actual data
  KDTreeFeatureMaps* data_;

  /// Index of this feature
  Size idx_;

private:

  /// Default constructor is not supposed to be called
  KDTreeFeatureNode();

};

}

