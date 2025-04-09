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

/// A node of the kD-tree with pointer to corresponding data and index
class OPENMS_DLLAPI KDTreeFeatureNode
{

public:

  /// Constructor
  KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx);

  /// Copy constructor - copy the pointer, use same data object
  KDTreeFeatureNode(const KDTreeFeatureNode& rhs);

  /// Assignment operator - copy the pointer, use same data object
  KDTreeFeatureNode& operator=(KDTreeFeatureNode const& rhs);

  /// Destructor
  virtual ~KDTreeFeatureNode();

  /// libkdtree++ needs this typedef
  typedef double value_type;

  /// Needed for 2D range queries using libkdtree++. [0] returns RT, [1] m/z.
  value_type operator[](Size i) const;

  /// Return index of corresponding feature in data_
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

