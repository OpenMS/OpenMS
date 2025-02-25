// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

namespace OpenMS
{
  /**@brief The purpose of this struct is to provide definitions of classes and typedefs which are used throughout all FeatureFinder classes.  */
  struct OPENMS_DLLAPI FeatureFinderDefs
  {
    /// Index to peak consisting of two UInts (scan index / peak index)
    typedef IsotopeCluster::IndexPair IndexPair;

    /// Index to peak consisting of two UInts (scan index / peak index) with charge information
    typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;

    /// A set of peak indices
    typedef IsotopeCluster::IndexSet IndexSet;

    /// Flags that indicate if a peak is already used in a feature
    enum Flag {UNUSED, USED};

    /// Exception that is thrown if a method an invalid IndexPair is given
    class OPENMS_DLLAPI NoSuccessor :
      public Exception::BaseException
    {
public:
      NoSuccessor(const char * file, int line, const char * function, const IndexPair & index) :
        BaseException(file, line, function, "NoSuccessor", String("there is no successor/predecessor for the given Index: ") + String(index.first) + "/" + String(index.second)),
        index_(index)
      {
        Exception::GlobalExceptionHandler::setMessage(what());
      }

      ~NoSuccessor() noexcept override = default;

protected:
      IndexPair index_; // index without successor/predecessor
    };
  };
}
