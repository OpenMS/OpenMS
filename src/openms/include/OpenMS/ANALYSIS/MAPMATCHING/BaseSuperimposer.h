// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
  /**
    @brief The base class of all superimposer algorithms.

    This class defines the basic interface for all superimposer algorithms. It
    works on several element maps and computes transformations that map the
    elements of the maps as near as possible to each other.
  */
  class OPENMS_DLLAPI BaseSuperimposer :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:

    /// Constructor
    BaseSuperimposer();

    /// Destructor
    ~BaseSuperimposer() override;

    /**
    @brief Estimates the transformation between input @p maps and returns the
    estimated transformations

    @exception IllegalArgument is thrown if the input maps are invalid.
    */
    virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation) = 0;



private:

    /// Copy constructor intentionally not implemented
    BaseSuperimposer(const BaseSuperimposer&);

    /// Assignment operator intentionally not implemented
    BaseSuperimposer& operator=(const BaseSuperimposer&);

  };

} // namespace OpenMS

