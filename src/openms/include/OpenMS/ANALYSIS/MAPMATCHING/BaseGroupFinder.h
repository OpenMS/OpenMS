// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
    @brief The base class of all element group finding algorithms.

    This class defines the basic interface for all element group finding
    algorithms.

    All derived algorithms take one or several consensus maps and find
    corresponding features across the maps (or within one map). They return one
    consensus map containing the found consensus features.

    The element indices of the result consensus features are the container
    access indices of the input maps. The map indices of the result consensus
    features are are the indices in the input map vector.
  */
  class OPENMS_DLLAPI BaseGroupFinder :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    BaseGroupFinder();

    /// Destructor
    ~BaseGroupFinder() override;

    /**
      @brief Run the algorithm

      @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    virtual void run(const std::vector<ConsensusMap> & input, ConsensusMap & result) = 0;



protected:

    /**
      @brief Checks if all file descriptions have disjoint map identifiers

      @exception Exception::IllegalArgument Is thrown if a file id is found twice
    */
    void checkIds_(const std::vector<ConsensusMap> & maps) const;

private:

    /// Copy constructor intentionally not implemented
    BaseGroupFinder(const BaseGroupFinder &);

    /// Assignment operator intentionally not implemented
    BaseGroupFinder & operator=(const BaseGroupFinder &);

  };

} // namespace OpenMS

