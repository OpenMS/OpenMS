// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{

  /**
  @brief Encapsulated weight queries to simplify mono vs average weight computation

  Supports EmpiricalFormula's and AASequence's getMonoWeight() and getAverageWeight()

  */
  class OPENMS_DLLAPI WeightWrapper
  {

public:

    enum WEIGHTMODE {AVERAGE = 0, MONO, SIZE_OF_WEIGHTMODE};

    /**
    @brief constructor
    */
    WeightWrapper();

    /**
    @brief constructor
    */
    explicit WeightWrapper(const WEIGHTMODE weight_mode);

    /**
    @brief destructor
    */
    virtual ~WeightWrapper();

    /**
    @brief copy constructor
    */
    WeightWrapper(const WeightWrapper & source);


    /**
    @brief Sets the weight mode (MONO or AVERAGE)

    Sets the mode in which getWeight() calls are answered.

    */
    void setWeightMode(const WEIGHTMODE mode);


    /**
    @brief Gets the weight mode (MONO or AVERAGE)

    Gets the mode in which getWeight() calls are answered.

    */
    WEIGHTMODE getWeightMode() const;


    /**
    @brief returns the weight of either mono or average value

    Which weight is returned depends on the current weight-mode.

    @return double weight in u
    */
    double getWeight(const AASequence & aa) const;

    /**
    @brief returns the weight of either mono or average value

    Which weight is returned depends on the current weight-mode.

    @return double weight in u
    */
    double getWeight(const EmpiricalFormula & ef) const;


    /**
    @brief returns the weight of either mono or average value

    Which weight is returned depends on the current weight-mode.

    @return double weight in u
    */
    double getWeight(const Residue & r, Residue::ResidueType res_type = Residue::Full) const;


private:

    WEIGHTMODE weight_mode_;         ///< one of WeightWrapper::WEIGHTMODE's values


  };
}
