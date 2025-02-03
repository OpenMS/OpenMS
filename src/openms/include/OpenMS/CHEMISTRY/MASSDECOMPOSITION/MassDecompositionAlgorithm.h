// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>

// ims includes
#ifdef OPENMS_COMPILER_MSVC
#pragma warning( push )
#pragma warning( disable : 4290 4267)
#endif

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
#ifdef OPENMS_COMPILER_MSVC
#pragma warning( pop )
#endif

#include <vector>

namespace OpenMS
{
  /**
    @brief Mass decomposition algorithm, given a mass it suggests possible compositions

    A mass decomposition algorithm decomposes a mass or a mass difference into
    possible amino acids and frequencies of them, which add up to the given mass.
    This class is a wrapper for the algorithm published in

    @htmlinclude OpenMS_MassDecompositionAlgorithm.parameters

    @ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI MassDecompositionAlgorithm :
    public DefaultParamHandler
  {
public:

    /**
      @name constructors and destructor
    */
    //@{
    /// Default constructor
    MassDecompositionAlgorithm();

    /// Destructor
    ~MassDecompositionAlgorithm() override;
    //@}

    /**
      @name Operators
    */
    //@{
    /// returns the possible decompositions given the weight
    void getDecompositions(std::vector<MassDecomposition> & decomps, double weight);
    //@}

protected:

    void updateMembers_() override;

    ims::IMSAlphabet * alphabet_;

    ims::RealMassDecomposer * decomposer_;

private:

    // will not be implemented
    /// Copy constructor
    MassDecompositionAlgorithm(const MassDecompositionAlgorithm & deco);

    /// assignment operator
    MassDecompositionAlgorithm & operator=(const MassDecompositionAlgorithm & rhs);
  };

} // namespace OpenMS

