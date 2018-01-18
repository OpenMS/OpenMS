// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_MASSDECOMPOSITIONALGORITHM_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_MASSDECOMPOSITIONALGORITHM_H

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

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_MASSDECOMPOSITIONALGORITHM_H
