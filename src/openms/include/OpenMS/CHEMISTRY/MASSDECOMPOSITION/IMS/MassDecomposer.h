// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_MASSDECOMPOSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_MASSDECOMPOSER_H

#include <vector>

namespace OpenMS
{

  namespace ims
  {
    /**
      @brief An interface to handle decomposing of integer values/masses
      over a set of integer weights (alphabet).

      An interface that addresses the following "mass decomposition" problems:

      - Existence Problem (whether the decomposition of the given mass exists),
      - One Decomposition Problem (returns one possible decomposition),
      - All Decompositions Problem (returns all possible decompositions),
      - Number of Decompositions Problem (returns the number of possible
       decompositions).

      Those problems are solved in integer arithmetic, i.e. only exact
      solutions are found with no error allowed.

      @param ValueType Type of values to be decomposed.
      @param DecompositionValueType Type of decomposition elements.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    template <typename ValueType,
              typename DecompositionValueType>
    class MassDecomposer
    {
public:
      /**
        Type of value to be decomposed.
      */
      typedef ValueType value_type;

      /**
        Type of decomposition value.
      */
      typedef DecompositionValueType decomposition_value_type;

      /**
        Type of decomposition container.
      */
      typedef std::vector<decomposition_value_type> decomposition_type;

      /**
        Type of container for many decompositions.
      */
      typedef std::vector<decomposition_type> decompositions_type;

      /**
        A virtual destructor.
      */
      virtual ~MassDecomposer(){}

      /**
        Returns true if the decomposition for the given @c mass exists, otherwise - false.

        @param mass Mass to be checked on decomposing.
        @return true, if the decomposition for @c mass exist, otherwise - false.
      */
      virtual bool exist(value_type mass) = 0;

      /**
        Returns one possible decomposition of the given @c mass.

        @param mass Mass to be decomposed.
        @return The decomposition of the @c mass, if one exists, otherwise - an empty container.
      */
      virtual decomposition_type getDecomposition(value_type mass) = 0;

      /**
        Returns all possible decompositions for the given @c mass.

        @param mass Mass to be decomposed.
        @return All possible decompositions of the @c mass, if there are any exist,
        otherwise - an empty container.
      */
      virtual decompositions_type getAllDecompositions(value_type mass) = 0;

      /**
        Returns the number of possible decompositions for the given @c mass.
        *
        @param mass Mass to be decomposed.
        @return The number of possible decompositions for the @c mass.
      */
      virtual decomposition_value_type getNumberOfDecompositions(value_type mass) = 0;

    };

  } // namespace ims
} // namespace OpenMS

#endif //OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_MASSDECOMPOSER_H
