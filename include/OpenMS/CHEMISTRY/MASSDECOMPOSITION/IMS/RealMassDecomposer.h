// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_REALMASSDECOMPOSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_REALMASSDECOMPOSER_H

#include <utility>
#include <map>
#include <memory>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>

namespace OpenMS
{
  namespace ims
  {

    /**
      @brief Handles decomposing of non-integer values/masses
      over a set of non-integer weights with an error allowed.

      Implements a decomposition of non-integer values with a certain error
      allowed. Exactness of decomposition can also be tuned by setting a
      precision factor for weights defining their scaling magnitude.

      Works in fact as a wrapper for classes that handle exact mass decomposing
      using integer arithmetics. Instead of decomposing a single value as done
      by integer mass decomposers, @c RealMassDecomposer defines a set of values
      that lie in the allowed range (defined by error and false negatives
      appeared due to rounding), scales those to integers, decomposes
      them using @c IntegerMassDecomposer, does some checks (i.e. on false
      positives appeared due to rounding) and collects decompositions together.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    class OPENMS_DLLAPI RealMassDecomposer
    {
public:

      /// Type of integer decomposer.
      typedef IntegerMassDecomposer<> integer_decomposer_type;

      /// Type of integer values that are decomposed.
      typedef integer_decomposer_type::value_type integer_value_type;

      /// Type of result decompositions from integer decomposer.
      typedef integer_decomposer_type::decompositions_type decompositions_type;

      /// Type of the number of decompositions.
      typedef unsigned long long number_of_decompositions_type;

      typedef std::map<unsigned int, std::pair<unsigned int, unsigned int> > constraints_type;

      /**
        Constructor with weights.

        @param weights Weights over which values/masses to be decomposed.
      */
      explicit RealMassDecomposer(const Weights & weights);

      /**
        Gets all decompositions for a @c mass with an @c error allowed.

        @param mass Mass to be decomposed.
        @param error Error allowed between given and result decomposition.
        @return All possible decompositions for a given mass and error.
      */
      decompositions_type getDecompositions(double mass, double error);

      decompositions_type getDecompositions(double mass, double error, const constraints_type & constraints);

      /**
       Gets a number of all decompositions for a @c mass with an @c error
       allowed. It's similar to the @c getDecompositions(double,double) function
       but less space consuming, since doesn't use container to store decompositions.

       @param mass Mass to be decomposed.
       @param error Error allowed between given and result decomposition.
       @return Number of all decompositions for a given mass and error.
      */
      number_of_decompositions_type getNumberOfDecompositions(double mass, double error);

private:
      /// Weights over which values/masses to be decomposed.
      Weights weights_;

      /// Minimal and maximal rounding errors.
      std::pair<double, double> rounding_errors_;

      /// Precision to scale double values to integer
      double precision_;

      /**
        Decomposer to be used for exact decomposing using
        integer arithmetics.
      */
      std::auto_ptr<integer_decomposer_type> decomposer_;
    };

  } // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_REALMASSDECOMPOSER_H
