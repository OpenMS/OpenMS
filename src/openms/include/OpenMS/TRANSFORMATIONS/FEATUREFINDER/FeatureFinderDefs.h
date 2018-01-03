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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H


#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

namespace OpenMS
{

  // forward declaration
  class FeatureFinder;

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
        BaseException(file, line, function, "NoSuccessor", "no successor/predecessor"),
        index_(index)
      {
        what_ = String("there is no successor/predecessor for the given Index: ") + index_.first + "/" + index_.second;
        Exception::GlobalExceptionHandler::getInstance().setMessage(what_);
      }

      ~NoSuccessor() throw() override
      {
      }

protected:
      IndexPair index_;            // index without successor/predecessor
    };
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H
