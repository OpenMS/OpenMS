// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H

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
    virtual ~BaseSuperimposer();

    /**
    @brief Estimates the transformation between input @p maps and returns the
    estimated transformations

    @exception IllegalArgument is thrown if the input maps are invalid.
    */
    virtual void run(const ConsensusMap & map_model, const ConsensusMap & map_scene, TransformationDescription & transformation) = 0;

    /// Register all derived classes here
    static void registerChildren();

private:

    /// Copy constructor intentionally not implemented
    BaseSuperimposer(const BaseSuperimposer &);

    /// Assignment operator intentionally not implemented
    BaseSuperimposer & operator=(const BaseSuperimposer &);

  };

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
