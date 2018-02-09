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
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H

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

    /// Register all derived classes here
    static void registerChildren();

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

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H
