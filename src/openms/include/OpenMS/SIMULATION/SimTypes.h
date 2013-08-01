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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_SIMTYPES_H
#define OPENMS_SIMULATION_SIMTYPES_H

#include <vector>
#include <utility>
#include <map>
#include <utility>

#include <tr1/memory>

#include <boost/random/mersenne_twister.hpp>

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

// GSL includes (random number generation)
#include "OpenMS/MATH/GSL_WRAPPER/gsl_wrapper.h"

namespace OpenMS
{
  //forward declarations
  class SimRandomNumberGenerator;

  /// Coordinate type in mz and rt dimension
  typedef Peak2D::CoordinateType SimCoordinateType;

  /// Abundance of proteins/peptides
  typedef Peak2D::IntensityType SimIntensityType;

  /// Charge of a peptide
  typedef Feature::ChargeType SimChargeType;

  /// Raw data point
  typedef Peak1D SimPointType;

  /// Container for FASTAEntry & abundance information
  class OPENMS_DLLAPI SampleProteins : public std::vector<std::pair<FASTAFile::FASTAEntry, MetaInfoInterface> > { }; 

  /// Container for multiple channels of SampleProteins
  class OPENMS_DLLAPI SampleChannels : public std::vector<SampleProteins> { }; 

  /// Sim FeatureMap
  typedef FeatureMap<> FeatureMapSim;

  /// Sim FeatureMap Vector
  typedef std::vector<FeatureMapSim> FeatureMapSimVector;

  /// Sim MSExperiment type
  typedef MSExperiment<SimPointType> MSSimExperiment;

  //Sim Shared Pointer type
  typedef boost::shared_ptr< SimRandomNumberGenerator > MutableSimRandomNumberGeneratorPtr;
  /**
    @brief Wrapper class for random number generators used by the simulation classes

    The random numbers are separated two sources of randomness:

    <ul>
      <li><em>technical random numbers</em> which should represent technical
          sources of variability like instrument noise and </li>
      <li><em>biological random numbers</em> which should represent biological
          sources of variability (e.g. between two samples of the same composition)</li>
    </ul>

    @ingroup Simulation
  */
  class SimRandomNumberGenerator
  {
  public:

    boost::random::mt19937_64& getBiologicalRng()
    {
      return biological_rng_;
    }

    boost::random::mt19937_64& getTechnicalRng()
    {
      return technical_rng_;
    }

    void setBiologicalRngSeed( unsigned long int seed )
    {
      biological_rng_.seed(seed);
    }

    void setTechnicalRngSeed( unsigned long int seed )
    {
      technical_rng_.seed(seed);
    }

    /// Initialize the RNGs
    void initialize(bool biological_random, bool technical_random)
    {
      // use 0 as default seed to get reproducible experiments
      if (biological_random)
      {
        biological_rng_ = boost::random::mt19937_64(std::time(0));
      }
      else
      {
        biological_rng_ = boost::random::mt19937_64(0);
      }

      if (technical_random)
      {
        technical_rng_ = boost::random::mt19937_64(std::time(0));
      }
      else
      {
          technical_rng_ = boost::random::mt19937_64(0);
      }
    }
  private:
    /// random number generator for biological variability
    boost::random::mt19937_64 biological_rng_;
    /// random number generator for technical variability
    boost::random::mt19937_64 technical_rng_;

  };

}

#endif
