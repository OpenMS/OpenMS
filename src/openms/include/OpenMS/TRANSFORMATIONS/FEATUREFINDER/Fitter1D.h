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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H

#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>

namespace OpenMS
{
  class InterpolationModel;

  /**
  @brief Abstract base class for all 1D-dimensional model fitter.

  Every derived class has to implement the static functions
  "T* create()" and "const String getProductName()" (see DefaultParamHandler for details)

  @htmlinclude OpenMS_Fitter1D.parameters

  @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI Fitter1D :
    public DefaultParamHandler,
    public FeatureFinderDefs
  {
public:

    /// IndexSet
    typedef IsotopeCluster::IndexSet IndexSet;
    /// IndexSet with charge information
    typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;
    /// Single coordinate
    typedef Feature::CoordinateType CoordinateType;
    /// Quality of a feature
    typedef Feature::QualityType QualityType;
    /// Raw data point type
    typedef Peak1D PeakType;
    /// Raw data container type using for the temporary storage of the input data
    typedef std::vector<PeakType> RawDataArrayType;
    /// Raw data iterator
    typedef RawDataArrayType::iterator PeakIterator;

    /// Default constructor.
    Fitter1D();

    /// copy constructor
    Fitter1D(const Fitter1D & source);

    /// destructor
    virtual ~Fitter1D()
    {
    }

    /// assignment operator
    virtual Fitter1D & operator=(const Fitter1D & source);

    /// return interpolation model
    virtual QualityType fit1d(const RawDataArrayType & /* range */, InterpolationModel * & /* model */);

    /// register all derived classes here
    static void registerChildren();

protected:

    /// standard derivation in bounding box
    CoordinateType tolerance_stdev_box_;
    /// minimum of the bounding box
    CoordinateType min_;
    /// maximum of the bounding box
    CoordinateType max_;
    /// standard derivation
    CoordinateType stdev1_;
    /// standard derivation
    CoordinateType stdev2_;
    /// basic statistics
    Math::BasicStatistics<> statistics_;
    /// interpolation step size
    CoordinateType interpolation_step_;

    virtual void updateMembers_();

  };

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FITTER1D_H
