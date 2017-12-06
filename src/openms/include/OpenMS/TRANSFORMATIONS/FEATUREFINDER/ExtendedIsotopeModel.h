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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{
  /**
        @brief Extended isotope distribution approximated using linear interpolation.

    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query).
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine

    Peak widening is achieved by a Gaussian shape.

        @htmlinclude OpenMS_ExtendedIsotopeModel.parameters
    */
  class OPENMS_DLLAPI ExtendedIsotopeModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef InterpolationModel::CoordinateType IntensityType;

    enum Averagines {C = 0, H, N, O, S, AVERAGINE_NUM};

    /// Default constructor
    ExtendedIsotopeModel();

    ///  copy constructor
    ExtendedIsotopeModel(const ExtendedIsotopeModel & source);

    /// destructor
    ~ExtendedIsotopeModel() override;

    /// assignment operator
    virtual ExtendedIsotopeModel & operator=(const ExtendedIsotopeModel & source);

    UInt getCharge();

    /// create new ExtendedIsotopeModel object (needed by Factory)
    static BaseModel<1> * create()
    {
      return new ExtendedIsotopeModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "ExtendedIsotopeModel";
    }

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over.
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model.
    */
    void setOffset(CoordinateType offset) override;

    CoordinateType getOffset();

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /** @brief get the monoisotopic mass of the Isotope model
    */
    CoordinateType getCenter() const override;

protected:
    CoordinateType isotope_stdev_;
    UInt charge_;
    CoordinateType monoisotopic_mz_;
    double averagine_[AVERAGINE_NUM];
    Int max_isotope_;
    double trim_right_cutoff_;
    double isotope_distance_;

    void updateMembers_() override;
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEMODEL_H
