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
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>

namespace OpenMS
{
  class EmpiricalFormula;

  /**
        @brief Isotope distribution approximated using linear interpolation.

    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query).
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine

    Peak widening is achieved by either a Gaussian or Lorentzian shape.

        @htmlinclude OpenMS_IsotopeModel.parameters
    */
  class OPENMS_DLLAPI IsotopeModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef InterpolationModel::CoordinateType IntensityType;

    enum Averagines {C = 0, H, N, O, S, AVERAGINE_NUM};

    /// Default constructor
    IsotopeModel();

    ///  copy constructor
    IsotopeModel(const IsotopeModel & source);

    /// destructor
    ~IsotopeModel() override;

    /// assignment operator
    virtual IsotopeModel & operator=(const IsotopeModel & source);

    UInt getCharge();

    /// create new IsotopeModel object (needed by Factory)
    static BaseModel<1> * create()
    {
      return new IsotopeModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "IsotopeModel";
    }

    /** @brief set the offset of the model

        The whole model will be shifted to the new offset without being computing all over.
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model.
    */
    void setOffset(CoordinateType offset) override;

    CoordinateType getOffset();

    /// return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
    EmpiricalFormula getFormula();

    /// set sample/supporting points of interpolation
    virtual void setSamples(const EmpiricalFormula & formula);

    /// set sample/supporting points of interpolation (from base class)
    using InterpolationModel::setSamples;

    /** @brief get the center of the Isotope model

         This is a m/z-value not necessarily the monoisotopic mass.
    */
    CoordinateType getCenter() const override;

    /** @brief the Isotope distribution (without widening) from the last setSamples() call

      Useful to determine the number of isotopes that the model contains and their position

    */
    const IsotopeDistribution & getIsotopeDistribution() const;


protected:
    CoordinateType isotope_stdev_;
    CoordinateType isotope_lorentz_fwhm_;

    UInt charge_;
    CoordinateType mean_;
    CoordinateType monoisotopic_mz_;
    double averagine_[AVERAGINE_NUM];
    Int max_isotope_;
    double trim_right_cutoff_;
    double isotope_distance_;
    IsotopeDistribution isotope_distribution_;

    void updateMembers_() override;

  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
