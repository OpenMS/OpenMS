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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
  /**
    @brief Isotope distribution fitter (1-dim.) approximated using linear interpolation.

    @htmlinclude OpenMS_IsotopeFitter1D.parameters
   */
  class OPENMS_DLLAPI IsotopeFitter1D :
    public MaxLikeliFitter1D
  {
public:

    /// Default constructor
    IsotopeFitter1D();

    /// copy constructor
    IsotopeFitter1D(const IsotopeFitter1D & source);

    /// destructor
    ~IsotopeFitter1D() override;

    /// assignment operator
    virtual IsotopeFitter1D & operator=(const IsotopeFitter1D & source);

    /// create new IsotopeFitter1D object (function needed by Factory)
    static Fitter1D * create()
    {
      return new IsotopeFitter1D();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "IsotopeFitter1D";
    }

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, InterpolationModel * & model) override;

protected:

    /// isotope charge
    CoordinateType charge_;
    /// standard derivation in isotope
    CoordinateType isotope_stdev_;
    /// maximum isotopic rank to be considered
    Int max_isotope_;

    void updateMembers_() override;
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEFITTER1D_H
