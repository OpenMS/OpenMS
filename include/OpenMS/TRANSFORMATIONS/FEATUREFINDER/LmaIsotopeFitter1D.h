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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>



namespace OpenMS
{
  /**
    @brief Isotope distribution fitter (1-dim.) approximated using Levenberg-Marquardt algorithm (GSL implementation) for parameter optimization.

    @htmlinclude OpenMS_LmaIsotopeFitter1D.parameters
   */
  class OPENMS_DLLAPI LmaIsotopeFitter1D :
    public LevMarqFitter1D
  {
public:

    enum Averagines {C = 0, H, N, O, S, AVERAGINE_NUM};

    /// Default constructor
    LmaIsotopeFitter1D();

    /// copy constructor
    LmaIsotopeFitter1D(const LmaIsotopeFitter1D & source);

    /// destructor
    virtual ~LmaIsotopeFitter1D();

    /// assignment operator
    virtual LmaIsotopeFitter1D & operator=(const LmaIsotopeFitter1D & source);

    /// create new LmaIsotopeFitter1D object (function needed by Factory)
    static Fitter1D * create()
    {
      return new LmaIsotopeFitter1D();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "LmaIsotopeFitter1D";
    }

    /// return interpolation model
    QualityType fit1d(const RawDataArrayType & range, InterpolationModel * & model);

protected:

    /// Helper struct (contains the size of an area, a raw data container, the relative abundance of i-th isotopic peak and the distance between consecutive isotopic peaks)
    struct Data
    {
      typedef Peak1D PeakType;
      typedef std::vector<PeakType> RawDataArrayType;
      typedef std::vector<double> ContainerType;
      typedef Feature::CoordinateType CoordinateType;

      Size n;
      RawDataArrayType set;
      ContainerType isotopes_exact;
      CoordinateType isotope_distance;
      // bool mono_known;
      // CoordinateType monoisotopic_mz;
      CoordinateType isotopes_stdev;
      CoordinateType sigma;
    };

    /// Compute start parameter
    void setInitialParameters_();

    /// Evaluation of the target function for nonlinear optimization
    static Int residual_(const gsl_vector * x, void * params, gsl_vector * f);

    /// Compute the Jacobian matrix, where each row of the matrix corresponds to a point in the data
    static Int jacobian_(const gsl_vector * x, void * params, gsl_matrix * J);

    /// Driver function for the evaluation of function and jacobian
    static Int evaluate_(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J);

    /** Display the intermediate state of the solution. The solver state contains
        the vector s->x which is the current position, and the vector s->f with
        corresponding function values
    */
    void printState_(Int iter, gsl_multifit_fdfsolver * s);

    /// isotope charge
    UInt charge_;
    /// standard derivation in isotope
    CoordinateType isotope_stdev_;
    /// total intensity (area under curve)
    CoordinateType total_intensity_;
    /// monoisotopic mass
    CoordinateType monoisotopic_mz_;
    /// maximum isotopic rank to be considered
    Int max_isotope_;
    /// cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered
    DoubleReal trim_right_cutoff_;
    /// distance between consecutive isotopic peaks
    DoubleReal isotope_distance_;
    /// Centroid m/z (as opposed to monoisotopic m/z)
    CoordinateType mean_;
    /// number of an atom per Dalton of mass
    DoubleReal averagine_[AVERAGINE_NUM];
    /// relative abundance of i-th isotopic peak
    ContainerType isotopes_exact_;
    /// The position of the monoisotopic mass is known(=1) or unknown(=0).
    bool monoisotopic_mass_known_;

    void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LMAISOTOPEFITTER1D_H
