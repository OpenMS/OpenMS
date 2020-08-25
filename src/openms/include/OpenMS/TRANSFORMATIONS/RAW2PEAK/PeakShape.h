// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#pragma once

#include <cmath>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /**
    @brief Internal representation of a peak shape (used by the PeakPickerCWT)

    It defines an asymmetric Lorentzian and asymmetric hyperbolic squared secan function.
  */
  class OPENMS_DLLAPI PeakShape
  {
  public:
    /**
      @brief Peak shape type (asymmetric Lorentzian or asymmetric hyperbolic secans squared).

      The peak shape can represent an asymmetric Lorentzian function, given by

      l(x) = height/(1.+pow(left_width*(x - mz_position), 2)) (x<=mz_position)

      l(x) = height/(1.+pow(right_width*(x - mz_position), 2)) (x>mz_position)

      or an asymmetric hyperbolic secans squared function

      s(x) = height/pow(cosh(left_width*(x-mz_position)), 2) (x<=mz_position)

      s(x) = height/pow(cosh(right_width*(x-mz_position)), 2) (x>mz_position)
    */
    enum Type
    {
      LORENTZ_PEAK,
      SECH_PEAK,
      UNDEFINED
    };

    /// Iterator to the raw data vector
    typedef MSSpectrum::const_iterator PeakIterator;

    /// Default constructor
    PeakShape();

    /// Constructor that sets most of the members
    PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, PeakIterator left_, PeakIterator right_, Type type_);

    /// Constructor that sets most of the members
    PeakShape(double height_, double mz_position_, double left_width_, double right_width_, double area_, Type type_);


    /// Copy constructor
    PeakShape(const PeakShape & rhs);

    /// Destructor
    virtual ~PeakShape();
    
    /// Assignment operator
    PeakShape & operator=(const PeakShape & rhs);

    //Equality operator
    bool operator==(const PeakShape & rhs) const;
    //Equality operator
    bool operator!=(const PeakShape & rhs) const;

    /// Compute the intensity of the peaks shape at position x
    double operator()(double x) const;
    /// Computes symmetry measure of the peak shape, which is corresponds to the ratio of the left and right width parameters.
    double getSymmetricMeasure() const;
    /// Estimates the full width at half maximum.
    double getFWHM() const;
    /// Check if endpoint iterators are provided
    bool iteratorsSet() const;

    PeakIterator getLeftEndpoint() const;
    void setLeftEndpoint(PeakIterator left_endpoint);

    PeakIterator getRightEndpoint() const;
    void setRightEndpoint(PeakIterator right_endpoint);
    /// Maximum intensity of the peak shape
    double height;
    /// Centroid position
    double mz_position;
    /// Left width parameter
    double left_width;
    /// Right width parameter
    double right_width;
    /// Area of the peak shape
    double area;
    /** @brief Correlation coefficient.

      It represents the squared Pearson correlation coefficient with the original data (0 <= r_value <= 1).
    */
    double r_value;
    /// The signal to noise ratio at the mz_position
    double signal_to_noise;

    ///peak shape type
    Type type;

    /**
         @brief  Comparison of mz_positions.
    */
    class OPENMS_DLLAPI PositionLess
    {
public:
      inline bool operator()(const PeakShape & a, const PeakShape & b)
      {
        return a.mz_position < b.mz_position;
      }

    };
protected:
    /// Left peak endpoint in the data
    PeakIterator left_endpoint_;
    /// Right peak endpoint in the data
    PeakIterator right_endpoint_;
    /// Needed for initialization of endpoint iterators
    MSSpectrum exp_;
    /// flag if left endpoint iterator differs from default value
    bool left_iterator_set_;
    /// flag if left endpoint iterator differs from default value
    bool right_iterator_set_;
  };
} // namespace OpenMS

