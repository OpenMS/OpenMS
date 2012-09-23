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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
#ifndef OPENMS_COMPARISON_SPECTRA_COMPAREFOURIERTRANSFORM_H
#define OPENMS_COMPARISON_SPECTRA_COMPAREFOURIERTRANSFORM_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{
  /**
      @brief Compare Discrete Cosines value from a Fourier transformation, also known as Discrete Cosines Transformation

      The Direct Cosines Transformation based on the theory of the Fourier transformation.
      In this class the Fast Fourier Transformation(FFT) algorithm of the gsl library is used. FFT has a run-time complexity of
      n (log n). To get the Direct Cosines Transformation from a FFT there is preparation necessary. First the
      input data has to be mirrored. This is necessary, because FFT needs data which has a periodic nature.
      After the computation of FFT only the cosine values are important and stored in the Meta Data Array. So an inverse transformation of these values
      to get the original spectrum is not available.
      The comparison is done between two Meta Data Arrays, which contain the stored cosine values of their individual spectrum.
      The advantage of this method is how the comparison works. There is only one sum which has to be count, no multiplication is needed.

      Attention: only use the compare function, if the Spectrum was transformed earlier, else an error is going to appear.
      Only use this method of transformation, if you are sure there exists enough free memory. This is a fast estimation, but it only gives one or
      zero back.

      @htmlinclude OpenMS_CompareFouriertransform.parameters

      @ingroup SpectraComparison

  */

  class OPENMS_DLLAPI CompareFouriertransform :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    CompareFouriertransform();

    /// copy constructor
    CompareFouriertransform(const CompareFouriertransform & source);

    /// destructor
    virtual ~CompareFouriertransform();
    // @}

    // @name Operators
    // @{
    /// assignment operator
    CompareFouriertransform & operator=(const CompareFouriertransform & source);

    /**
      @brief Dummy function

          This function only returns 0 for any given PeakSpectrum, please use the other compare operator function
    */
    double operator()(const PeakSpectrum &) const;
    /**
        @brief compare two PeakSpectrum by their Discrete Cosines Transformation.

        This function compares two given PeakSpectrum about their  Discrete Cosines Transformation.
        First, a transformation has to be calculated. Please use the function transform() in this class, before calling this
            function. The comparison works by summing the subtractions of each coefficient for all elements of both transformations. sum(_i=1)
            ^n x_i-y_i. If the sum is zero, both Spectrums are identical in the real part and one is emited, otherwise a zero.
        */
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const;

    ///
    static PeakSpectrumCompareFunctor * create() { return new CompareFouriertransform(); }

    ///Returns the name used in the factory
    static const String getProductName()
    {
      return "CompareFouriertransform";
    }

    /**
    @brief calculate the Discrete Cosines Fourier Transformation.

    This function transforms a given PeakSpectrum to a Discrete Cosines Fourier Transformation. It stores only the part of the cosines                  of the FFT in
    the FloatDataArray which is a container from the PeakSpectrum. Only call this function, if you are sure there is no other               transformation done earlier over the same PeakSpectrum, because it isn't checked if there already exists a transformation.
*/
    void transform(PeakSpectrum & spec);
protected:
    /**
        @brief Search in the PeakSpectrum, if a Discrete Fourier transformation occurs, if not an error is going to be thrown, else the index               of the occurrence is returned.

        This function gives back the position, where the transformation was saved in a FloatDataArray. If there is no entry, an error is thrown to indicate that a transformation has to be calculated before calling this comparison operator.
    */
    UInt searchTransformation_(const PeakSpectrum & spec) const;


  };

}
#endif /*OPENMS_COMPARISON_SPECTRA_COMPAREFOURIERTRANSFORM_H*/
