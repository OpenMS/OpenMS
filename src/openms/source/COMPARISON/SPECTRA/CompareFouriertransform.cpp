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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>

#include <cmath>
#include <OpenMS/MATH/gsl_wrapper.h>


namespace OpenMS
{
  CompareFouriertransform::CompareFouriertransform() :
    PeakSpectrumCompareFunctor()
  {
    setName(CompareFouriertransform::getProductName());
    defaults_.setValue("epsilon", 0.2, "defines the absolute error of the mass spectrometer");
    defaultsToParam_();
  }

  CompareFouriertransform::CompareFouriertransform(const CompareFouriertransform & source) :
    PeakSpectrumCompareFunctor(source)
  {
  }

  CompareFouriertransform::~CompareFouriertransform()
  {
  }

  CompareFouriertransform & CompareFouriertransform::operator=(const CompareFouriertransform & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double CompareFouriertransform::operator()(const PeakSpectrum &) const
  {
    return 0;
  }

  double CompareFouriertransform::operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const
  {
    const MSSpectrum<>::FloatDataArrays & temp1 = spec1.getFloatDataArrays();
    if (temp1.size() == 0)
    {

      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input need to be a fouriertransformation, try first transform ()");
      //transform(spec1);
    }

    UInt i = searchTransformation_(spec1);

    const MSSpectrum<>::FloatDataArrays & temp2 = spec2.getFloatDataArrays();
    if (temp2.size() == 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input need to be a fouriertransformation, try first transform ()");
    }
    UInt j = searchTransformation_(spec2);

    if (temp1[i].size() != temp2[j].size())
    {
      // std::cout<< temp1[i].size() << temp2[j].size() << std::endl;
      return 0.0;
    }
    else
    {

      Real sum = 0;
      for (Size k = 0; k < temp1[i].size(); ++k)
      {

        sum = sum + temp1[i][k] - temp2[j][k];
      }
      // std::cout << sum << " summe " << std::endl;
      if (sum != 0)
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }
  }

  UInt CompareFouriertransform::searchTransformation_(const PeakSpectrum & spec) const
  {
    const MSSpectrum<>::FloatDataArrays & temp = spec.getFloatDataArrays();
    UInt i = 0;
    while (i < temp.size())
    {
      if (temp[i].getName() == "Fouriertransformation")
      {
        break;
      }
      else
      {
        ++i;
      }
    }
    //find entry
    if (i == temp.size() || temp[i].getName() != "Fouriertransformation")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input need to be a fouriertransformation, try first transform ()");            //transform_(spec);
      //return i+1;

    }
    else
    {
      return i;
    }
  }

  void CompareFouriertransform::transform(PeakSpectrum & spec)
  {

    double * data = new double[2 * spec.size()];
    //normalize first the intensity!!!
    DoubleReal int_sum = 0;
    for (Size p = 0; p < spec.size(); ++p)
    {
      int_sum += spec[p].getIntensity();
    }
    //copy the peaks two times
    Size i = 0;
    for (Size p = 0; p < spec.size(); ++p)
    {
      data[i] = spec[p].getIntensity() / int_sum;
      ++i;
    }
    for (SignedSize p = spec.size() - 1; p >= 0; --p)
    {
      data[i] = spec[p].getIntensity() / int_sum;
      ++i;
    }

    deprecated_gsl_fft_real_wavetable * real;
    deprecated_gsl_fft_real_workspace * work;
    work = deprecated_gsl_fft_real_workspace_alloc(spec.size());
    real = deprecated_gsl_fft_real_wavetable_alloc(spec.size());
    deprecated_gsl_fft_real_transform(data, 1, spec.size(), real, work);
    deprecated_gsl_fft_real_wavetable_free(real);
    deprecated_gsl_fft_real_workspace_free(work);

    MSSpectrum<>::FloatDataArrays & temp = spec.getFloatDataArrays();
    i = temp.size();
    temp.resize(i + 1);
    temp[i].setName("Fouriertransformation");
    UInt j = 0;
    while (j < spec.size())
    {
      temp[i].push_back(data[j]);
      if (j == 0)
        ++j;          //we only intress in the real part of FFT
      else
        j = j + 2;
    }
    delete[] data;
  }

}
