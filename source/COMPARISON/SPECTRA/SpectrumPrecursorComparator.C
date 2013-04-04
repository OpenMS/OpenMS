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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumPrecursorComparator::SpectrumPrecursorComparator() :
    PeakSpectrumCompareFunctor()
  {
    setName(SpectrumPrecursorComparator::getProductName());
    defaults_.setValue("window", 2, "Allowed deviation between precursor peaks.");
    defaultsToParam_();
  }

  SpectrumPrecursorComparator::SpectrumPrecursorComparator(const SpectrumPrecursorComparator & source) :
    PeakSpectrumCompareFunctor(source)
  {
  }

  SpectrumPrecursorComparator::~SpectrumPrecursorComparator()
  {
  }

  SpectrumPrecursorComparator & SpectrumPrecursorComparator::operator=(const SpectrumPrecursorComparator & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double SpectrumPrecursorComparator::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double SpectrumPrecursorComparator::operator()(const PeakSpectrum & x, const PeakSpectrum & y) const
  {
    double window = (double)param_.getValue("window");

    DoubleReal mz1 = 0.0;
    if (!x.getPrecursors().empty())
      mz1 = x.getPrecursors()[0].getMZ();
    DoubleReal mz2 = 0.0;
    if (!y.getPrecursors().empty())
      mz2 = y.getPrecursors()[0].getMZ();

    if (fabs(mz1 - mz2) > window)
    {
      return 0;
    }

    return window - fabs(mz1 - mz2);
  }

}
