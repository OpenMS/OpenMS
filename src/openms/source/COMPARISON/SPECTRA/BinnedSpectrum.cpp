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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>


using namespace std;

namespace OpenMS
{
  BinnedSpectrum::BinnedSpectrum() :
    bin_spread_(1), bin_size_(2.0), bins_(), raw_spec_()
  {
  }

  BinnedSpectrum::BinnedSpectrum(float size, UInt spread, PeakSpectrum ps) :
    bin_spread_(spread), bin_size_(size), bins_(), raw_spec_(ps)
  {
    setBinning();
  }

  BinnedSpectrum::BinnedSpectrum(const BinnedSpectrum& source) :
    bin_spread_(source.getBinSpread()), bin_size_(source.getBinSize()), bins_(source.getBins()), raw_spec_(source.raw_spec_)
  {
  }

  BinnedSpectrum::~BinnedSpectrum()
  {
  }

  //accessors and operators see .h file


  void BinnedSpectrum::setBinning()
  {
    if (raw_spec_.empty())
    {
      throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    bins_.clear();

    //make all necessary bins accessible
    raw_spec_.sortByPosition();
    bins_ = SparseVector<float>((UInt)ceil(raw_spec_.back().getMZ() / bin_size_) + bin_spread_, 0, 0);

    //put all peaks into bins
    UInt bin_number;
    for (Size i = 0; i < raw_spec_.size(); ++i)
    {
      //bin_number counted form 0 -> floor
      bin_number = (UInt)floor(raw_spec_[i].getMZ() / bin_size_);
      //e.g. bin_size_ = 1.5: first bin covers range [0,1.5] so peak at 1.5 falls in first bin (index 0)

      if (raw_spec_[i].getMZ() / bin_size_ == (double)bin_number)
      {
        --bin_number;
      }

      //add peak to corresponding bin
      bins_[bin_number] = bins_.at(bin_number) + raw_spec_[i].getIntensity();

      //add peak to neighboring binspread many
      for (Size j = 0; j < bin_spread_; ++j)
      {
        bins_[bin_number + j + 1] = bins_.at(bin_number + j + 1) + raw_spec_[i].getIntensity();
        // we are not in one of the first bins (0 to bin_spread)
        //not working:  if (bin_number-j-1 >= 0)
        if (bin_number >= j + 1)
        {
          bins_[bin_number - j - 1] = bins_.at(bin_number - j - 1) + raw_spec_[i].getIntensity();
        }
      }
    }

  }

  //yields false if given BinnedSpectrum size or spread differs from this one (comparing those might crash)
  bool BinnedSpectrum::checkCompliance(const BinnedSpectrum& bs) const
  {
    return (this->bin_size_ == bs.getBinSize()) &&
           (this->bin_spread_ == bs.getBinSpread());
  }

  BinnedSpectrum::NoSpectrumIntegrated::NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message) throw() :
    BaseException(file, line, function, "BinnedSpectrum::NoSpectrumIntegrated", message)
  {
  }

  BinnedSpectrum::NoSpectrumIntegrated::~NoSpectrumIntegrated() throw()
  {
  }

  const PeakSpectrum& BinnedSpectrum::getRawSpectrum() const
  {
    return raw_spec_;
  }

}
