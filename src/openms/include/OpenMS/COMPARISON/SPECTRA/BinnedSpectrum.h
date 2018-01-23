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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>

namespace OpenMS
{

  /**
    @brief This is a binned representation of a PeakSpectrum

    @param sz the size of the bins and
    @param sp number of neighboring bins to both sides affected by a peak contribution
    @param ps the peakspectrum, that shall be represented

    sz denotes the size of a bin in @p Th, thereby deciding the number of bins(all of size sz) the spectrum is discretized to.
    Each bin will represent a certain @p Th range and the peaks will be put in the respective bins and sum up inside.
    sp denotes the number of neighboring bins to the left and the number of neighboring bins to the right a peak is also added to.
    E.g. a BinnedSpectrum with binsize of 0.5 @p Th will have a peak at 100 @p Th in bin no. 200, a peak at 100.1 @p Th will be in bin no. 201.
    If the binspread is 1, the peak at 100 Th will be added to bin no. 199, 200 and 201.
    If the binspread is 2, the peak at 100 @p Th will also be added to bin no. 198 and 202, and so on.

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSpectrum
  {

private:

    UInt bin_spread_;
    float bin_size_;
    /// The computed bins
    SparseVector<float> bins_;
    /// The original raw spectrum
    PeakSpectrum raw_spec_;

public:

    /**
        @brief  Exception which is thrown if BinnedSpectrum bins are accessed and no PeakSpektrum has been
                integrated yet i.e. bins_ is empty
    */
    class OPENMS_DLLAPI NoSpectrumIntegrated :
      public Exception::BaseException
    {
public:
      NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message = "BinnedSpectrum hasn't got a PeakSpectrum to base on yet") throw();

      ~NoSpectrumIntegrated() throw() override;
    };

    typedef SparseVector<float>::const_iterator const_bin_iterator;
    typedef SparseVector<float>::iterator bin_iterator;

    /// default constructor
    BinnedSpectrum();

    /// detailed constructor
    BinnedSpectrum(float size, UInt spread, PeakSpectrum ps);

    /// copy constructor
    BinnedSpectrum(const BinnedSpectrum& source);

    /// destructor
    virtual ~BinnedSpectrum();

    /// assignment operator
    BinnedSpectrum& operator=(const BinnedSpectrum& source)
    {
      if (&source != this)
      {
        setBinSize(source.getBinSize());
        setBinSpread(source.getBinSpread());
        bins_ = source.getBins();
        raw_spec_ = source.raw_spec_;
      }
      return *this;
    }

    /// assignment operator for PeakSpectra
    BinnedSpectrum& operator=(const PeakSpectrum& source)
    {
      if (raw_spec_ != source)
      {
        raw_spec_ = source;
        setBinning();
      }
      return *this;
    }

    /// equality operator
    bool operator==(const BinnedSpectrum& rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return raw_spec_ == rhs.raw_spec_ &&
             rhs.getBinSize() == this->bin_size_ &&
             rhs.getBinSpread() == this->bin_spread_;

#pragma clang diagnostic pop
    }

    /// inequality operator
    bool operator!=(const BinnedSpectrum& rhs) const
    {
      return !(operator==(rhs));
    }

    /// equality operator for PeakSpectra
    bool operator==(const PeakSpectrum& rhs) const
    {
      return raw_spec_ == rhs;
    }

    /// inequality operator for PeakSpectra
    bool operator!=(const PeakSpectrum& rhs) const
    {
      return !(operator==(rhs));
    }

    /// get the BinSize
    inline double getBinSize() const
    {
      return this->bin_size_;
    }

    /// get the BinSpread
    inline UInt getBinSpread() const
    {
      return this->bin_spread_;
    }

    /// get the BinNumber, number of Bins
    inline UInt getBinNumber() const
    {
      return (UInt) this->bins_.size();
    }

    /// get the FilledBinNumber, number of filled Bins
    inline UInt getFilledBinNumber() const
    {
      return (UInt) this->bins_.nonzero_size();
    }

    /** immutable access to the Bincontainer

            @throw NoSpectrumIntegrated is thrown if no spectrum was integrated
    */
    inline const SparseVector<float>& getBins() const
    {
      if (bins_.empty())
      {
        throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      return bins_;
    }

    /** mutable access to the Bincontainer

            @throw NoSpectrumIntegrated is thrown if no spectrum was integrated
    */
    inline SparseVector<float>& getBins()
    {
      if (bins_.empty())
      {
        try
        {
          this->setBinning();
        }
        catch (...)
        {
          throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
      }
      return bins_;
    }

    /// returns the const begin iterator of the container
    inline const_bin_iterator begin() const
    {
      return bins_.begin();
    }

    /// returns the const end iterator of the container
    inline const_bin_iterator end() const
    {
      return bins_.end();
    }

    /// returns the begin iterator of the container
    inline bin_iterator begin()
    {
      return bins_.begin();
    }

    /// returns the end iterator of the container
    inline bin_iterator end()
    {
      return bins_.end();
    }

    /** sets the BinSize_ (and re-bins)

            @param s defines the size of the bins
            @throw NoSpectrumIntegrated is thrown if no spectrum is integrated
    */
    inline void setBinSize(double s)
    {
      if (this->bin_size_ != s)
      {
        this->bin_size_ = s;
        try
        {
          this->setBinning();
        }
        catch (...)
        {
          throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
      }
    }

    /** sets the BinSpread_ (and re-bins)

            @param s defines the binning spread, given as positive integer
            @throw NoSpectrumIntegrated is thrown if no spec was integrated into the instance
    */
    inline void setBinSpread(UInt s)
    {
      if (this->bin_spread_ != s)
      {
        this->bin_spread_ = s;
        try
        {
          this->setBinning();
        }
        catch (...)
        {
          throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
      }
    }

    /** makes the binning: all Peaks of the containing PeakSpectrum are summed up in the bins corresponding to @p m/z ranges

            @throw NoSpectrumIntegrated is thrown if no spectrum was integrated before
    */
    void setBinning();

    /// function to check comparability of two BinnedSpectrum objects, i.e. if they have equal bin size and spread
    bool checkCompliance(const BinnedSpectrum& bs) const;

    /// Gives access to the underlying raw spectrum
    const PeakSpectrum& getRawSpectrum() const;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
