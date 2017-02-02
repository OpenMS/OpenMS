// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/HashedSpectrum.h>

#include <vector>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  HashedSpectrum::HashedSpectrum(const MSSpectrum<Peak1D>& spectrum, const double mz_bin, const bool mz_unit_ppm)
  : mz_bin_(mz_bin), mz_unit_ppm_(mz_unit_ppm), spectrum_(spectrum)
  {
    if (spectrum_.size() <= 2)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, spectrum_.size());
    }
 
    mz_min_ = spectrum_.begin()->getMZ();
      
    int last_index = 0;
    HashedSpectrum::MzInterval bin;
    bin.begin = spectrum_.begin();
    bin.end = spectrum_.begin();
    for (MSSpectrum<Peak1D>::ConstIterator it = spectrum_.begin(); it != spectrum_.end(); ++it)
    {  
      int index = getIndex_(it->getMZ());
  
      if (index > last_index)
      {
        // close previous interval
        bins_.push_back(bin);
   
        // add empty interval
        bin.begin = spectrum_.end();
        bin.end = spectrum_.end();
        for (int i=0; i < (index - last_index -1); ++i)
        {
          bins_.push_back(bin);
        }
    
        // open new interval
        bin.begin = it;
        bin.end = it;

        last_index = index;
      }
      else
      {
        bin.end = it;
      }   
    }
    
    // close last interval
    bins_.push_back(bin);              
  }

  HashedSpectrum::~HashedSpectrum()
  {
  }

  double HashedSpectrum::getMzBin() const
  {
    return mz_bin_;
  }
  
  bool HashedSpectrum::getMzUnitPpm() const
  {
    return mz_unit_ppm_;
  }
  
  int HashedSpectrum::getIndex_(const double mz) const
  {
    if (mz_unit_ppm_)
    {
      // m/z intervals in ppm
      return (int) floor(log(mz/mz_min_)/log(1.0 + mz_bin_*1.0e-6));
    }
    else
    {
      // m/z intervals in Da
     return (int) floor((mz - mz_min_)/mz_bin_);
    }
  }
  
  double HashedSpectrum::getDistance_(const double mz1, const double mz2) const
  {
    if (mz_unit_ppm_)
    {
      // m/z tolerance in ppm
     return std::abs(mz1 - mz2)/mz1*1.0e+6;
    }
    else
    {
      // m/z tolerance in Da
      return std::abs(mz1 - mz2);
    }
  }
  
  MSSpectrum<Peak1D>::ConstIterator HashedSpectrum::findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const
  {
    // find range of bins in which possible peaks may lie
    int index_lower_bound;
    int index_upper_bound;
    if (mz_unit_ppm)
    {
      // m/z tolerance in ppm
      index_lower_bound = getIndex_(mz * (1 - mz_tolerance*1.0e-6));
      index_upper_bound = getIndex_(mz * (1 + mz_tolerance*1.0e-6));
    }
    else
    {
      // m/z tolerance in Da
      index_lower_bound = getIndex_(mz - mz_tolerance);
      index_upper_bound = getIndex_(mz + mz_tolerance);
    }
   
    if (index_upper_bound < 0 || index_lower_bound >= (int) bins_.size())
    {
      return spectrum_.end();
    }   
   
    if (index_lower_bound < 0)
    {
      index_lower_bound = 0;    
    }
   
    if (index_upper_bound >= (int) bins_.size())
    {
      index_upper_bound = bins_.size() - 1;    
    }
      
    //std::cout << "index start = " << index_lower_bound << "\n";
    //std::cout << "index end   = " << index_upper_bound << "\n\n";
   
    // find first peak (i.e. the first non-empty bin)
    MSSpectrum<Peak1D>::ConstIterator begin = spectrum_.end();
    for (int index = index_lower_bound; index <= index_upper_bound; ++index)
    {
      if (bins_[index].begin != spectrum_.end())
      {
        begin = bins_[index].begin;
        break;
      }
    }
    if (begin == spectrum_.end())
    {
      return spectrum_.end();
    }
    MSSpectrum<Peak1D>::ConstIterator first_peak(begin);
   
    // find last peak (i.e. the last non-empty bin)
    MSSpectrum<Peak1D>::ConstIterator end = spectrum_.end();
    for (int index = index_upper_bound; index >= index_lower_bound; --index)
    {
      if (bins_[index].end != spectrum_.end())
      {
        end = bins_[index].end;
        break;
      }
    }
    if (end == spectrum_.end())
    {
      return spectrum_.end();
    }
    MSSpectrum<Peak1D>::ConstIterator last_peak(end);
   
    //std::cout << "m/z start = " << first_peak->getMZ() << "\n";
    //std::cout << "m/z end   = " << last_peak->getMZ() << "\n\n";
   
    // loop over relevant peaks
    MSSpectrum<Peak1D>::ConstIterator it_mz_closest = first_peak;
    double distance_closest = getDistance_(mz, first_peak->getMZ());
    for (MSSpectrum<Peak1D>::ConstIterator it_mz = first_peak; it_mz <= last_peak; ++it_mz)
    {    
      double distance = getDistance_(mz, it_mz->getMZ());
    
      if (distance <= distance_closest)
      {
        it_mz_closest = it_mz;
        distance_closest = distance;
      }
      else
      {
        break;
      }
    }
   
    return it_mz_closest;   
  }  
  
}
