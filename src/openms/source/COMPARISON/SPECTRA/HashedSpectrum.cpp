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
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/COMPARISON/SPECTRA/HashedSpectrum.h>
#include <OpenMS/MATH/MISC/Spline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  HashedSpectrum::HashedSpectrum(MSSpectrum<Peak1D>& raw_spectrum, double mz_bin, double mz_tolerance, bool mz_unit_ppm)
  : mz_bin_(mz_bin), mz_tolerance_(mz_tolerance), mz_unit_ppm_(mz_unit_ppm)
  {
	if (raw_spectrum.size() <= 2)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, raw_spectrum.size());
    }
	
	mz_min_ = raw_spectrum.begin()->getMZ();
		    
    std::vector<double> mz;
    std::vector<double> intensity;
    int index = 0;
    int last_index = 0;
    HashedSpectrum::MzInterval bin;
    bin.first = &(*raw_spectrum.begin());
    bin.last = &(*raw_spectrum.begin());
    for (MSSpectrum<Peak1D>::Iterator it = raw_spectrum.begin(); it != raw_spectrum.end(); ++it)
    {		
		index = getIndex_(it->getMZ());
		
		if (index > last_index)
		{
			// close previous interval
			bins_.push_back(bin);
			
			// add empty interval
			bin.first = NULL;
			bin.last = NULL;
			for (int i=0; i < (index - last_index -1); ++i)
			{
				bins_.push_back(bin);
			}
		  
			// open new interval
			bin.first = &*it;
			bin.last = &*it;
			
			last_index = index;
		}
		else
		{
			bin.last = &*it;
		}
	  
      mz.push_back(it->getMZ());
      intensity.push_back(it->getIntensity());
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
  
  double HashedSpectrum::getMzTolerance() const
  {
    return mz_tolerance_;
  }
  
  bool HashedSpectrum::getMzUnitPpm() const
  {
    return mz_unit_ppm_;
  }
  
  int HashedSpectrum::getIndex_(double mz)
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
  
  double HashedSpectrum::getDistance_(double mz1, double mz2)
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
  
  MSSpectrum<Peak1D>::pointer HashedSpectrum::getPeak(double mz)
  {
	  // find range of bins in which possible peaks may lie
	  int index_lower_bound;
	  int index_upper_bound;
	  if (mz_unit_ppm_)
	  {
		// m/z tolerance in ppm
		index_lower_bound = getIndex_(mz * (1 - mz_tolerance_*1.0e-6));
		index_upper_bound = getIndex_(mz * (1 + mz_tolerance_*1.0e-6));
	  }
	  else
	  {
		// m/z tolerance in Da
		index_lower_bound = getIndex_(mz - mz_tolerance_);
		index_upper_bound = getIndex_(mz + mz_tolerance_);
	  }
	  
	  if (index_upper_bound < 0 || index_lower_bound >= (int) bins_.size())
	  {
		  return NULL;
	  }	  
	  
	  if (index_lower_bound < 0)
	  {
		  index_lower_bound = 0;		  
	  }
	  
	  if (index_upper_bound >= (int) bins_.size())
	  {
		  index_upper_bound = bins_.size() - 1;		  
	  }
	  	  
	  //std::cout << "index start = " << index_lower_bound << "    m/z start = " << bins_[index_lower_bound].first->getMZ() << "\n";
	  //std::cout << "index end   = " << index_upper_bound << "    m/z end   = " << bins_[index_upper_bound].last->getMZ() << "\n\n";
	  
	  MSSpectrum<Peak1D>::Iterator first_peak(bins_[index_lower_bound].first);
	  MSSpectrum<Peak1D>::Iterator last_peak(bins_[index_upper_bound].last);
	  
	  MSSpectrum<Peak1D>::Iterator it_mz_closest = first_peak;
	  double distance_closest = getDistance_(mz, first_peak->getMZ());
	  for (MSSpectrum<Peak1D>::Iterator it_mz = first_peak; it_mz <= last_peak; ++it_mz)
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
	  
	  if (distance_closest < mz_tolerance_)
	  {
		  return &(*it_mz_closest);
	  }
	  else
	  {
		  return NULL;
	  }
	  
  }  
  
}
