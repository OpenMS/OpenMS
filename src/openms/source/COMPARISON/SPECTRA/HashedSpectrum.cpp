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
	
	std::cout << "number of data points   = " << raw_spectrum.size() << "\n\n";
	/*std::cout << "m/z bin size            = " << mz_bin_size_ << "\n\n";
	if (mz_bin_unit_ppm_)
	{
		std::cout << "m/z unit is ppm.\n";
	}
	else
	{
		std::cout << "m/z unit is Da.\n";
	}*/
	    
    std::vector<double> mz;
    std::vector<double> intensity;
    int index = 0;
    int last_index = 0;
    HashedSpectrum::MzInterval interval;
    interval.first = &*raw_spectrum.begin();
    interval.last = &*raw_spectrum.begin();
    for (MSSpectrum<Peak1D>::Iterator it = raw_spectrum.begin(); it != raw_spectrum.end(); ++it)
    {		
		index = getIndex_(it->getMZ());
		
		if (index > last_index)
		{
			// close previous interval
			intervals_.push_back(interval);
			
			// add empty intervals
			interval.first = NULL;
			interval.last = NULL;
			for (int i=0; i < (index - last_index -1); ++i)
			{
				intervals_.push_back(interval);
			}
		  
			// open new interval
			interval.first = &*it;
			interval.last = &*it;
			
			last_index = index;
		}
		else
		{
			interval.last = &*it;
		}
	  
	  //std::cout << "m/z = " << it->getMZ() << "    index = " << index << "\n";
		
      mz.push_back(it->getMZ());
      intensity.push_back(it->getIntensity());
    }
    // close last interval
    intervals_.push_back(interval);
    
    //std::cout << "\n";
    
    //std::cout << "number of intervals = " << intervals_.size() << "\n\n";
           
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
		// m/z bins in ppm
		return (int) floor(log(mz/mz_min_)/log(1.0 + mz_bin_*1.0e-6));
	}
	else
	{
		// m/z bins in Da
		return (int) floor((mz - mz_min_)/mz_bin_);
	}
  }
  
  MSSpectrum<Peak1D>::pointer HashedSpectrum::getPeak(double mz)
  {
	  std::cout << "m/z min = " << mz_min_ << "\n\n";
	  
	  int index = getIndex_(mz);
	  std::cout << "index = " << index << "\n";
	  
	  if (index < 0 || index >= (int) intervals_.size())
	  {
		  return NULL;
	  }
	  
	  std::cout << "m/z start = " << intervals_[index].first->getMZ() << "\n";
	  std::cout << "m/z end   = " << intervals_[index].last->getMZ() << "\n\n";
	  
	  return NULL;
  }  
  
}
