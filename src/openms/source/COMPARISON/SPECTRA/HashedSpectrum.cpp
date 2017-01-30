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

  HashedSpectrum::HashedSpectrum(MSSpectrum<Peak1D>& raw_spectrum, double mz_bin_size, bool mz_bin_unit_ppm)
  : mz_bin_size_(mz_bin_size), mz_bin_unit_ppm_(mz_bin_unit_ppm)
  {
	if (raw_spectrum.size() <= 2)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, raw_spectrum.size());
    }
	  
	std::cout << "number of data points   = " << raw_spectrum.size() << "\n\n";
	//std::cout << "m/z bin size            = " << mz_bin_size_ << "\n\n";
	//std::cout << "m/z bin size unit ppm   = " << mz_bin_unit_ppm_ << "\n\n";
	    
    std::vector<double> mz;
    std::vector<double> intensity;
    int index = 0;
    int last_index = 0;
    HashedSpectrum::MzInterval interval;
    interval.first = &*raw_spectrum.begin();
    interval.last = &*raw_spectrum.begin();
    for (MSSpectrum<Peak1D>::Iterator it = raw_spectrum.begin(); it != raw_spectrum.end(); ++it)
    {
		index = (int) floor((it->getMZ() - raw_spectrum.begin()->getMZ())/mz_bin_size_);			
		
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
			
			//std::cout << "\n";
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
    
    mz_min_ = mz.front();
    mz_max_ = mz.back();
           
  }

  HashedSpectrum::~HashedSpectrum()
  {
  }

  double HashedSpectrum::getMzMin() const
  {
    return mz_min_;
  }

  double HashedSpectrum::getMzMax() const
  {
    return mz_max_;
  }

  double HashedSpectrum::getMzBinSize() const
  {
    return mz_bin_size_;
  }
  
  bool HashedSpectrum::getMzBinUnitPpm() const
  {
    return mz_bin_unit_ppm_;
  }
  
}
