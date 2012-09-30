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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLERALIGN_H
#define OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLERALIGN_H

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

namespace OpenMS
{

	/**
		@brief Linear Resampling of raw data with alignment.
		
		This class can be used to generate uniform data from non-uniform raw data (e.g. ESI-TOF or MALDI-TOF experiments).
		Therefore the intensity at every position x in the input raw data is spread to the two
		adjacent resampling points.
		This method preserves the area of the input signal and also the centroid position of a peak.
		Therefore it is recommended for quantitation as well as for ProteinIdentification experiments.

    In addition to the LinearResampler, this class also allows to fix the
    points at which resampling will occur. This is useful if the resampling
    points are known in advance, e.g. if one needs to resample a chromatogram
    at the positions of another chromatogram.
		
	*/
  class LinearResamplerAlign
    : public LinearResampler

  {

public:

	/** 
		@brief Applies the resampling algorithm to an MSSpectrum.
	*/
  template <template <typename> class MSSpectrum, typename PeakType>
  void raster(MSSpectrum<PeakType>& spectrum)
  {
    //return if nothing to do
    if (spectrum.empty()) return;
    
    typename MSSpectrum<PeakType>::iterator first = spectrum.begin();
    typename MSSpectrum<PeakType>::iterator last = spectrum.end();
    
    double end_pos = (last-1)->getMZ();
    double start_pos = first->getMZ();
    int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing_ + 1));

    typename std::vector<PeakType> resampled_peak_container;
    resampled_peak_container.resize(number_resampled_points);

    // generate the resampled peaks at positions origin+i*spacing_
    typename std::vector<PeakType>::iterator it = resampled_peak_container.begin();
    for (int i=0; i < number_resampled_points; ++i)
    {
      it->setMZ( start_pos + i*spacing_);
      ++it;
    }

    raster(spectrum.begin(), spectrum.end(), resampled_peak_container.begin(), resampled_peak_container.end());
    
    resampled_peak_container.swap(spectrum);
  }

	/** 
		@brief Applies the resampling algorithm to an MSSpectrum but it will be aligned between start_pos and end_pos
	*/
  template <template <typename> class MSSpectrum, typename PeakType>
  void raster_align(MSSpectrum<PeakType>& spectrum, double start_pos, double end_pos)
  {
    //return if nothing to do
    if (spectrum.empty()) return;
    if (end_pos < start_pos)
    {
      MSSpectrum<PeakType> empty;
      empty.swap(spectrum);
      return;
    }
    
    typename MSSpectrum<PeakType>::iterator first = spectrum.begin();
    typename MSSpectrum<PeakType>::iterator last = spectrum.end();

    // get the iterators just before / after the two points start_pos / end_pos
    while (first != spectrum.end() && (first)->getMZ() < start_pos) {first++;}
    while (last != first && (last-1)->getMZ() > end_pos) {last--;}

    int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing_ + 1));

    typename std::vector<PeakType> resampled_peak_container;
    resampled_peak_container.resize(number_resampled_points);

    // generate the resampled peaks at positions origin+i*spacing_
    typename std::vector<PeakType>::iterator it = resampled_peak_container.begin();
    for (int i=0; i < number_resampled_points; ++i)
    {
        it->setMZ( start_pos + i*spacing_);
        ++it;
    }

    raster(first, last, resampled_peak_container.begin(), resampled_peak_container.end());
    
    resampled_peak_container.swap(spectrum);
  }

	/** 
		@brief Applies the resampling algorithm to an MSSpectrum.
	*/
  template < typename PeakTypeIterator, typename ConstPeakTypeIterator>
  void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end)
  {
    PeakTypeIterator resample_start = resample_it;

    // need to get the raw iterator between two resampled iterators of the raw data
    while(raw_it->getMZ() < resample_it->getMZ() && raw_it != raw_end) 
    {
      resample_it->setIntensity( resample_it->getIntensity() + raw_it->getIntensity() );
      raw_it++;
    } 

    while(raw_it != raw_end)
    {
      //advance the resample iterator until our raw point is between two resampled iterators 
      while(resample_it->getMZ() < raw_it->getMZ() && resample_it != resample_end) {resample_it++;} 
      if (resample_it != resample_start) {resample_it--;}

      // if we have the last datapoint we break
      if ((resample_it+1) == resample_end) {break;}

      double dist_left =  fabs(raw_it->getMZ() - resample_it->getMZ());
      double dist_right = fabs(raw_it->getMZ() - (resample_it+1)->getMZ());

      // distribute the intensity of the raw point according to the distance to resample_it and resample_it+1
      resample_it->setIntensity(resample_it->getIntensity() + raw_it->getIntensity() * dist_right / (dist_left+dist_right));
      (resample_it+1)->setIntensity((resample_it+1)->getIntensity() + raw_it->getIntensity() * dist_left / (dist_left+dist_right));

      raw_it++;
    }

    // add the final intensity to the right
    while(raw_it != raw_end) 
    {
      resample_it->setIntensity( resample_it->getIntensity() + raw_it->getIntensity() );
      raw_it++;
    } 
  }

	/** 
		@brief Applies the resampling algorithm using a linear interpolation
	*/
  template < typename PeakTypeIterator>
  void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator it, PeakTypeIterator resampled_end)
  {
    PeakTypeIterator raw_start = raw_it;

    // need to get the resampled iterator between two iterators of the raw data
    while(it->getMZ() < raw_it->getMZ() && it != resampled_end) {it++;} 

    while(it != resampled_end)
    {
      //advance the raw_iterator until our current point we want to interpolate is between them
      while(raw_it->getMZ() < it->getMZ() && raw_it != raw_end) {raw_it++;} 
      if (raw_it != raw_start) {raw_it--;}

      // if we have the last datapoint we break
      if ((raw_it+1) == raw_end) {break;}

      // use a linear interpolation between raw_it and raw_it+1
      double m = ((raw_it+1)->getIntensity() - raw_it->getIntensity() ) / ((raw_it+1)->getMZ() - raw_it->getMZ() );
      it->setIntensity( raw_it->getIntensity() + (it->getMZ() - raw_it->getMZ())*m );
      it++;
    }

  }

  };

}

#endif
