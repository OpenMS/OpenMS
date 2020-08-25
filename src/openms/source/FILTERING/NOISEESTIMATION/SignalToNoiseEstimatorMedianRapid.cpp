// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>

#include <numeric>

// array_wrapper needs to be included before it is used
// only in boost1.64+. See issue #2790
#if OPENMS_BOOST_VERSION_MINOR >= 64
#include <boost/serialization/array_wrapper.hpp>
#endif
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace OpenMS
{

  void SignalToNoiseEstimatorMedianRapid::computeNoiseInWindows_(
      const std::vector<double>& mz_array, std::vector<double> int_array, 
      std::vector<double> & result, double mz_start)
  {
    // PRECONDITION
    assert(mz_array.size() == int_array.size());
    assert(mz_array.size() > 2);

    // compute mean and standard deviation
    double sum = std::accumulate(int_array.begin(), int_array.end(), 0.0);
    double int_mean = sum / int_array.size();
    double sq_sum = std::inner_product(int_array.begin(), int_array.end(), int_array.begin(), 0.0);
    double int_stdev = std::sqrt(sq_sum / int_array.size() - int_mean * int_mean);

    std::vector<double>::const_iterator mz_start_it = mz_array.begin();
    std::vector<double>::const_iterator mz_end_it;
    std::vector<double>::iterator int_start_win = int_array.begin();
    std::vector<double>::iterator int_end_win = int_array.begin();
    for (size_t i = 0; i < result.size(); i++)
    {
      // Compute the the correct windows in m/z
      double mz_end = mz_start + window_length_;
      mz_end_it = std::lower_bound(mz_start_it, (std::vector<double>::const_iterator)mz_array.end(), mz_end);

      // Compute the the correct windows in intensity
      std::iterator_traits< std::vector<double>::const_iterator >::difference_type iterator_pos = std::distance(mz_start_it, mz_end_it);
      std::advance(int_end_win, iterator_pos);

      // compute median of all data between intensity start and intensity end
      double median = computeMedian_(int_start_win, int_end_win);
      result[i] = median;

      // Deal with a median of zero
      //
      // If we find a zero here, try to impute some value that might make sense as noise value ...
      // alternatively, one could also remove all zeros and compute the median on that 
      if (result[i] == 0)
      {
        // Legacy implementation from SignalToNoiseEstimatorMedian
        //
        // max_intensity_ = gauss_global.mean + std::sqrt(gauss_global.variance) * auto_max_stdev_Factor_;
        // From the maximum intensity we can compute the value of the lowest
        // bin in the histogram of the SignalToNoiseEstimatorMedian algorithm:
        // maximum intensity divided by 60 
        result[i] = (int_mean  + 3.0 * int_stdev) / 60;
      }  

      mz_start_it = mz_end_it;
      int_start_win = int_end_win;
      mz_start += window_length_;
    }
  }

  double SignalToNoiseEstimatorMedianRapid::computeMedian_(std::vector<double>::iterator & first, std::vector<double>::iterator & last)
  {
    std::iterator_traits< std::vector<double>::const_iterator >::difference_type iterator_pos = std::distance(first, last);
    std::nth_element(first, first + iterator_pos / 2, last);

    double median;
    if (iterator_pos == 0)
    {
      median = 0.0;
    }
    else if (iterator_pos % 2 == 0)
    {
      // even case
      // compute the arithmetic mean between the two middle elements
      double f = *(first + iterator_pos / 2);
      std::nth_element(first, first + iterator_pos / 2 -1, last);
      double s = *(first + iterator_pos / 2 - 1);
      median = (f+s)/2.0; 
    }
    else
    {
      // odd case
      median = *(first + iterator_pos / 2);
    }
    return median;
  }

} // namespace OpenMS

