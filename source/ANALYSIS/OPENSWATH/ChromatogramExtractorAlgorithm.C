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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  void ChromatogramExtractorAlgorithm::extract_value_tophat(
      const std::vector<double>::const_iterator& mz_start, 
            std::vector<double>::const_iterator& mz_it,
      const std::vector<double>::const_iterator& mz_end,
            std::vector<double>::const_iterator& int_it,
      const double& mz, double& integrated_intensity, double& mz_extraction_window, bool ppm)
  {
    integrated_intensity = 0;
    if (mz_start == mz_end)
    {
      return;
    }

    // calculate extraction window
    double left, right;
    if (ppm)
    {
      left  = mz - mz * mz_extraction_window / 2.0 * 1.0e-6;
      right = mz + mz * mz_extraction_window / 2.0 * 1.0e-6;
    }
    else
    {
      left  = mz - mz_extraction_window / 2.0;
      right = mz + mz_extraction_window / 2.0;
    }

    std::vector<double>::const_iterator mz_walker;
    std::vector<double>::const_iterator int_walker;

    // advance the mz / int iterator until we hit the m/z value of the next transition
    while (mz_it != mz_end && (*mz_it) < mz)
    {
      mz_it++; int_it++;
    }

    // walk right and left and add to our intensity
    mz_walker  = mz_it;
    int_walker = int_it;

    // if we moved past the end of the spectrum, we need to try the last peak of the spectrum (it could still be within the window)
    if (mz_it == mz_end)
    {
      mz_walker--; int_walker--;
    }

    // add the current peak if it is between right and left
    if ((*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker);
    }

    // walk to the left until we go outside the window, then start walking to the right until we are outside the window
    mz_walker  = mz_it;
    int_walker = int_it;
    if (mz_it != mz_start)
    {
      mz_walker--;
      int_walker--;
    }
    while (mz_walker != mz_start && (*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker); mz_walker--; int_walker--;
    }
    mz_walker  = mz_it;
    int_walker = int_it;
    if (mz_it != mz_end)
    {
      mz_walker++;
      int_walker++;
    }
    while (mz_walker != mz_end && (*mz_walker) > left && (*mz_walker) < right)
    {
      integrated_intensity += (*int_walker); mz_walker++; int_walker++;
    }
  }

  void ChromatogramExtractorAlgorithm::extractChromatograms(const OpenSwath::SpectrumAccessPtr input,
      std::vector< OpenSwath::ChromatogramPtr >& output, 
      std::vector<ExtractionCoordinates> extraction_coordinates, double& mz_extraction_window,
      bool ppm, double rt_extraction_window, String filter)
  {
    Size input_size = input->getNrSpectra();
    if (input_size < 1)
    {
      return;
    }

    if (output.size() != extraction_coordinates.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "Output and extraction coordinates need to have the same size");
    }

    int used_filter = get_filter_nr(filter);
    // assert that they are sorted!
    if (std::adjacent_find(extraction_coordinates.begin(), extraction_coordinates.end(), 
          ExtractionCoordinates::SortExtractionCoordinatesReverseByMZ) != extraction_coordinates.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "Input to extractChromatogram needs to be sorted by m/z");
    }

    //go through all spectra
    startProgress(0, input_size, "Extracting chromatograms");
    for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
    {
      setProgress(scan_idx);

      OpenSwath::SpectrumPtr sptr = input->getSpectrumById(scan_idx);
      OpenSwath::SpectrumMeta s_meta = input->getSpectrumMetaById(scan_idx);

      OpenSwath::BinaryDataArrayPtr mz_arr = sptr->getMZArray();
      OpenSwath::BinaryDataArrayPtr int_arr = sptr->getIntensityArray();
      std::vector<double>::const_iterator mz_start = mz_arr->data.begin();
      std::vector<double>::const_iterator mz_end = mz_arr->data.end();
      std::vector<double>::const_iterator mz_it = mz_arr->data.begin();
      std::vector<double>::const_iterator int_it = int_arr->data.begin();

      if (sptr->getMZArray()->data.size() == 0)
        continue;

      // go through all transitions / chromatograms which are sorted by
      // ProductMZ. We can use this to step through the spectrum and at the
      // same time step through the transitions. We increase the peak counter
      // until we hit the next transition and then extract the signal.
      for (Size k = 0; k < extraction_coordinates.size(); ++k)
      {
        double integrated_intensity = 0;
        double current_rt = s_meta.RT;
        if (rt_extraction_window > 0 && 
            (current_rt < extraction_coordinates[k].rt - rt_extraction_window / 2.0 || 
             current_rt > extraction_coordinates[k].rt + rt_extraction_window / 2.0) )
        {
          continue;
        }

        if (used_filter == 1)
        {
          extract_value_tophat( mz_start, mz_it, mz_end, int_it,
                  extraction_coordinates[k].mz, integrated_intensity, mz_extraction_window, ppm);
        }
        else if (used_filter == 2)
        {
          throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        }

        // Time is first, intensity is second
        output[k]->binaryDataArrayPtrs[0]->data.push_back(current_rt);
        output[k]->binaryDataArrayPtrs[1]->data.push_back(integrated_intensity);
      }
    }
    endProgress();
  }

  int ChromatogramExtractorAlgorithm::get_filter_nr(String filter)
  {
    if (filter == "tophat")
    {
      return 1;
    }
    else if (filter == "bartlett")
    {
      return 2;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Filter either needs to be tophat or bartlett");
    }
  }

}
