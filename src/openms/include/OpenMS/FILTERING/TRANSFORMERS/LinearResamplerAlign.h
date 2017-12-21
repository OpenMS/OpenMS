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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLERALIGN_H
#define OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLERALIGN_H

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/CONCEPT/Macros.h>

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
  class OPENMS_DLLAPI LinearResamplerAlign :
    public LinearResampler
  {

public:

    LinearResamplerAlign()
    {
      defaults_.setValue("spacing", 0.05, "Spacing of the resampled output peaks.");
      defaults_.setValue("ppm", "false", "Whether spacing is in ppm or Th");
      defaultsToParam_();
    }

    /**
        @brief Applies the resampling algorithm to a container (MSSpectrum or MSChromatogram).

        The container will be resampled at equally spaced points between the
        start and end of the container. The resampling frequency can be
        controlled by the "spacing" parameter.

        @param container The container to be resampled
    */
    template <class SpecT>
    void raster(SpecT& container)
    {
      //return if nothing to do
      if (container.empty()) return;

      typename SpecT::iterator first = container.begin();
      typename SpecT::iterator last = container.end();

      double end_pos = (last - 1)->getMZ();
      double start_pos = first->getMZ();
      int number_resampled_points = (int)(ceil((end_pos - start_pos) / spacing_ + 1));

      std::vector<typename SpecT::PeakType> resampled_peak_container;
      populate_raster_(resampled_peak_container, start_pos, end_pos, number_resampled_points);

      raster(container.begin(), container.end(), resampled_peak_container.begin(), resampled_peak_container.end());

      container.swap(resampled_peak_container);
    }

    /**
        @brief Applies the resampling algorithm to a container (MSSpectrum or MSChromatogram) with fixed coordinates.

        The container will be resampled at equally spaced points between the
        supplied start and end positions. The resampling frequency can be
        controlled by the "spacing" parameter.

        This allows the user to specify the grid for alignment explicitly.
        This is especially useful if multiple spectra or chromatograms need to
        be resampled according to the same raster.

        @param container The container to be resampled
        @param start_pos The start position to be used for resampling
        @param end_pos The end position to be used for resampling
    */
    template <typename SpecT>
    void raster_align(SpecT& container, double start_pos, double end_pos)
    {
      //return if nothing to do
      if (container.empty()) return;

      if (end_pos < start_pos)
      {
        std::vector<typename SpecT::PeakType> empty;
        container.swap(empty);
        return;
      }

      typename SpecT::iterator first = container.begin();
      typename SpecT::iterator last = container.end();

      // get the iterators just before / after the two points start_pos / end_pos
      while (first != container.end() && (first)->getMZ() < start_pos) {++first;}
      while (last != first && (last - 1)->getMZ() > end_pos) {--last;}

      int number_resampled_points = (int)(ceil((end_pos - start_pos) / spacing_ + 1));

      std::vector<typename SpecT::PeakType> resampled_peak_container;
      populate_raster_(resampled_peak_container, start_pos, end_pos, number_resampled_points);

      raster(first, last, resampled_peak_container.begin(), resampled_peak_container.end());

      container.swap(resampled_peak_container);
    }

    /**
        @brief Applies the resampling algorithm.

        This will use the raster provided by the output container to resample
        the data provided in the input container. The intensities will be added
        to the intensities in the output container (which in most cases will be
        zero).

        The intensities will be distributed between the two closest resampling
        points, thus conserving the sum of intensity over the whole container.

        Note that all intensity in the input data container that is in peaks
        outside the range of the output container will simply be added to the
        first or last data point.

        @param raw_it Start of the input container to be resampled (containing the data)
        @param raw_end End of the input container to be resampled (containing the data)
        @param resample_it Iterator pointing to start of the output spectrum range (m/z need to be populated, intensities should be zero)
        @param resample_it Iterator pointing to end of the output spectrum range (m/z need to be populated, intensities should be zero)

    */
    template <typename PeakTypeIterator, typename ConstPeakTypeIterator>
    void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end)
    {
      OPENMS_PRECONDITION(resample_it != resample_end, "Output iterators cannot be identical") // as we use +1
      // OPENMS_PRECONDITION(raw_it != raw_end, "Input iterators cannot be identical")

      PeakTypeIterator resample_start = resample_it;

      // need to get the raw iterator between two resampled iterators of the raw data
      while (raw_it != raw_end && raw_it->getMZ() < resample_it->getMZ())
      {
        resample_it->setIntensity(resample_it->getIntensity() + raw_it->getIntensity());
        raw_it++;
      }

      while (raw_it != raw_end)
      {
        //advance the resample iterator until our raw point is between two resampled iterators
        while (resample_it != resample_end && resample_it->getMZ() < raw_it->getMZ()) {resample_it++;}
        if (resample_it != resample_start) {resample_it--;}

        // if we have the last datapoint we break
        if ((resample_it + 1) == resample_end) {break;}

        double dist_left =  fabs(raw_it->getMZ() - resample_it->getMZ());
        double dist_right = fabs(raw_it->getMZ() - (resample_it + 1)->getMZ());

        // distribute the intensity of the raw point according to the distance to resample_it and resample_it+1
        resample_it->setIntensity(resample_it->getIntensity() + raw_it->getIntensity() * dist_right / (dist_left + dist_right));
        (resample_it + 1)->setIntensity((resample_it + 1)->getIntensity() + raw_it->getIntensity() * dist_left / (dist_left + dist_right));

        raw_it++;
      }

      // add the final intensity to the right
      while (raw_it != raw_end)
      {
        resample_it->setIntensity(resample_it->getIntensity() + raw_it->getIntensity());
        raw_it++;
      }
    }

    /**
        @brief Applies the resampling algorithm.

        This will use the raster provided by the output container to resample
        the data provided in the input container. The intensities will be added
        to the intensities in the output container (which in most cases will be
        zero).

        The intensities will be distributed between the two closest resampling
        points, thus conserving the sum of intensity over the whole container.

        Note that all intensity in the input data container that is in peaks
        outside the range of the output container will simply be added to the
        first or last data point.

        @param mz_raw_it Start of the input container to be resampled (containing the m/z data)
        @param mz_raw_end End of the input container to be resampled (containing the m/z data)
        @param int_raw_it Start of the input container to be resampled (containing the intensity data)
        @param int_raw_end End of the input container to be resampled (containing the intensity data)
        @param mz_resample_it Iterator pointing to start of the output spectrum range (m/z which need to be populated)
        @param mz_resample_it Iterator pointing to end of the output spectrum range (m/z which need to be populated)
        @param int_resample_it Iterator pointing to start of the output spectrum range (intensities)
        @param int_resample_it Iterator pointing to end of the output spectrum range (intensities)

    */
    template <typename PeakTypeIterator, typename ConstPeakTypeIterator>
#ifdef OPENMS_ASSERTIONS
    void raster(ConstPeakTypeIterator mz_raw_it, ConstPeakTypeIterator mz_raw_end,
        ConstPeakTypeIterator int_raw_it, ConstPeakTypeIterator int_raw_end,
        PeakTypeIterator mz_resample_it, PeakTypeIterator mz_resample_end,
        PeakTypeIterator int_resample_it, PeakTypeIterator int_resample_end)
#else
    void raster(ConstPeakTypeIterator mz_raw_it, ConstPeakTypeIterator mz_raw_end,
        ConstPeakTypeIterator int_raw_it, ConstPeakTypeIterator /* int_raw_end */,
        PeakTypeIterator mz_resample_it, PeakTypeIterator mz_resample_end,
        PeakTypeIterator int_resample_it, PeakTypeIterator /* int_resample_end */)
#endif
    {
      OPENMS_PRECONDITION(mz_resample_it != mz_resample_end, "Output iterators cannot be identical") // as we use +1
      OPENMS_PRECONDITION(std::distance(mz_resample_it, mz_resample_end) == std::distance(int_resample_it, int_resample_end),
          "Resample m/z and intensity iterators need to cover the same distance")
      OPENMS_PRECONDITION(std::distance(mz_raw_it, mz_raw_end) == std::distance(int_raw_it, int_raw_end),
          "Raw m/z and intensity iterators need to cover the same distance")
      // OPENMS_PRECONDITION(raw_it != raw_end, "Input iterators cannot be identical")

      PeakTypeIterator mz_resample_start = mz_resample_it;

      // need to get the raw iterator between two resampled iterators of the raw data
      while (mz_raw_it != mz_raw_end && (*mz_raw_it) < (*mz_resample_it) )
      {
        (*int_resample_it) = *int_resample_it + *int_raw_it;
        ++mz_raw_it;
        ++int_raw_it;
      }

      while (mz_raw_it != mz_raw_end)
      {
        //advance the resample iterator until our raw point is between two resampled iterators
        while (mz_resample_it != mz_resample_end && *mz_resample_it < *mz_raw_it)
        {
          ++mz_resample_it; ++int_resample_it;
        }
        if (mz_resample_it != mz_resample_start)
        {
          --mz_resample_it; --int_resample_it;
        }

        // if we have the last datapoint we break
        if ((mz_resample_it + 1) == mz_resample_end) {break;}

        double dist_left =  fabs(*mz_raw_it - *mz_resample_it);
        double dist_right = fabs(*mz_raw_it - *(mz_resample_it + 1));

        // distribute the intensity of the raw point according to the distance to resample_it and resample_it+1
        *(int_resample_it) = *int_resample_it + (*int_raw_it) * dist_right / (dist_left + dist_right);
        *(int_resample_it + 1) = *(int_resample_it + 1) + (*int_raw_it) * dist_left / (dist_left + dist_right);

        ++mz_raw_it;
        ++int_raw_it;
      }

      // add the final intensity to the right
      while (mz_raw_it != mz_raw_end)
      {
        *int_resample_it = *int_resample_it + (*int_raw_it);
        ++mz_raw_it;
        ++int_raw_it;
      }
    }

    /**
        @brief Applies the resampling algorithm using a linear interpolation

        This will use the raster provided by the output container to resample
        the data provided in the input container. The intensities will be added
        to the intensities in the output container (which in most cases will be
        zero).

        The intensities at the resampling point is computed by a linear
        interpolation between the two closest resampling points.

        @param raw_it Start of the input (raw) spectrum to be resampled
        @param raw_end End of the input (raw) spectrum to be resampled
        @param resample_it Iterator pointing to start of the output spectrum range (m/z need to be populated, intensities should be zero)
        @param resample_it Iterator pointing to end of the output spectrum range (m/z need to be populated, intensities should be zero)
    */
    template <typename PeakTypeIterator>
    void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resampled_end)
    {
      // OPENMS_PRECONDITION(resample_it != resampled_end, "Output iterators cannot be identical")
      OPENMS_PRECONDITION(raw_it != raw_end, "Input iterators cannot be identical") // as we use +1

      PeakTypeIterator raw_start = raw_it;

      // need to get the resampled iterator between two iterators of the raw data
      while (resample_it != resampled_end && resample_it->getMZ() < raw_it->getMZ()) {resample_it++;}

      while (resample_it != resampled_end)
      {
        //advance the raw_iterator until our current point we want to interpolate is between them
        while (raw_it != raw_end && raw_it->getMZ() < resample_it->getMZ()) {raw_it++;}
        if (raw_it != raw_start) {raw_it--;}

        // if we have the last datapoint we break
        if ((raw_it + 1) == raw_end) {break;}

        // use a linear interpolation between raw_it and raw_it+1
        double m = ((raw_it + 1)->getIntensity() - raw_it->getIntensity()) / ((raw_it + 1)->getMZ() - raw_it->getMZ());
        resample_it->setIntensity(raw_it->getIntensity() + (resample_it->getMZ() - raw_it->getMZ()) * m);
        resample_it++;
      }

    }

protected:

    /// Spacing of the resampled data
    bool ppm_;

    void updateMembers_() override
    {
      spacing_ =  param_.getValue("spacing");
      ppm_ =  (bool)param_.getValue("ppm").toBool();
    }

    /// Generate raster for resampled peak container
    template <typename PeakType>
    void populate_raster_(std::vector<PeakType>& resampled_peak_container,
        double start_pos, double end_pos, int number_resampled_points)
    {
      if (!ppm_)
      {
        // generate the resampled peaks at positions origin+i*spacing_
        resampled_peak_container.resize(number_resampled_points);
        typename std::vector<PeakType>::iterator it = resampled_peak_container.begin();
        for (int i = 0; i < number_resampled_points; ++i)
        {
          it->setMZ(start_pos + i * spacing_);
          ++it;
        }
      }
      else
      {
        // generate resampled peaks with ppm distance (not fixed)
        double current_mz = start_pos;
        while (current_mz < end_pos)
        {
          PeakType p;
          p.setIntensity(0);
          p.setMZ(current_mz);
          resampled_peak_container.push_back(p);

          // increment current_mz
          current_mz += current_mz * (spacing_ / 1e6);
        }
      }
    }
  };

}

#endif // OPENMS_FILTERING_TRANSFORMERS_LINEARRESAMPLERALIGN_H

