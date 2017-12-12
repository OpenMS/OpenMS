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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H
#define OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RichPeak2D.h>

#include <vector>
#include <set>

namespace OpenMS
{

    /**
      @brief A helper class, holding all calibration points.

      Calibration points can be filled from Peptide IDs (using FeatureMaps or vector<PeptideIds>)
      or from lock masses in Raw data (MSExperiment).

      The m/z error can be queried using getError(). The unit of error is either ppm or Th, depending on
      usePPM().

      Each calibration point can be assigned to a peak group. This should be done for calibration points
      derived from lock masses, to enable querying for a medianized representation of
      a lock mass trace in a certain RT range (see median()). For calibration points derived from 
      peptide IDs, this does not make sense.

      From this data, a calibration function can be computed (see MZTrafoModel class).

    */
    class OPENMS_DLLAPI CalibrationData
    {
    public:

      typedef RichPeak2D CalDataType;
      typedef std::vector<CalDataType>::const_iterator const_iterator;
      typedef std::vector<CalDataType>::value_type value_type;

      /**
        @brief Default constructor
      */
      CalibrationData();

      /**
        @brief Retrieve the observed m/z of the i'th calibration point
      */
      CalDataType::CoordinateType getMZ(Size i) const;

      /**
        @brief Retrieve the observed RT of the i'th calibration point
      */
      CalDataType::CoordinateType getRT(Size i) const;

      /**
        @brief Retrieve the intensity of the i'th calibration point
      */
      CalDataType::CoordinateType getIntensity(Size i) const;

      /**
        @brief Begin iterator for calibration points
      */
      const_iterator begin() const;
      
      /**
        @brief Past-the-end iterator for calibration points
      */
      const_iterator end() const;

      /**
        @brief Number of calibration points
      */
      Size size() const;

      /**
        @brief Do we have any calibration points
      */
      bool empty() const;

      /**
        @brief Remove all calibration points
      */
      void clear();

      /**
        @brief When calling getError(), should ppm error or m/z error be returned?
      */
      void setUsePPM(bool usePPM);

      /**
        @brief Current error unit (ppm or Th)
      */
      bool usePPM() const;

      /**
        @brief Add a new calibration point

        @param rt Retention time
        @param mz_obs Observed m/z
        @param intensity Intensity (useful for weighted model fitting)
        @param mz_ref Theoretical m/z
        @param weight Weight of calibration point (useful for weighted model fitting)
        @param group Peak group of this calibration point. Using -1 will not assign any peak group. See also: median()
      */
      void insertCalibrationPoint(CalDataType::CoordinateType rt, CalDataType::CoordinateType mz_obs, CalDataType::IntensityType intensity, 
                                  CalDataType::CoordinateType mz_ref, double weight,
                                  int group = -1);

      /**
        @brief Number of peak groups (can be 0).
      */
      Size getNrOfGroups() const;

      /**
        @brief Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())
      */
      CalDataType::CoordinateType getError(Size i) const;

      /**
        @brief Retrieve the theoretical m/z of the i'th calibration point
      */
      CalDataType::CoordinateType getRefMZ(Size i) const;

      /**
        @brief Retrieve the weight of the i'th calibration point
      */
      CalDataType::CoordinateType getWeight(Size i) const;

      /**
        @brief Retrieve the group of the i'th calibration point.

        @param i Index
        @return Group; returns -1 if peak has no group.
      */
      int getGroup(Size i) const;

      /**
        @brief List of meta-values which are used internally (for conversion to PeakMap).
      */
      static StringList getMetaValues();

      /**
        @brief Compute the median in the given RT range for every peak group

        This is usually applied on calibration data obtained from lock masses,
        where each lock mass has its own peak group.
        Median() then computes an 'medianized' observed(!) lock mass within a certain RT range
        and returns calibration data with one calibration point
        per group. Also intensity is 'medianized'. The theoretical m/z is expected to be
        identical for all calibration points in a peak group.

        Groups must be specified during insertCalibrationPoint().
        If no groups are present, the result is empty.

        The container must be sorted by RT (see sortByRT())!

        @param rt_left Left border of RT range to medianize
        @param rt_right Right border of RT range to medianize
        @return New container, containing median representation for each peak group

      */
      CalibrationData median(double rt_left, double rt_right) const;

      /**
        @brief Sort calibration points by RT, to allow for valid RT chunking
      */
      void sortByRT();


    private:
      std::vector<RichPeak2D> data_; ///< calibration points
      bool use_ppm_; ///< return ppm values as y-values for the model instead of absolute delta in [Th]
      std::set<int> groups_; ///< peak groups present in this data
    };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H
