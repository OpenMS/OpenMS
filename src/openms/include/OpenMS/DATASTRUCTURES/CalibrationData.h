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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H
#define OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RichPeak2D.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <algorithm>

namespace OpenMS
{

    /**
      @brief A helper class, holding all calibration points.

      Calibration points can be filled from Peptide IDs (using FeatureMaps or vector<PeptideIds>)
      or from Raw data (MSExperiment).

      From this data, a calibration function can be computed.

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

      CalDataType::CoordinateType getMZ(Size i) const;

      CalDataType::CoordinateType getRT(Size i) const;
      CalDataType::CoordinateType getIntensity(Size i) const;

      const_iterator begin() const;

      const_iterator end() const;

      Size size() const;

      bool empty() const;

      void clear();

      void setUsePPM(bool usePPM);

      bool usePPM() const;

      void insertCalibrationPoint(CalDataType::CoordinateType rt, CalDataType::CoordinateType mz_obs, CalDataType::IntensityType intensity, 
                                  CalDataType::CoordinateType mz_ref, double weight,
                                  int group = -1);

      Size getNrOfGroups() const;

      CalDataType::CoordinateType getError(Size i) const;

      CalDataType::CoordinateType getRefMZ(Size i) const;

      CalDataType::CoordinateType getWeight(Size i) const;

      int getGroup(Size i) const;

      static StringList getMetaValues();

      /**
        @brief Compute the median in the given RT range for every peak group
      */
      CalibrationData median(double rt_left, double rt_right) const;

      /**
        @brief Sort calibration points by RT, to allow for valid RT chunking
      */
      void sortByRT();


    private:
      MSSpectrum<RichPeak2D> data_; //< calibration points
      bool use_ppm_; //< return ppm values as y-values for the model instead of absolute delta in [Th]
      std::set<int> groups_; //< peak groups present in this data
    };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CALIBRATIONDATA_H