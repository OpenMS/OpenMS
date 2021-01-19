// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

namespace OpenMS
{
  /**
   * @brief Wrapper struct for all the structs needed by the FLASHDeconv
   *
   * @see FLASHDeconv
   * @reference: FeatureFinderAlgorithmPickedHelperStructs
   */

  struct OPENMS_DLLAPI FLASHDeconvHelperStructs
  {
    /// This struct contains the averagine patterns precalulated for speed up. Other variables are also calculated for fast cosine calculation
    struct OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes;
      /// L2 norms for masses
      std::vector<double> norms;
      /// mass differences between average mass and monoisotopic mass
      std::vector<double> average_mono_mass_difference;
      /// Isotope start indices: isotopes of the indices less than them have very low intensities
      std::vector<Size> isotope_start_indices;
      /// Isotope end indices: isotopes of the indices larger than them have very low intensities
      std::vector<Size> isotope_end_indices;
      /// max isotope index
      int max_isotope_index;

      /// mass interval for calculation
      double mass_interval;
      /// min mass for calculation
      double min_mass;
    public:
      /// default constructor
      PrecalculatedAveragine() = default;

      /**
       @brief constructor with parameters such as mass ranges and interval ( delta ).
       @param m min_mass_
       @param M max_mass
       @param delta mass interval between m and M
       @param generator gen
       \]erator by which the calulation is done
       @param useRNAavg if set, nucleotide patters are calculated
    */
      PrecalculatedAveragine(const double min_mass,
                             const double max_mass,
                             const double delta,
                             CoarseIsotopePatternGenerator *generator,
                             const bool use_RNA_averagine);

      /// get distribution for input mass
      IsotopeDistribution get(const double mass) const;

      /// get max isotope index
      int getMaxIsotopeIndex() const;

      /// get max isotope index
      void setMaxIsotopeIndex(const int index);

      /// get norm
      double getNorm(const double mass) const;

      /// get isotope start index
      Size getIsotopeStartIndex(const double mass) const;

      /// get isotope end index
      Size getIsotopeEndIndex(const double mass) const;

      /// get mass difference between avg and mono masses
      double getAverageMassDelta(const double mass) const;

    };

    /// log transformed peak. After deconvolution, other information such as charge and isotope index are stored.
    struct OPENMS_DLLAPI LogMzPeak
    {
      /// original peak mz
      double mz = 0;
      /// original peak intensity
      double intensity = 0;
      /// log transformed mz
      double logMz = -1000;
      /// deteremined mass after deconvolution. NOT monoisotopic but only decharged
      double mass = .0;
      /// charge
      int charge = 0;
      /// isotope index
      int isotopeIndex = -1;

      /// default constructor
      LogMzPeak() = default;

      /**
        //       @brief constructor from Peak1D.
        //       @param positive determines the charge carrier mass*/
      explicit LogMzPeak(const Peak1D& peak, const bool positive);

      /// copy constructor
      LogMzPeak(const LogMzPeak& ) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT monoisotopic mass. Valid only when charge is set
      double getUnchargedMass();

      bool operator<(const LogMzPeak& a) const;

      bool operator>(const LogMzPeak& a) const;

      bool operator==(const LogMzPeak& other) const;

    };

    /**
        //       @brief calculate averagines
        //       @param maxMass max mass
        //       @param useRNAavg if set, nucleotides averagines are calculated */
    static PrecalculatedAveragine calculateAveragines(const double max_mass, const bool use_RNA_averagine);

    /**
        //       @brief calculate log mzs from mzs
        //       @param mz mz
        //       @param positive determines the charge carrier mass
       */
    static double getLogMz(const double mz, const bool positive);

    /**
            //       @brief get charge carrier mass
            //       @param positive determines the charge carrier mass*/
    static double getChargeMass(const bool positive);
  };
}