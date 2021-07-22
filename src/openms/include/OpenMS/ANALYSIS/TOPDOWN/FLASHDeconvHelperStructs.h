// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
    /// This struct contains the averagine patterns pre-calulated for speed up. Other variables are also calculated for fast cosine calculation
    struct OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes_;
      /// L2 norms_ for masses
      std::vector<double> norms_;
      /// mass differences between average mass and monoisotopic mass
      std::vector<double> average_mono_mass_difference_;
      /// Isotope start indices: isotopes of the indices less than them have very low intensities
      std::vector<int> left_count_from_apex_;
      /// Isotope end indices: isotopes of the indices larger than them have very low intensities
      std::vector<int> right_count_from_apex_;
      /// most abundant isotope index
      std::vector<Size> apex_index_;

      /// max isotope index
      Size max_isotope_index_;
      /// mass interval for calculation
      double mass_interval_;
      /// min mass for calculation
      double min_mass_;
    public:
      /// default constructor
      PrecalculatedAveragine() = default;

      /**
       @brief constructor with parameters such as mass ranges and bin size.
       @param min_mass the averagine distributions will be calculated from this min_mass
       @param max_mass to the max_mass
       @param delta with the bin size delta
       @param generator this generates (calculates) the distributions
       @param use_RNA_averagine if set, nucleotide-based isotope patters are calculated
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

      /// get isotope start index
      Size getLeftCountFromApex(const double mass) const;

      /// get isotope end index
      Size getRightCountFromApex(const double mass) const;

      /// get index of most abundant isotope
      Size getApexIndex(const double mass) const;


      /// get mass difference between avg and mono masses
      double getAverageMassDelta(const double mass) const;

    };

    /// struct for TopPIC identification (both PrSMs and proteoforms)
    struct OPENMS_DLLAPI TopPicItem
    {
    public:
      TopPicItem() = default;

      /// parse a single line of TopPIC output
      explicit TopPicItem(String in);

      /// the line string
      String str_;
      /// information from each column
      int prsm_id_;
      int spec_id_;
      int scan_;
      double rt_;
      int peak_count_;
      int charge_;
      double precursor_mass_;
      double adj_precursor_mass_;
      int proteform_id_ = -1;
      String protein_acc_ = "";
      int first_residue_;
      int last_residue_;
      std::vector<double> unexp_mod_;
      int matched_peaks_;
      int matched_frags_;
      double e_value_;
      double spec_q_value_;
      double proteofrom_q_value_;
      double intensity_;

      bool operator<(const TopPicItem &a) const;

      bool operator>(const TopPicItem &a) const;

      bool operator==(const TopPicItem &other) const;

    };


    /// log transformed peak. After deconvolution, all necessary information from deconvolution such as charge and isotope index is stored.
    struct OPENMS_DLLAPI LogMzPeak
    {
    public:
      /// original peak mz
      double mz = 0;
      /// original peak intensity
      double intensity = 0;
      /// log transformed mz
      double logMz = -1000;
      /// deteremined mass after deconvolution. NOT monoisotopic but only decharged
      double mass = .0;
      /// absolute charge (in case negative, is_positive is set to false
      int abs_charge = 0;
      /// is positive mode
      bool is_positive = true;
      /// isotope index
      int isotopeIndex = -1;

      /// default constructor
      LogMzPeak() = default;

      /**
               @brief constructor from Peak1D.
              @param positive determines the charge carrier mass
        */
      explicit LogMzPeak(const Peak1D &peak, const bool positive);

      /// copy constructor
      LogMzPeak(const LogMzPeak &) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT monoisotopic mass. Valid only when charge is set
      double getUnchargedMass();

      bool operator<(const LogMzPeak &a) const;

      bool operator>(const LogMzPeak &a) const;

      bool operator==(const LogMzPeak &other) const;

    };

    /**
        //       @brief Static function to calculate averagines. PrecalculatedAveragine class is constructed inside for the calculation.
        //       @param max_mass max mass
        //       @param use_RNA_averagine if set, nucleotide-based averagines are calculated */
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
