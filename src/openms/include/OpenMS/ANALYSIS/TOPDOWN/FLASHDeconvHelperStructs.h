// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

namespace OpenMS
{
  /**
   * @brief Wrapper struct for all the structs needed by the FLASHDeconv
   * Three structures are defined: PrecalculatedAveragine, TopPicItem, and LogMzPeak
   * i) PrecalculatedAveragine - to match observed isotopic envelope against theoretical one, theoretical envelope from
   * averagine model should be quickly calculated. To do so, precalculate averagines for different masses at the beginning of FLASHDeconv runs
   * ii) TopPicItem - represent TopPic identification. Currently used for QScore training.
   * iii) LogMzPeak - Log transformed peak from original peak. Contains information such as charge, isotope index, and uncharged mass.
   * @see FLASHDeconvAlgorithm
   * @reference: FeatureFinderAlgorithmPickedHelperStructs
   */

  struct OPENMS_DLLAPI FLASHDeconvHelperStructs
  {
    /// @brief Averagine patterns pre-calculated for speed up. Other variables are also calculated for fast cosine calculation
    class OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes_;
      /// L2 norms_ for masses
      std::vector<double> norms_;
      /// mass differences between average mass and monoisotopic mass
      std::vector<double> average_mono_mass_difference_;
      /// mass differences between most abundant mass and monoisotopic mass
      std::vector<double> abundant_mono_mass_difference_;
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
      /// calculate the mass bin index from mass
      Size massToIndex_(double const mass) const;
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

      /// get distribution for input mass. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      IsotopeDistribution get(const double mass) const;

      /// get max isotope index
      int getMaxIsotopeIndex() const;

      /// set max isotope index
      void setMaxIsotopeIndex(const int index);

      /// get isotope distance (from apex to the left direction) to consider. This is specified in the constructor. index. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getLeftCountFromApex(const double mass) const;

      /// get isotope distance (from apex to the rigth direction) to consider. This is specified in the constructor. index. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getRightCountFromApex(const double mass) const;

      /// get index of most abundant isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getApexIndex(const double mass) const;

      /// get index of last isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getLastIndex(const double mass) const;

      /// get mass difference between avg and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      double getAverageMassDelta(const double mass) const;

      /// get mass difference between most abundant mass and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      double getMostAbundantMassDelta(const double mass) const;
    };

    /// struct for TopPIC identification (both PrSMs and proteoforms)
    struct OPENMS_DLLAPI TopPicItem
    {
    public:
      TopPicItem() = default;

      /// parse a single line of TopPIC output. Header includes:
      /// Data file name	Prsm ID	Spectrum ID	Fragmentation	Scan(s)	Retention time	#peaks	Charge	Precursor mass
      /// Adjusted precursor mass	Proteoform ID	Feature intensity	Feature score	Protein accession	Protein description
      /// First residue	Last residue	Proteoform	#unexpected modifications	MIScore	#variable PTMs	#matched peaks
      /// #matched fragment ions	E-value	Spectrum-level Q-value	Proteoform-level Q-value
      explicit TopPicItem(String in);

      /// the line string
      String str;
      /// information from each column
      int prsm_id;
      int spec_id;
      int scan;
      double rt;
      int peak_count;
      int charge;
      double precursor_mass;
      double adj_precursor_mass;
      int proteform_id = -1;
      String protein_acc = "";
      int first_residue;
      int last_residue;
      std::vector<double> unexp_mod;
      int matched_peaks;
      int matched_frags;
      double e_value;
      double spec_q_value;
      double proteofrom_q_value;
      double intensity;

      /// scan numbers are compared
      bool operator<(const TopPicItem& a) const;

      bool operator>(const TopPicItem& a) const;

      bool operator==(const TopPicItem& other) const;
    };

    /// Mass feature (Deconvolved masses in spectra are traced by Mass tracing to generate mass features - like LC-MS features).
    struct OPENMS_DLLAPI MassFeature
    {
    public:
      MassTrace mt;
      std::vector<float> per_charge_intensity;
      std::vector<float> per_isotope_intensity;
      int iso_offset;
      int scan_number, rep_charge;
      double avg_mass;
      int min_charge, max_charge, charge_count;
      double isotope_score, qscore;
      double rep_mz;
    };

    /// log transformed peak. After deconvolution, all necessary information from deconvolution such as charge and isotope index is stored.
    class OPENMS_DLLAPI LogMzPeak
    {
    public:
      /// original peak mz
      double mz = 0;
      /// original peak intensity
      double intensity = 0;
      /// log transformed mz
      double logMz = -1000;
      /// determined mass after deconvolution. NOT monoisotopic but only decharged
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
      explicit LogMzPeak(const Peak1D& peak, const bool positive);

      /// copy constructor
      LogMzPeak(const LogMzPeak& ) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT a monoisotopic mass of a PeakGroup, rather a monoisotopic mass of each LogMzPeak. Returns 0 if no charge set
      double getUnchargedMass();

      /// log mz values are compared
      bool operator<(const LogMzPeak& a) const;

      bool operator>(const LogMzPeak& a) const;

      bool operator==(const LogMzPeak& other) const;

    };

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
