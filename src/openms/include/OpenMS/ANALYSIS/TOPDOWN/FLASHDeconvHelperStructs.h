// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
   * ii) TopPicItem - represent TopPic identification. Currently used for setQscore training. TopPic is the top-down proteomics identification tool (https://www.toppic.org/).
   * iii) LogMzPeak - Log transformed peak from original peak. Contains information such as charge, isotope index, and uncharged mass.
   * @see SpectralDeconvolution
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

      std::vector<double> snr_mul_factor_;
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
      Size massToIndex_(double mass) const;

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
      PrecalculatedAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine);

      /// copy constructor
      PrecalculatedAveragine(const PrecalculatedAveragine&) = default;

      /// move constructor
      PrecalculatedAveragine(PrecalculatedAveragine&& other) noexcept = default;

      /// copy assignment operator
      PrecalculatedAveragine& operator=(const PrecalculatedAveragine& pc) = default;

      /// move assignment operator
      PrecalculatedAveragine& operator=(PrecalculatedAveragine&& pc) noexcept = default;

      /// destructor
      ~PrecalculatedAveragine() = default;


      /// get distribution for input mass. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      IsotopeDistribution get(double mass) const;

      /// get max isotope index
      size_t getMaxIsotopeIndex() const;

      /// set max isotope index
      void setMaxIsotopeIndex(int index);

      /// get isotope distance (from apex to the left direction) to consider. This is specified in the constructor. index. If input mass exceeds the maximum mass (specified in constructor), output for
      /// the maximum mass
      Size getLeftCountFromApex(double mass) const;

      /// get isotope distance (from apex to the rigth direction) to consider. This is specified in the constructor. index. If input mass exceeds the maximum mass (specified in constructor), output
      /// for the maximum mass
      Size getRightCountFromApex(double mass) const;

      /// get index of most abundant isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getApexIndex(double mass) const;

      /// get index of last isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getLastIndex(double mass) const;

      /// get mass difference between avg and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      double getAverageMassDelta(double mass) const;

      /// get mass difference between most abundant mass and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      double getMostAbundantMassDelta(double mass) const;

      double getSNRMultiplicationFactor(double mass) const;
    };

    /// Mass feature (Deconvolved masses in spectra are traced by Mass tracing to generate mass features - like LC-MS features).
    struct OPENMS_DLLAPI MassFeature
    {
    public:
      /// feature index;
      uint index;
      MassTrace mt;
      std::vector<float> per_charge_intensity;
      std::vector<float> per_isotope_intensity;
      int iso_offset;
      int scan_number, rep_charge;
      double avg_mass;
      int min_charge, max_charge, charge_count;
      double isotope_score, qscore;
      double rep_mz;
      bool is_decoy;
      uint ms_level;
      /// features are compared
      bool operator<(const MassFeature& a) const
      {
        return avg_mass < a.avg_mass;
      }
      bool operator>(const MassFeature& a) const
      {
        return avg_mass > a.avg_mass;
      }
      bool operator==(const MassFeature& other) const
      {
        return avg_mass == other.avg_mass;
      }
    };

    /// Isobaric quantities.
    struct OPENMS_DLLAPI IsobaricQuantities {
    public:
      int scan;
      double rt;
      double precursor_mz;
      double precursor_mass;
      std::vector<double> quantities;
      std::vector<double> merged_quantities;

      bool empty() const;
    };

    /// log transformed peak. After deconvolution, all necessary information from deconvolution such as charge and isotope index is stored.
    class OPENMS_DLLAPI LogMzPeak
    {
    public:
      /// original peak mz
      double mz = 0;
      /// original peak intensity
      float intensity = 0;
      /// log transformed mz
      double logMz = -1000;
      /// determined mass after deconvolution. NOT monoisotopic but only decharged
      double mass = .0;
      /// absolute charge (in case negative, is_positive is set to false)
      int abs_charge = 0;
      /// is positive mode
      bool is_positive = true;
      /// isotope index
      int isotopeIndex = -1;

      /// default constructor
      LogMzPeak() = default;

      /**
        @brief constructor from Peak1D.
        @param peak the original spectral peak
        @param positive determines the charge carrier mass. Can be obtained by getChargeMass(true) for positive mode (Constants::PROTON_MASS_U) and getChargeMass(false) for negative mode
        (-Constants::PROTON_MASS_U)
      */
      explicit LogMzPeak(const Peak1D& peak, bool positive);

      /// copy constructor
      LogMzPeak(const LogMzPeak&) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT a monoisotopic mass of a PeakGroup, rather a monoisotopic mass of each LogMzPeak. Returns 0 if no charge set
      double getUnchargedMass() const;

      /// log mz values are compared
      bool operator<(const LogMzPeak& a) const;
      bool operator>(const LogMzPeak& a) const;
      bool operator==(const LogMzPeak& other) const;
    };

    /// Sequence tag. No mass gap is allowed in the seq. The mass gap containing tag should be enumerated into multiple Tag instances from outside.
    class OPENMS_DLLAPI Tag
    {
    public:
      /// constructor
      explicit Tag(String  seq, double n_mass, double c_mass, std::vector<int>& scores, std::vector<double>& mzs);

      /// copy constructor
      Tag(const Tag&) = default;

      /// destructor
      ~Tag() = default;

      bool operator<(const Tag& a) const;
      bool operator>(const Tag& a) const;
      bool operator==(const Tag& other) const;

      String getSequence() const;
      double getNtermMass() const;
      double getCtermMass() const;
      Size getLength() const;
      int getScore() const;
      int getScore(int pos) const;
      String toString() const;
      const std::vector<double>& getMzs() const;

    private:
      String seq_;
      double n_mass_ = -1, c_mass_ = -1;
      std::vector<int> scores_;
      std::vector<double> mzs_;
      Size length_;
    };

    /**
       @brief calculate log mzs from mzs
       @param mz mz
       @param positive determines the charge carrier mass
     */
    static double getLogMz(double mz, bool positive);

    /**
       @brief get charge carrier mass : positive mode mass of (Constants\::PROTON_MASS_U) and negative mode mass of (-Constants\::PROTON_MASS_U)
       @param positive_ioniziation_mode Determines the charge carrier mass (true = positive or false = negative)
    */
    static float getChargeMass(bool positive_ioniziation_mode);
  };
} // namespace OpenMS
