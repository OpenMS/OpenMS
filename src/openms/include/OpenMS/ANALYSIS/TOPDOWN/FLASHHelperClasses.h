// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#ifdef _MSC_VER
using uint = unsigned int; // POSIX uint not available on MSVC; was provided transitively via Qt
#endif
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <boost/dynamic_bitset.hpp>
#include <functional>

namespace OpenMS
{
  /**
   * @brief Wrapper struct for all the structs needed by the FLASHDeconv
   * The following structures/classes are defined: PrecalculatedAveragine, MassFeature, IsobaricQuantities, LogMzPeak,
   * Tag, and DAG for tagging.
   * i) PrecalculatedAveragine - to match observed isotopic envelope against theoretical one, theoretical envelope from
   * averagine model should be quickly calculated. To do so, precalculate averagines for different masses at the beginning of FLASHDeconv runs
   * ii) MassFeature - the feature of the deconvolved masses (instead of m/zs for MS features).
   * iii) IsobaricQuantities - the quantities from isobaric quantification (implemented within FLASHDeconv)
   * iv) LogMzPeak - Log transformed peak from original peak. Contains information such as charge, isotope index, and
   * uncharged mass.
   * v) Tag - the sequence tag generated on deconvolved spectra in FLASHTaggerAlgorithm
   * vi) DAG - the graph for sequence tagging in FLASHTaggerAlgorithm and for proteoform characterization in FLASHExtenderAlgorithm
   * @see SpectralDeconvolution
   */

  struct OPENMS_DLLAPI FLASHHelperClasses
  {
    /// @brief Averagine patterns pre-calculated for speed up. Other variables are also calculated for fast cosine calculation
    class OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes_;
      /// L2 norms_ for masses - for quick isotope cosine calculation
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
       @param[in] min_mass the averagine distributions will be calculated from this min_mass
       @param[in] max_mass to the max_mass
       @param[in] delta with the bin size delta
       @param[in] generator this generates (calculates) the distributions
       @param[in] use_RNA_averagine if set, nucleotide-based isotope patters are calculated
       @param[in] decoy_iso_distance if set to a positive value, nonsensical isotope patterns are generated - the distance between isotope = decoy_iso_distance * normal distance.
    */
      PrecalculatedAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance = -1);

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

      /// get isotope distance (from apex to the left direction) to consider. This is specified in the constructor. index. If input mass exceeds the
      /// maximum mass (specified in constructor), output for the maximum mass
      Size getLeftCountFromApex(double mass) const;

      /// get isotope distance (from apex to the right direction) to consider. This is specified in the constructor. index. If input mass exceeds the
      /// maximum mass (specified in constructor), output for the maximum mass
      Size getRightCountFromApex(double mass) const;

      /// get index of most abundant isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getApexIndex(double mass) const;

      /// get index of last isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getLastIndex(double mass) const;

      /// get mass difference between avg and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum
      /// mass
      double getAverageMassDelta(double mass) const;

      /// get mass difference between most abundant mass and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for
      /// the maximum mass
      double getMostAbundantMassDelta(double mass) const;

      /// get the SNR multiplication factor - used for quick SNR calculation
      double getSNRMultiplicationFactor(double mass) const;
    };

    /// Mass feature (Deconvolved masses in spectra are traced by Mass tracing to generate mass features - like LC-MS features).
    struct OPENMS_DLLAPI MassFeature
    {
    public:
      /// feature index;
      uint index;
      /// the trace calculated from the masses
      MassTrace mt;
      /// per charge and isotope intensities
      std::vector<float> per_charge_intensity;
      std::vector<float> per_isotope_intensity;
      /// isotope offset between deconvolved masses and mass feature
      int iso_offset;
      /// representative scan number, min and max scan numbers, and representative charge
      int scan_number, min_scan_number, max_scan_number, rep_charge;
      /// average mass
      double avg_mass;
      /// minimum maximum charges, and number of charges
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
    struct OPENMS_DLLAPI IsobaricQuantities
    {
    public:
      int scan;
      double rt;
      double precursor_mz;
      double precursor_mass;
      std::vector<double> quantities;
      std::vector<double> merged_quantities;
      /// return true if no isobaric quantities have been stored
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
        @param[in] peak the original spectral peak
        @param[in] positive determines the charge carrier mass. Can be obtained by getChargeMass(true) for positive mode (Constants::PROTON_MASS_U) and
        getChargeMass(false) for negative mode
        (-Constants::PROTON_MASS_U)
      */
      explicit LogMzPeak(const Peak1D& peak, bool positive);

      /// copy constructor
      LogMzPeak(const LogMzPeak&) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT a monoisotopic mass of a PeakGroup, rather the mass of each LogMzPeak. Returns 0 if no
      /// charge is set
      double getUnchargedMass() const;

      /// log mz values are compared
      bool operator<(const LogMzPeak& a) const;
      bool operator>(const LogMzPeak& a) const;
      bool operator==(const LogMzPeak& other) const;
    };

    /**
       @brief calculate log mzs from mzs
       @param[in] mz mz
       @param[in] positive determines the charge carrier mass
     */
    static double getLogMz(double mz, bool positive);

    /**
       @brief get charge carrier mass : positive mode mass of (Constants\::PROTON_MASS_U) and negative mode mass of (-Constants\::PROTON_MASS_U)
       @param[in] positive_ioniziation_mode Determines the charge carrier mass (true = positive or false = negative)
    */
    static float getChargeMass(bool positive_ioniziation_mode);
  };
} // namespace OpenMS

namespace std
{
  /// @brief Hash specialization for FLASHHelperClasses::MassFeature
  template<>
  struct hash<OpenMS::FLASHHelperClasses::MassFeature>
  {
    std::size_t operator()(const OpenMS::FLASHHelperClasses::MassFeature& mf) const noexcept
    {
      // Hash based on avg_mass (the field used in operator==)
      return OpenMS::hash_float(mf.avg_mass);
    }
  };

  /// @brief Hash specialization for FLASHHelperClasses::LogMzPeak
  template<>
  struct hash<OpenMS::FLASHHelperClasses::LogMzPeak>
  {
    std::size_t operator()(const OpenMS::FLASHHelperClasses::LogMzPeak& peak) const noexcept
    {
      // Hash based on logMz and intensity (the fields used in operator==)
      std::size_t seed = OpenMS::hash_float(peak.logMz);
      OpenMS::hash_combine(seed, OpenMS::hash_float(peak.intensity));
      return seed;
    }
  };
} // namespace std