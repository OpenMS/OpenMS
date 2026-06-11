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
    /**
      @brief Pre-binned averagine isotope-distribution cache used by FLASHDeconv for fast cosine-similarity scoring.

      Builds a table keyed by a mass bin index — one @c IsotopeDistribution plus its derived
      quantities (apex index, left/right counts of significant isotopes, average / most-abundant
      mass offsets, an SNR multiplication factor) per mass bin from @c min_mass to @c max_mass
      at step @c delta. The constructor performs all the per-mass work up front; the accessors
      are O(1) table lookups via @c massToIndex_(mass) (rounded to the nearest bin).

      Bin lookup silently clamps: any mass below @c min_mass maps to bin @c 0 and any mass above
      @c max_mass maps to the last bin, so the accessors never throw on out-of-range input.
    */
    class OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// Per-bin isotope distribution. After trimming and L2 normalisation in the constructor.
      std::vector<IsotopeDistribution> isotopes_;
      /// Per-bin L2 norm — populated and consumed alongside @ref isotopes_ for fast cosine scoring (vestigial: not consumed by the public API in this header).
      std::vector<double> norms_;
      /// Per-bin difference between the average mass and the monoisotopic mass of the trimmed distribution.
      std::vector<double> average_mono_mass_difference_;
      /// Per-bin difference between the most-abundant-isotope mass and the monoisotopic mass.
      std::vector<double> abundant_mono_mass_difference_;
      /// Per-bin SNR multiplication factor: @c (sum-of-normalised-intensities)^2 of the trimmed distribution.
      std::vector<double> snr_mul_factor_;
      /// Per-bin count of significant isotopes on the left side of the apex (running maximum across bins; see the @c PrecalculatedAveragine parameterised constructor for the unusual running-max semantics).
      std::vector<int> left_count_from_apex_;
      /// Per-bin count of significant isotopes on the right side of the apex (running maximum across bins; same caveat as @ref left_count_from_apex_).
      std::vector<int> right_count_from_apex_;
      /// Per-bin index of the most-abundant isotope inside the trimmed distribution.
      std::vector<Size> apex_index_;

      /// Cap on the isotope index used externally by FLASHDeconv. NOT initialised by the constructor — call @ref setMaxIsotopeIndex before reading via @ref getMaxIsotopeIndex.
      Size max_isotope_index_;
      /// Bin step size (the @c delta passed to the constructor).
      double mass_interval_;
      /// Lower bound of the cached mass range (the @c min_mass passed to the constructor).
      double min_mass_;
      /// Convert @p mass to a bin index. Clamps to @c [0, isotopes_.size()-1] for out-of-range inputs (no throw).
      Size massToIndex_(double mass) const;

    public:
      /// Default-construct an empty cache (no bins populated; the accessors are not safe until the parameterised ctor is used).
      PrecalculatedAveragine() = default;

      /**
        @brief Precompute averagine isotope distributions across @p min_mass .. @p max_mass in steps of @p delta.

        For each bin mass @c m @c = @c i*delta in the range, calls
        @c CoarseIsotopePatternGenerator::estimateFromPeptideMonoWeight(m) (or
        @c estimateFromRNAMonoWeight(m) when @p use_RNA_averagine is @c true) and applies the
        following processing:
          - **Trimming**: iteratively drops the least-intense end peak (left or right) until
            the trimmed-out power is below @c 0.0001 of the total, with a minimum kept length
            of @c 2 peaks. Trimmed peaks have their intensity zeroed out, not removed from the
            distribution.
          - **L2 normalisation**: divides every remaining peak's intensity by
            @c sqrt(total_pwr) of the kept set.
          - **Decoy mode**: when @p decoy_iso_distance is @c > @c 0, every peak's m/z is
            multiplied by @p decoy_iso_distance before trimming, then the distribution is
            re-sorted by mass and renormalised. This produces a "stretched" pattern whose
            inter-isotope spacing is @p decoy_iso_distance times the natural spacing.
          - **Left/right counts**: clamped to a floor of @c 2 each, and then accumulated as a
            running maximum across bins — so all bins end up with the same width (the maximum
            seen so far during construction).

        Per bin, the constructor records the trimmed distribution, the apex index, the
        running-max left/right counts, the average-mono and most-abundant-mono mass deltas,
        and the SNR multiplication factor @c (sum-of-normalised-intensities)^2.

        @param[in]     min_mass            Lower bound of the cached mass range; bins are built starting at the first multiple of @p delta that is @c >= @c min_mass.
        @param[in]     max_mass            Upper bound of the cached mass range; bin generation stops at the first multiple of @p delta that exceeds @c max_mass.
        @param[in]     delta               Bin step size.
        @param[in,out] generator           Isotope-pattern generator used to compute each per-mass distribution.
        @param[in]     use_RNA_averagine   If @c true, use the RNA-mass averagine; otherwise use peptide averagine.
        @param[in]     decoy_iso_distance  @c > @c 0 enables decoy mode (per-peak m/z multiplied by this value); @c <= @c 0 disables it. Default @c -1.
      */
      PrecalculatedAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance = -1);

      /// Copy constructor.
      PrecalculatedAveragine(const PrecalculatedAveragine&) = default;

      /// Move constructor.
      PrecalculatedAveragine(PrecalculatedAveragine&& other) noexcept = default;

      /// Copy assignment.
      PrecalculatedAveragine& operator=(const PrecalculatedAveragine& pc) = default;

      /// Move assignment.
      PrecalculatedAveragine& operator=(PrecalculatedAveragine&& pc) noexcept = default;

      /// Destructor.
      ~PrecalculatedAveragine() = default;


      /**
        @brief Return the cached trimmed + L2-normalised isotope distribution for the bin containing @p mass.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Cached isotope distribution for the resolved bin.
      */
      IsotopeDistribution get(double mass) const;

      /**
        @brief Return the externally-set isotope-index cap.

        @warning The constructor does not initialise @c max_isotope_index_; reading this
                 before @ref setMaxIsotopeIndex has been called returns an uninitialised value.

        @return The most recent value passed to @ref setMaxIsotopeIndex.
      */
      size_t getMaxIsotopeIndex() const;

      /**
        @brief Set the isotope-index cap consulted by external callers. Does not influence the cached distributions or other accessors.

        @param[in] index New isotope-index cap value.
      */
      void setMaxIsotopeIndex(int index);

      /**
        @brief Return the number of significant isotopes to the left of the apex for the bin containing @p mass.

        @note The stored value is a running maximum accumulated across bins during construction,
              so it grows monotonically with @p mass (see the parameterised constructor's
              left/right-count semantics).

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Running-max left count of significant isotopes.
      */
      Size getLeftCountFromApex(double mass) const;

      /**
        @brief Return the number of significant isotopes to the right of the apex for the bin containing @p mass.

        @note Same running-maximum caveat as @ref getLeftCountFromApex.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Running-max right count of significant isotopes.
      */
      Size getRightCountFromApex(double mass) const;

      /**
        @brief Return the index of the most-abundant isotope inside the trimmed distribution for the bin containing @p mass.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Position of the apex peak inside the trimmed distribution.
      */
      Size getApexIndex(double mass) const;

      /**
        @brief Return @c apex_index @c + @c right_count_from_apex for the bin containing @p mass — the index of the rightmost significant isotope.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Index of the rightmost significant isotope.
      */
      Size getLastIndex(double mass) const;

      /**
        @brief Return @c (average_mass @c - @c monoisotopic_mass) for the bin's trimmed distribution.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Mass difference between the trimmed distribution's average mass and its monoisotopic mass.
      */
      double getAverageMassDelta(double mass) const;

      /**
        @brief Return @c (most_abundant_mass @c - @c monoisotopic_mass) for the bin's trimmed distribution.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Mass difference between the trimmed distribution's most-abundant-isotope mass and its monoisotopic mass.
      */
      double getMostAbundantMassDelta(double mass) const;

      /**
        @brief Return the SNR multiplication factor for the bin containing @p mass.

        @c (sum-of-normalised-intensities)^2 of the trimmed distribution. Used by FLASHDeconv
        for fast SNR computation.

        @param[in] mass Input mass to query. Out-of-range values are silently clamped to @c [min_mass, max_mass].
        @return Squared sum of normalised peak intensities for the bin.
      */
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