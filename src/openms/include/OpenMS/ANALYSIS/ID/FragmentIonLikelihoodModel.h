// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS
{
  class AASequence;

  /**
    @brief Per-run likelihood model of fragment ion presence, intensity rank and mass error, learned from confident PSMs.

    Every theoretical ion of a peptide-spectrum match is placed in a context and gets an outcome. Two context sets
    are available (see ContextSet):
    - BASIC: terminal series (prefix a/b/c or suffix x/y/z), precursor charge (1-2, 3, at least 4), fragment charge
      (1 or higher) and relative position of the cleavage along the backbone (10 bins of fragment length / peptide
      length). Outcome: absent, or the intensity rank of the matched peak within the spectrum (7 rank bins: 1-2, 3-5,
      6-10, 11-20, 21-40, 41-80, above 80).
    - RICH: BASIC plus the residues at the cleavage site (bond X|P, bond D/E|X, any other bond) and whether the
      complementary singly charged ion was matched. Outcome: absent, or rank bin x mass-error bin (|error| / tolerance
      below 0.25, below 0.5, up to 1). This is a self-trained surrogate of the model-free part of the rich-ion
      likelihood of ANDES (ion context, complementary-ion coherence, matched rank and mass error).

    Two sets of counts are kept: the ions of confident peptides ("signal") and the ions of their reversed sequences
    matched against the same spectra ("noise"). finalize() turns both into smoothed probabilities with pseudo-count
    back-off over a context hierarchy (BASIC: full context -> series, precursor and fragment charge -> series -> flat
    prior; RICH: full context -> series, precursor and fragment charge, cleavage site, complement -> series and
    fragment charge -> flat prior), so that a few thousand PSMs suffice and unseen contexts stay finite.

    A scored PSM yields three additive features: the summed log-likelihood ratio of its outcomes (present ions are
    credited, absent ions debited by what the run has taught about them), the fraction of the predicted ion presence
    that was observed, and the fraction of the ions most likely to be present that were observed. The model carries no
    pre-trained tables; it is meant to be trained on the run it scores, with the PSMs it scores held out from its
    training data (see ProSEAlgorithm, annotate:self_trained_ion_priors, which trains one model per half of the
    spectra and scores each half with the model of the other half).

    The peaks are read from PeakLists, a compact copy of the spectra (m/z and intensity-rank bin per peak). Matching
    is closest-peak within the tolerance, as in HyperScore.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FragmentIonLikelihoodModel
  {
  public:
    /// Ion contexts and outcomes the model distinguishes
    enum class ContextSet
    {
      BASIC, ///< series, precursor charge, fragment charge, position; outcome: intensity rank
      RICH   ///< BASIC plus cleavage-site residues and complementary-ion presence; outcome: intensity rank x mass error
    };

    static constexpr Size SERIES = 2;            ///< prefix (a/b/c) and suffix (x/y/z) ions
    static constexpr Size PRECURSOR_BUCKETS = 3; ///< precursor charge 1-2, 3, at least 4
    static constexpr Size FRAGMENT_CHARGES = 2;  ///< fragment charge 1, at least 2
    static constexpr Size POSITION_BINS = 10;    ///< fragment length / peptide length
    static constexpr Size CLEAVAGE_SITES = 3;    ///< RICH: bond before P, bond after D or E, any other bond
    static constexpr Size COMPLEMENT_STATES = 2; ///< RICH: complementary singly charged ion not matched / matched
    static constexpr Size RANK_BINS = 7;         ///< matched peak intensity rank bins
    static constexpr Size ERROR_BINS = 3;        ///< RICH: |mass error| / tolerance below 0.25, below 0.5, up to 1

    /**
      @brief Compact peak lists of the spectra of one run: per peak only its m/z and its intensity-rank bin.

      Five bytes per peak (m/z as float, rank bin as byte) instead of a copy of the spectra, so the peak lists of a
      whole run can be kept from preprocessing until the PSMs are annotated. The float m/z deviates from the double by
      at most 0.06 ppm, far below any fragment tolerance.
    */
    class OPENMS_DLLAPI PeakLists
    {
    public:
      /// Remove all peak lists and make room for @p spectra empty ones
      void reset(Size spectra);

      /// Remove all peak lists and release their memory
      void clear();

      /// Number of peak lists (spectra)
      Size size() const { return mz_.size(); }

      /**
        @brief Store the peaks of @p spectrum as list @p index, replacing what was there.

        The intensity ranks are taken among the peaks of @p spectrum (1 = most intense; ties by ascending m/z).
        Different indices may be assigned concurrently.

        @throws Exception::IndexOverflow if @p index is not below size()
        @throws Exception::IllegalArgument if @p spectrum is not sorted by m/z
      */
      void assign(Size index, const MSSpectrum& spectrum);

      /// Number of peaks of list @p index
      Size peaks(Size index) const { return mz_[index].size(); }

      /// m/z of peak @p peak of list @p index
      double mz(Size index, Size peak) const { return mz_[index][peak]; }

      /// Intensity-rank bin (see rankBin()) of peak @p peak of list @p index
      Size rankBin(Size index, Size peak) const { return rank_bins_[index][peak]; }

      /**
        @brief The peak of list @p index nearest to @p mz, if it lies within [mz - tolerance, mz + tolerance].

        Same rule as MSSpectrum::findNearest(mz, tolerance): of two equally near peaks the lower one.

        @return the peak index, or -1
      */
      Int findNearest(Size index, double mz, double tolerance) const;

      /// Number of peaks over all lists
      Size totalPeaks() const;

      /// Bytes held by the peak lists (peak data, list headers)
      Size memoryUsage() const;

    private:
      std::vector<std::vector<float>> mz_;        ///< m/z per peak, ascending, per spectrum
      std::vector<std::vector<std::uint8_t>> rank_bins_; ///< intensity-rank bin per peak, per spectrum
    };

    /// One theoretical ion after matching
    struct Ion
    {
      UInt32 context; ///< flat context index, below contexts()
      UInt32 outcome; ///< outcome index, below outcomes(); absentOutcome() if unmatched
    };

    /// Features of one PSM
    struct Features
    {
      double log_likelihood_ratio = 0.0;   ///< sum over all theoretical ions of ln P(outcome | signal) / P(outcome | noise)
      double explained_presence = 0.0;     ///< sum of P(present) over matched ions / sum of P(present) over all ions
      double top_predicted_observed = 0.0; ///< fraction of the top_k ions with the highest P(present) that matched
      Size matched_ions = 0;               ///< theoretical ions with a peak within tolerance
      Size theoretical_ions = 0;           ///< theoretical ions with a recognised ion name
    };

    /**
      @brief Untrained model.

      @param context_set Contexts and outcomes to distinguish.
      @param pseudo_count Pseudo-count of the back-off smoothing (must be finite and positive).
      @throws Exception::InvalidParameter if the pseudo-count is not finite and positive
    */
    explicit FragmentIonLikelihoodModel(ContextSet context_set = ContextSet::RICH, double pseudo_count = 20.0);

    /// The context set of this model
    ContextSet contextSet() const { return context_set_; }

    /// Number of ion contexts
    Size contexts() const { return contexts_; }

    /// Number of outcomes, absence included
    Size outcomes() const { return outcomes_; }

    /// Outcome of an unmatched ion
    Size absentOutcome() const { return outcomes_ - 1; }

    /// Rank bin (0 .. RANK_BINS - 1) of a matched peak with the given 1-based intensity rank
    static Size rankBin(Size rank);

    /// Error bin (0 .. ERROR_BINS - 1) of a match with the given |mass error| / tolerance (at most 1)
    static Size errorBin(double relative_error);

    /**
      @brief Series and ordinal of an ion name as written by TheoreticalSpectrumGenerator, e.g. "b5+", "y3++", "z.4+".

      Names carrying the ion type after a '$' (cross-link annotations) are supported. Returns false
      for names without a terminal series letter or ordinal.
    */
    static bool parseIonName(const std::string& name, bool& prefix, Size& ordinal);

    /**
      @brief Match the named ions of a theoretical spectrum against peak list @p index: one Ion per recognised ion.

      @param peaks Peak lists of the run.
      @param index Peak list of the PSM's spectrum.
      @param theoretical Sorted theoretical spectrum with ion names (first StringDataArray) and charges (first
             IntegerDataArray), as TheoreticalSpectrumGenerator writes them with add_metainfo. Without charges every
             ion counts as singly charged.
      @param peptide The peptide of the theoretical spectrum (length; residues at the cleavage sites for RICH).
      @param precursor_charge Precursor charge of the PSM.
      @param tolerance Positive finite matching tolerance.
      @param ppm Whether the tolerance is in ppm rather than Da.
      @param[out] ions The matched ions, in the order of the theoretical spectrum.
      @throws Exception::InvalidValue if the theoretical spectrum lacks ion names
      @throws Exception::InvalidParameter if the tolerance is not finite and positive
    */
    void matchIons(const PeakLists& peaks,
                   Size index,
                   const MSSpectrum& theoretical,
                   const AASequence& peptide,
                   int precursor_charge,
                   double tolerance,
                   bool ppm,
                   std::vector<Ion>& ions) const;

    /// Add the matched ions of one PSM to the signal (@p noise false) or noise counts
    void addObservations(const std::vector<Ion>& ions, bool noise);

    /// Turn the counts into smoothed log-probabilities. Required before scoring; observations added later need another call.
    void finalize();

    /// Whether finalize() was called after the last observation
    bool isTrained() const { return finalized_; }

    /// Number of confident PSMs added
    Size signalPsms() const { return signal_psms_; }

    /// Number of reversed PSMs added
    Size noisePsms() const { return noise_psms_; }

    /// Pseudo-count of the back-off smoothing
    double pseudoCount() const { return pseudo_count_; }

    /**
      @brief ln P(outcome | signal, context) - ln P(outcome | noise, context)
      @throws Exception::Precondition if the model is not finalized.
    */
    double logLikelihoodRatio(Size context, Size outcome) const;

    /**
      @brief Probability that an ion of this context is present in the spectrum of its peptide.
      @throws Exception::Precondition if the model is not finalized.
    */
    double presenceProbability(Size context) const;

    /**
      @brief Features of one PSM from its matched ions (see matchIons()); @p top_k bounds top_predicted_observed.
      @throws Exception::Precondition if the model is not finalized.
    */
    Features score(const std::vector<Ion>& ions, Size top_k = 6) const;

  private:
    /// Contexts of the coarser back-off levels: level 1 is the coarsest above the flat prior
    Size level1Context_(Size context) const;
    Size level2Context_(Size context) const;

    /// Smoothed log-probabilities of one count table, backing off over the context hierarchy
    std::vector<double> smooth_(const std::vector<double>& counts) const;

    /// Throws unless finalize() was called after the last observation
    void requireTrained_() const;

    ContextSet context_set_;
    Size sites_;      ///< cleavage-site classes of the context set (1 for BASIC)
    Size complement_; ///< complement states of the context set (1 for BASIC)
    Size errors_;     ///< error bins of the context set (1 for BASIC)
    Size contexts_;
    Size outcomes_;
    std::vector<double> signal_counts_; ///< contexts_ * outcomes_ counts of confident ions
    std::vector<double> noise_counts_;  ///< contexts_ * outcomes_ counts of reversed ions
    std::vector<double> signal_logp_;   ///< smoothed ln P(outcome | signal, context)
    std::vector<double> noise_logp_;    ///< smoothed ln P(outcome | noise, context)
    Size signal_psms_ = 0;
    Size noise_psms_ = 0;
    double pseudo_count_ = 20.0;
    bool finalized_ = false;
  };
} // namespace OpenMS
