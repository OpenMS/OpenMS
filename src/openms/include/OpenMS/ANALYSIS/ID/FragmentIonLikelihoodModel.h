// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Per-run likelihood model of fragment ion presence and intensity rank, learned from confident PSMs.

    Every theoretical ion of a peptide-spectrum match is placed in a context made of its terminal
    series (prefix a/b/c or suffix x/y/z), the precursor charge (1-2, 3, at least 4), the fragment
    charge (1 or higher) and the relative position of the cleavage along the backbone (10 bins of
    fragment length / peptide length). Its outcome is either absence or the intensity rank of the
    matched peak within the spectrum (7 rank bins: 1-2, 3-5, 6-10, 11-20, 21-40, 41-80, above 80).

    Two sets of counts are kept: the ions of confident peptides ("signal") and the ions of their
    reversed sequences matched against the same spectra ("noise"). finalize() turns both into
    smoothed probabilities with pseudo-count back-off from the full context to the series and
    charge context, to the series and to a flat prior, so that a few thousand PSMs suffice and
    unseen contexts stay finite.

    A scored PSM yields three additive features: the summed log-likelihood ratio of its outcomes
    (present ions are credited, absent ions debited by what the run has taught about them), the
    fraction of the predicted ion presence that was observed, and the fraction of the ions most
    likely to be present that were observed. The model carries no pre-trained tables; it is
    meant to be trained on the run it scores, see ProSEAlgorithm (annotate:self_trained_ion_priors).
    Matching is closest-peak within a tolerance, as in HyperScore.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FragmentIonLikelihoodModel
  {
  public:
    static constexpr Size SERIES = 2;            ///< prefix (a/b/c) and suffix (x/y/z) ions
    static constexpr Size PRECURSOR_BUCKETS = 3; ///< precursor charge 1-2, 3, at least 4
    static constexpr Size FRAGMENT_CHARGES = 2;  ///< fragment charge 1, at least 2
    static constexpr Size POSITION_BINS = 10;    ///< fragment length / peptide length
    static constexpr Size RANK_BINS = 7;         ///< matched peak intensity rank bins
    static constexpr Size ABSENT = RANK_BINS;    ///< outcome of an unmatched ion
    static constexpr Size OUTCOMES = RANK_BINS + 1;
    static constexpr Size CONTEXTS = SERIES * PRECURSOR_BUCKETS * FRAGMENT_CHARGES * POSITION_BINS;

    /// Context of one theoretical ion
    struct Context
    {
      Size series = 0;          ///< 0 prefix, 1 suffix
      Size precursor_bucket = 0; ///< 0: charge 1-2, 1: charge 3, 2: charge >= 4
      Size fragment_charge = 0;  ///< 0: charge 1, 1: charge >= 2
      Size position_bin = 0;     ///< POSITION_BINS * fragment length / peptide length, capped
    };

    /// Features of one PSM
    struct Features
    {
      double log_likelihood_ratio = 0.0; ///< sum over all theoretical ions of ln P(outcome | signal) / P(outcome | noise)
      double explained_presence = 0.0;   ///< sum of P(present) over matched ions / sum of P(present) over all ions
      double top_predicted_observed = 0.0; ///< fraction of the top_k ions with the highest P(present) that matched
      Size matched_ions = 0;             ///< theoretical ions with a peak within tolerance
      Size theoretical_ions = 0;         ///< theoretical ions with a recognised ion name
    };

    /// Model with the default pseudo-count (20) of the back-off smoothing
    FragmentIonLikelihoodModel();

    /// Model with an explicit pseudo-count (must be positive)
    explicit FragmentIonLikelihoodModel(double pseudo_count);

    /// Intensity ranks (1 = most intense; ties by ascending m/z) of the peaks of @p spectrum, in peak order
    static std::vector<Size> intensityRanks(const MSSpectrum& spectrum);

    /// Outcome (rank bin) of a matched peak with the given 1-based intensity rank
    static Size rankOutcome(Size rank);

    /**
      @brief Context of a theoretical ion.

      @param prefix Whether the ion is an N-terminal (a/b/c) rather than a C-terminal (x/y/z) ion.
      @param precursor_charge Precursor charge of the PSM.
      @param fragment_charge Charge of the ion.
      @param fragment_length Number of residues of the fragment (the ion's ordinal).
      @param peptide_length Number of residues of the peptide.
    */
    static Context contextOf(bool prefix, int precursor_charge, int fragment_charge, Size fragment_length, Size peptide_length);

    /**
      @brief Series and ordinal of an ion name as written by TheoreticalSpectrumGenerator, e.g. "b5+", "y3++", "z.4+".

      Names carrying the ion type after a '$' (cross-link annotations) are supported. Returns false
      for names without a terminal series letter or ordinal.
    */
    static bool parseIonName(const std::string& name, bool& prefix, Size& ordinal);

    /**
      @brief Add the ions of one PSM to the signal or noise counts.

      @param spectrum Experimental spectrum, sorted by m/z.
      @param ranks intensityRanks(spectrum).
      @param theoretical Sorted theoretical spectrum with ion names (first StringDataArray) and charges
             (first IntegerDataArray), as TheoreticalSpectrumGenerator writes them with add_metainfo.
      @param peptide_length Number of residues of the peptide.
      @param precursor_charge Precursor charge of the PSM.
      @param tolerance Positive finite matching tolerance.
      @param ppm Whether the tolerance is in ppm rather than Da.
      @param noise Whether the ions belong to a reversed (noise) rather than a confident peptide.
      @throws Exception::InvalidValue if the theoretical spectrum lacks ion names.
      @throws Exception::InvalidParameter if the tolerance is not finite and positive.
    */
    void addObservations(const MSSpectrum& spectrum,
                         const std::vector<Size>& ranks,
                         const MSSpectrum& theoretical,
                         Size peptide_length,
                         int precursor_charge,
                         double tolerance,
                         bool ppm,
                         bool noise);

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
    double logLikelihoodRatio(const Context& context, Size outcome) const;

    /**
      @brief Probability that an ion of this context is present in the spectrum of its peptide.
      @throws Exception::Precondition if the model is not finalized.
    */
    double presenceProbability(const Context& context) const;

    /**
      @brief Features of one PSM. Parameters as in addObservations(); @p top_k bounds top_predicted_observed.
      @throws Exception::Precondition if the model is not finalized.
      @throws Exception::InvalidValue if the theoretical spectrum lacks ion names.
      @throws Exception::InvalidParameter if the tolerance is not finite and positive.
    */
    Features score(const MSSpectrum& spectrum,
                   const std::vector<Size>& ranks,
                   const MSSpectrum& theoretical,
                   Size peptide_length,
                   int precursor_charge,
                   double tolerance,
                   bool ppm,
                   Size top_k = 6) const;

  private:
    /// One theoretical ion after matching
    struct Ion_
    {
      Context context;
      Size outcome; ///< rank bin or ABSENT
    };

    /// Flat index of (context, outcome) into the count and probability tables
    static Size index_(const Context& context, Size outcome);

    /// Match every named theoretical ion against the spectrum
    static void matchIons_(const MSSpectrum& spectrum,
                           const std::vector<Size>& ranks,
                           const MSSpectrum& theoretical,
                           Size peptide_length,
                           int precursor_charge,
                           double tolerance,
                           bool ppm,
                           std::vector<Ion_>& ions);

    /// Smoothed log-probabilities of one count table, backing off over the context hierarchy
    static std::vector<double> smooth_(const std::vector<double>& counts, double pseudo_count);

    /// Throws unless finalize() was called after the last observation
    void requireTrained_() const;

    std::vector<double> signal_counts_; ///< CONTEXTS * OUTCOMES counts of confident ions
    std::vector<double> noise_counts_;  ///< CONTEXTS * OUTCOMES counts of reversed ions
    std::vector<double> signal_logp_;   ///< smoothed ln P(outcome | signal, context)
    std::vector<double> noise_logp_;    ///< smoothed ln P(outcome | noise, context)
    Size signal_psms_ = 0;
    Size noise_psms_ = 0;
    double pseudo_count_ = 20.0;
    bool finalized_ = false;
  };
} // namespace OpenMS
