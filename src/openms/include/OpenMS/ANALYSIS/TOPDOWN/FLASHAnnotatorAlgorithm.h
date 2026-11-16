// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <array>
#include <unordered_map>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace OpenMS
{
/**
@brief Annotates proteoform-spectrum matches (PrSMs) on the raw (centroided, not deconvolved) MS2 spectrum.

The theoretical isotope envelope of every terminal fragment of the proteoform is placed at every plausible charge on
the centroided spectrum, the way ExD Viewer annotates, and kept if it matches:
 - isotope peaks: averagine (FLASHHelperClasses::PrecalculatedAveragine, 25 Da bins) peaks of at least
   min_rel_isotope_intensity of the most abundant one, at mono + k * C13C12_MASSDIFF_U;
 - charges 1 .. min(precursor charge, floor(fragment mass / min_mass_per_charge));
 - each predicted peak is matched to the nearest observed peak within tol ppm (after calibration);
 - an envelope is accepted when at least min_peaks peaks, the most abundant one included, are matched and the cosine and
   the Pearson correlation between predicted and observed intensities (unmatched = 0) pass min_cos and min_pearson.

Blind modifications: a fragment carries a modification mass only if it contains the whole modification range;
fragments that split a range are not annotated.

Calibration (optional): a first pass at calibration_tol over the b/y fragments takes, for every envelope with at least
min_peaks matched peaks, the ppm error of the monoisotopic peak (if kept and matched) or else of the most abundant one;
the median error corrects the observed m/z. Fewer than 5 such envelopes: no correction.

Score: a raw match count grows with the number of fragments tried and with the density of the spectrum. The chance
match rate p0 is therefore measured per PrSM, in its own spectrum, with the same rule and charges on decoy fragment sets
(every target envelope, i.e. its isotope pattern and charges, moved by each of decoy_shifts):
p0 = (decoy matched hypotheses + 1) / (K * n + 2), with K the number of shifts and n the fragment hypotheses. A hypothesis
is an ion type and a length with at least min_peaks predicted isotope peaks and at least one charge whose most abundant
peak falls in the spectrum m/z range; ion types whose masses differ by one to three isotope spacings (z, z+1 and z+2) put
their envelopes on the same peaks and form one hypothesis. A hypothesis counts once, whatever the number of its matched
charges and ion types, and decoys are counted only for the n target hypotheses. With k matched hypotheses,
 - score = -log10 P(X >= k), X ~ Binomial(n, p0): the length enters only through n, so the score does not grow with the
   proteoform length. It is a p-value only as far as the decoy fragments match as often as random fragments do. The
   default shifts are whole multiples of 1.000495 Da, the spacing of the peptide mass bands, so that the decoy masses lie
   in the bands like the masses of random fragments; shifts between the bands underestimate p0;
 - llr = k ln(p1 / p0) + (n - k) ln((1 - p1) / (1 - p0)).
Both are evidence scores of one PrSM in its spectrum, not identification p-values: a candidate selected by a search
(tags, blind modifications) matches more than a random proteoform even when it is wrong. For identification, use
target-decoy competition, e.g. the scores of decoy-protein PrSMs.

@ingroup Topdown
*/
class OPENMS_DLLAPI FLASHAnnotatorAlgorithm : public DefaultParamHandler
{
public:
  /// a centroided spectrum as flat arrays, sorted by m/z, without zero-intensity peaks
  struct OPENMS_DLLAPI PeakArrays
  {
    std::vector<double> mzs;
    std::vector<double> intensities;
    /// per m/z bucket b (bucket_width Th wide), the index of the first peak at or above b * bucket_width (indexBuckets; empty: not built).
    /// A peak search then looks only inside the bucket of its lower end instead of the whole spectrum.
    std::vector<int> bucket_start;
    static constexpr double bucket_width = 1.0;
    /// builds bucket_start from mzs (ascending)
    void indexBuckets();
  };

  /// a blind modification: its mass and its 0-based inclusive residue range in the proteoform
  struct OPENMS_DLLAPI Modification
  {
    double mass = 0;
    int start = 0;
    int end = 0;
  };

  /// a PrSM to annotate
  struct OPENMS_DLLAPI PrSM
  {
    /// proteoform sequence (one-letter codes, without modifications)
    std::string sequence;
    std::vector<Modification> modifications;
    /// fragment charges are capped by this charge (<= 0: only by min_mass_per_charge)
    int precursor_charge = 0;
    /// isolation target m/z of the precursor (0: unknown). 1000 / precursor_mz is its charge density in charges per kDa, whatever
    /// its oligomer state; above internal_charge_switch the internal-ion thresholds after Asp and before Pro are exchanged
    double precursor_mz = 0;
    /// if false (undetermined start or end position), prefix or suffix fragments are not annotated
    bool known_n_term = true;
    bool known_c_term = true;
    /// index of the spectrum in the spectra given to annotate(); negative if the spectrum is missing
    int spectrum_index = -1;
  };

  /// a terminal fragment hypothesis
  struct OPENMS_DLLAPI Fragment
  {
    /// index into getIonTypes()
    int ion_type = 0;
    /// number of residues
    int length = 0;
    /// the bond (1 .. L-1) the fragment supports
    int bond = 0;
    /// neutral monoisotopic mass
    double mono_mass = 0;
    /// highest charge tried
    int max_charge = 0;
  };

  /// an accepted isotope envelope
  struct OPENMS_DLLAPI Envelope
  {
    /// index into getIonTypes(), or for an internal fragment into getInternalIonTypes()
    int ion_type = 0;
    /// first residue (0-based) of an internal fragment; -1 for a terminal fragment
    int first = -1;
    int length = 0;
    int bond = 0;
    int charge = 0;
    double mono_mass = 0;
    /// theoretical m/z of the most abundant isotope peak
    double most_abundant_mz = 0;
    int matched_peaks = 0;
    int predicted_peaks = 0;
    double cosine = 0;
    double pearson = 0;
    /// summed intensity of the matched peaks
    double intensity = 0;
    /// ppm error of the most abundant peak after calibration
    double ppm_error = 0;
    /// fraction of the spectrum's peaks that are not more intense than the matched most abundant peak
    double intensity_rank = 0;
    /// observed peaks per Th within 1 Th of the most abundant peak (the local density sets the chance of a match)
    double local_density = 0;
    /// 0 for a target envelope; the mass shift of a decoy envelope (written only with report_decoy_envelopes)
    double decoy_shift = 0;
  };

  /**
    annotation of a PrSM. Values are -1 when not computed (no spectrum, unknown residue; for the scores also no decoy shifts
    or no fragment hypothesis), except llr and calibration_ppm, which can be negative and are then NaN.
  */
  struct OPENMS_DLLAPI Result
  {
    bool annotated = false;
    std::vector<Envelope> envelopes;
    int covered_bonds = -1;
    /// covered bonds / (L-1) * 100
    double coverage = -1;
    /// k: fragment hypotheses (ion type group, length) with an accepted envelope at any charge
    int matched_fragments = -1;
    /// n: fragment hypotheses (ion type group, length) that can be matched (enough predicted peaks, a charge in the spectrum m/z range)
    int fragment_hypotheses = -1;
    /// matched decoy hypotheses, summed over the decoy shifts
    int decoy_matched_fragments = -1;
    /// p0: the chance rate of the spectrum, from the decoy fragments of all PrSMs annotated on it in the same annotate() call
    double chance_rate = -1;
    /// -log10 P(X >= k), X ~ Binomial(n, p0)
    double score = -1;
    /// -log10 of the combined p-value of the calibrated terminal and internal evidence (Fisher); -1 if not computed
    double fa_score = -1;
    double llr = std::numeric_limits<double>::quiet_NaN();
    double calibration_ppm = std::numeric_limits<double>::quiet_NaN();
    /// number of envelopes the calibration was computed from
    int calibration_envelopes = 0;
    /// internal fragments (with internal_ions): pieces with both ends inside the proteoform that are tested (n), counted (k), counted
    /// in the decoy sets (summed; only with report_decoy_envelopes, else -1), the calibrated chance of at least one counted piece, the
    /// calibrated -log10 P(K >= k | n), and the authentic internal cleavages (with a ladder)
    int internal_hypotheses = -1;
    int internal_matched = -1;
    int internal_decoy_matched = -1;
    double internal_chance_rate = -1;
    double internal_score = -1;
    int authentic_cleavages = -1;
    /// bonds covered by terminal and accepted internal fragments (both ends of an internal fragment), and their percentage
    int covered_bonds_with_internal = -1;
    double coverage_with_internal = -1;
    /// bonds that the internal pieces of a decoy shift would add to the terminal coverage, mean over the shifts that measure chance
    /// (with internal_alpha the half that did not choose the tag lengths; else all, only with report_decoy_envelopes); -1 if not computed
    double internal_chance_bonds = -1;
    /// pieces a tag needs after Asp, before Pro and elsewhere for this PrSM (0: no cleavage of the class is accepted)
    int internal_tag_pieces_asp = -1;
    int internal_tag_pieces_pro = -1;
    int internal_tag_pieces_other = -1;
    /// tolerance (ppm) of the internal pieces around the line offset + slope * m/z (ppm) fitted to the terminal envelopes
    double internal_tol_ppm = -1;
    double internal_line_offset = 0;
    double internal_line_slope = 0;
  };

  /// constructor
  FLASHAnnotatorAlgorithm();

  /// destructor
  ~FLASHAnnotatorAlgorithm() override = default;

  /// sorted flat arrays of the peaks of @p spec with positive intensity
  static void toPeakArrays(const MSSpectrum& spec, PeakArrays& peaks);

  /**
    @brief Annotate @p prsms on @p spectra (centroided) in parallel over the PrSMs.

    The averagine table is computed at the first call, up to the heaviest fragment (at most max_mass), and extended only
    if a later call needs heavier fragments.
    @param spectra spectra referred to by PrSM::spectrum_index
    @param prsms the PrSMs
    @param results one result per PrSM, in the same order
  */
  void annotate(const std::vector<PeakArrays>& spectra, const std::vector<PrSM>& prsms, std::vector<Result>& results);

  /**
    @brief The terminal fragments of @p prsm that are annotated; false if the sequence contains an unknown residue.
    @param span_fragments if given, the terminal fragments that contain only part of a modification range (their mass is unknown): with the
    modifications they contain whole, and without the ones they split. They are not annotated, but count as fragment hypotheses.
  */
  bool getFragments(const PrSM& prsm, std::vector<Fragment>& fragments, std::vector<Fragment>* span_fragments = nullptr) const;

  /// theoretical m/z and relative intensities (most abundant = 1) of the kept isotope peaks of a fragment
  void getIsotopeEnvelope(double mono_mass, int charge, std::vector<double>& mzs, std::vector<double>& intensities);

  /// the ion types, in the order Fragment::ion_type and Envelope::ion_type refer to
  const std::vector<std::string>& getIonTypes() const;

  /// the internal ion types (C-end type/N-end type, e.g. "b/y"), in the order Envelope::ion_type of an internal fragment refers to
  const std::vector<std::string>& getInternalIonTypes() const;

  /// neutral monoisotopic mass of the proteoform (residues + water + modifications); -1 for an unknown residue
  static double getProteoformMass(const std::string& sequence, const std::vector<Modification>& modifications);

  /// -log10 P(X >= k) for X ~ Binomial(n, p), computed in log space and capped at max_score_
  static double binomialTailScore(int n, int k, double p);

  /**
    @brief the MSn spectra of the given scan numbers of a raw mzML file, read in one pass

    Only those spectra are kept in memory. Scan numbers come from the native IDs; when they carry none, they are the spectrum
    index + 1 (the file is then read a second time, all spectra). The spectra are as in the file (profile or centroided).
    @return the number of spectra read
  */
  static Size loadSpectra(const std::string& mzml_file, const std::set<int>& scans, std::map<int, MSSpectrum>& spectra);

  /// the peak arrays of a raw spectrum, centroided with PeakPickerHiRes if it is a profile spectrum, at most the 30000 most intense peaks
  /// (FLASHDeconv's cap for centroided spectra); true if it was centroided here
  static bool toCentroidedPeakArrays(const MSSpectrum& spec, PeakArrays& peaks);

  /**
    @brief the charge that caps the fragment charges of a PrSM

    The charge of the raw precursor if @p precursor_mass has an isotope peak in the isolation window at that charge (or no mass is
    given, or the charges agree); else the charge implied by the mass and the isolation target m/z if that one does; else the
    larger of the two. Without a raw charge, the implied charge; without a mass either, 0 (no cap beyond min_mass_per_charge).
    @param reason 0 raw charge, 1 implied charge, 2 the larger of both, 3 implied charge (no raw charge), 4 no cap
  */
  int precursorChargeCap(const MSSpectrum& raw, double precursor_mass, int& reason);

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  /// the flattened averagine table: the kept isotope peaks of every 25 Da mass bin
  struct IsotopeTable
  {
    /// bins cover masses min_mass_ + b * bin_width_
    double max_mass = 0;
    /// begin of the kept peaks of a bin in iso_index / iso_rel; size = number of bins + 1
    std::vector<int> offset;
    /// isotope index k of a kept peak
    std::vector<int> iso_index;
    /// intensity relative to the most abundant peak
    std::vector<double> iso_rel;
    /// per bin: position (in the kept peaks of the bin) of the most abundant peak and of the monoisotopic peak (-1 if not kept)
    std::vector<int> apex_pos;
    std::vector<int> mono_pos;
    /// per bin: sum and sum of squares of the relative intensities (cosine and Pearson)
    std::vector<double> rel_sum;
    std::vector<double> rel_sum_sq;
  };

  /// match statistics of one envelope
  struct EnvelopeMatch
  {
    int matched = 0;
    int apex_peak = -1;
    /// indices of the matched observed peaks
    std::vector<int> peaks;
    double cosine = 0;
    double pearson = 0;
    double intensity = 0;
  };

  /// the 1/128 Th bins of a spectrum that hold a peak: a window without a peak is not searched (most windows hold none)
  struct PeakOccupancy
  {
    /// one bit per bin, from the bin of the first peak; empty: not built (every window may hold a peak)
    std::vector<UInt64> bits;
    double first_mz = 0;
    double last_mz = 0;
    /// floor(first_mz * bins_per_th)
    double first_bin = 0;
    static constexpr double bins_per_th = 128.0;
    void build(const PeakArrays& spec);
    /// false only if no peak of the spectrum lies in [low, high]
    bool mayHavePeak(double low, double high) const;
  };

  /// make sure the isotope table covers @p mass (computed in parallel, once per needed mass range)
  void ensureIsotopeTable_(double mass);
  /// averagine mass bin of a monoisotopic mass
  int massToBin_(double mass) const;
  /**
    Test the envelope of mass @p mono_mass at charge @p z with the isotope peaks of bin @p bin, the predicted m/z scaled by
    @p cal_factor (1 + calibration ppm * 1e-6). True if accepted; @p match is then filled. @p occupancy is that of @p spec.
  */
  bool matchEnvelope_(const PeakArrays& spec, const PeakOccupancy& occupancy, double mono_mass, int bin, int z, double cal_factor,
                      EnvelopeMatch& match) const;
  /// matchEnvelope_ once the window of the most abundant peak may hold a peak; @p base = mass + z * proton, @p scale = cal_factor / z
  bool searchEnvelope_(const PeakArrays& spec, const PeakOccupancy& occupancy, double base, double scale, int bin, EnvelopeMatch& match) const;
  /// the median ppm error of the calibration envelopes, or 0 if there are fewer than min_calibration_envelopes_
  double calibrate_(const PeakArrays& spec, const std::vector<Fragment>& fragments, int& count) const;
  /// annotate one PrSM (without the scores, which are computed afterwards); @p sorted_intensities: the spectrum's intensities, sorted
  void annotateOne_(const PeakArrays& spec, const std::vector<double>& sorted_intensities, const PrSM& prsm, Result& result) const;

  /**
    The internal pieces of a proteoform: residues [first, last] with 1 <= first < last <= L - 2 (both ends inside) and at most
    internal_max_length residues, for every internal ion type. Index (first, last, offset) = base[first] + (last - first - 1) *
    offsets + offset. A piece that splits a modification range, or has no charge in the spectrum's m/z range, has mono_mass -1.
  */
  struct InternalPieces
  {
    int length = 0;
    int offsets = 0;
    int max_length = 0;
    std::vector<int> base;
    std::vector<double> mono_mass;
    std::vector<int> max_charge;
    int index(int first, int last, int offset) const { return base[first] + (last - first - 1) * offsets + offset; }
    /// the last residue of the longest piece starting at first
    int lastMax(int first) const { return std::min(length - 2, first + max_length - 1); }
  };

  /// internal pieces of one world (the target, or the pieces moved by a decoy shift): which are accepted, and their runs
  struct InternalWorld
  {
    std::vector<char> accepted;
    /// per piece, the size of its run of consecutive lengths among the pieces that share its last / first residue
    std::vector<int> run_end, run_start;
    /// per cleavage a (after residue a), the longest run of pieces ending with a / starting right after a
    std::vector<int> site_end, site_start;
    /// envelopes of the accepted pieces (kept only on request)
    std::unordered_map<int, std::vector<Envelope>> envelopes;
  };

  /**
    Internal fragments of one world: every piece is tested at all its charges. An envelope is used if its most abundant peak is not
    a peak of an accepted terminal envelope (@p terminal_peaks) and lies within @p tol_ppm of the line
    @p line_offset + @p line_slope * m/z (ppm) that the terminal envelopes define. Then, for the cleavage after residue a and one
    internal ion type, the runs of accepted pieces of consecutive lengths that all end with a, or all start right after it.
  */
  void matchInternal_(const PeakArrays& spec, const PeakOccupancy& occupancy, const std::vector<double>& sorted_intensities, const InternalPieces& pieces,
                      const std::vector<char>& terminal_peaks, double cal_factor, double shift, double line_offset, double line_slope, double tol_ppm,
                      bool keep_envelopes, InternalWorld& world) const;

  /**
    The counted pieces of a world. The cleavage after residue a is authentic if one of its runs has at least @p site_min_run[a]
    pieces (a tag); the pieces of such runs count. Beyond them, accepted pieces with an end at an authentic cleavage (of this world
    or given in @p given_authentic: a decoy world gets the target's, so that chance matches of the extension at real cleavages are
    measured too) count if they lie in a run of at least @p extension_min_run pieces there (1: every piece, 2: a tag of length one;
    0: none). Returns the counted pieces as (first, last); their envelopes are added to @p envelopes if @p record; @p authentic
    (one flag per residue) receives the authentic cleavages of this world.
  */
  std::set<std::pair<int, int>> countInternal_(const InternalPieces& pieces, const InternalWorld& world, const std::vector<int>& site_min_run,
                                               int extension_min_run, const std::vector<char>& given_authentic, bool record,
                                               std::vector<Envelope>& envelopes, std::vector<char>& authentic) const;

  /**
    Tag lengths from the decoy shifts (internal_alpha): for the cleavages after Asp (0), before Pro (1) and elsewhere (2), the
    smallest number of pieces k at which the decoy worlds in @p choosing give at most alpha of the cleavages the target has with a
    run of at least k pieces. Beyond the longest run the decoys show, their count is continued geometrically. 0: no length passes.
  */
  std::array<int, 3> chooseTagLengths_(const std::string& sequence, const std::vector<InternalWorld>& worlds, const std::vector<Size>& choosing) const;

  /// matched peaks an envelope with @p predicted kept isotope peaks needs; -1 if it cannot be annotated
  int requiredPeaks_(int predicted) const;
  /// -log10 of the calibrated p-value of a RawScore
  double calibratedTerminalLog10P_(double score) const;
  /// -log10 P(K >= k | n) of the internal evidence, from the calibration model
  double internalLog10P_(int k, int n) const;

  /// observed peaks per Th within 1 Th of @p mz
  static double densityAt_(const PeakArrays& spec, double mz);
  /// fraction of @p sorted_intensities not above @p intensity
  static double rankOf_(const std::vector<double>& sorted_intensities, double intensity);

  /// parameters
  double tol_ = 10.0;
  double cal_tol_ = 20.0;
  bool calibrate_on_ = true;
  int min_peaks_ = 3;
  int min_peaks_small_ = 2;
  std::vector<double> terminal_calibration_ {0.479, 0.697};
  std::vector<double> internal_calibration_ {-11.531, 2.062, -0.471};
  double min_rel_iso_ = 0.1;
  double min_cos_ = 0.8;
  double min_pearson_ = 0.7;
  double min_mass_per_charge_ = 300.0;
  double max_mass_ = 200000.0;
  double p1_ = 0.3;
  bool report_decoys_ = false;
  bool internal_on_ = true;
  int internal_max_length_ = 150;
  /// pieces in a ladder that make a cleavage authentic: after Asp, elsewhere, and before Pro (tag length + 1)
  int min_run_asp_ = 4;
  int min_run_other_ = 6;
  int min_run_pro_ = 5;
  /// charges per kDa of the precursor above which the thresholds after Asp and before Pro are exchanged (0: never)
  double charge_switch_ = 0.45;
  /// only internal pieces of the b/y type (two collisional cleavages)
  bool internal_by_only_ = true;
  /// share of chance cleavages allowed when the tag lengths are taken from the decoy shifts (0: fixed lengths)
  double internal_alpha_ = 0.05;
  /// tolerance (ppm) of internal pieces: 0 automatic (from the terminal envelopes), > 0 fixed, < 0 the input tolerance without a line
  double internal_tol_ = 0;
  /// minimum run of consecutive lengths for pieces beyond the ladders: 0 none, 1 all, 2 a tag of length one
  int extension_min_run_ = 0;
  /// internal ion types: mass added to the residue sum (C-end prefix shift + N-end suffix shift - H2O), and their names
  std::vector<double> internal_offsets_;
  std::vector<std::string> internal_ion_types_;
  std::vector<double> decoy_shifts_;
  std::vector<std::string> ion_types_;
  /// mass added to the residue sum of a fragment, per ion type; and whether the ion type is a prefix (N-terminal) one
  std::vector<double> ion_shifts_;
  std::vector<bool> ion_is_prefix_;
  /// used for calibration (b and y if present, else all)
  std::vector<bool> ion_for_calibration_;
  /// per ion type: its group of ion types whose envelopes fall on the same peaks (z, z+1, z+2), and the number of groups
  std::vector<int> ion_group_;
  int ion_group_count_ = 0;

  IsotopeTable table_;
  double table_min_rel_iso_ = -1;

  /// averagine bins as in SpectralDeconvolution::calculateAveragine
  static constexpr double min_mass_ = 50.0;
  static constexpr double bin_width_ = 25.0;
  /// fewer calibration envelopes than this: no correction
  static constexpr int min_calibration_envelopes_ = 5;
  /// cap of the score
  static constexpr double max_score_ = 1000.0;
  /// ion types of the same terminus whose shifts differ by 1 .. max_isotope_offset_ isotope spacings within this mass (Da)
  /// match the same envelopes (z, z+1 and z+2 differ by 4.5 mDa per spacing) and are one fragment hypothesis
  static constexpr double isotope_equivalence_tol_ = 0.01;
  static constexpr int max_isotope_offset_ = 3;
};
} // namespace OpenMS
