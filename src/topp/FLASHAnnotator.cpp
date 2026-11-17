// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
@page TOPP_FLASHAnnotator FLASHAnnotator

  @brief FLASHAnnotator annotates the PrSMs reported by FLASHTnT on the raw MS2 spectra.

  PrSMs of other tools are read as ProForma 2.0 strings: a tsv file with the columns Scan and ProForma (and optionally
  PrecursorMass) and no ProteinSequence column. Modifications are given by name, accession, formula or mass, on residues, ranges or
  termini; one at an unknown position spans the sequence and a labile one is left out of the fragments (see ProForma for what is
  supported). All other columns are written unchanged.

  For each PrSM of a FLASHTnT PrSM file (-out_prsm of FLASHTnT), the theoretical isotope envelope of every terminal
  fragment of the reported proteoform is placed at every plausible charge on the centroided raw MS2 spectrum (profile
  spectra are centroided with PeakPickerHiRes), the way ExD Viewer annotates, and the envelopes that match are kept
  (see FLASHAnnotatorAlgorithm). FLASHAnnotator reports the bond coverage and a PrSM score that is calibrated against
  the chance match rate of mass-shifted decoy fragments in the same spectrum, so that it does not grow with the
  proteoform length.

  The PrSM rows are written unchanged, followed by the columns RawScore (-log10 of the binomial tail probability of the
  matched fragment hypotheses), RawLLR (log-likelihood ratio), RawMatchedFragments, RawFragmentHypotheses, RawChanceRate,
  RawCoveredBonds, RawCoverage(%) and RawCalibrationPPM. -1 marks values that could not be computed (scan not found,
  unknown residue; for RawScore and RawChanceRate also no decoy shifts or no fragment hypothesis); RawLLR and
  RawCalibrationPPM, which can be negative, are NaN then. RawFragmentHypotheses = -1 marks a PrSM that was not annotated.
  -out_annotation lists the accepted envelopes.

  RawScore and RawLLR are evidence scores of a PrSM in its own spectrum, not identification p-values: FLASHTnT selects
  its candidates by tags and blind modifications fitted to the spectrum, so even a wrong candidate matches more fragments
  than a random proteoform. For identification, use target-decoy competition (e.g. the scores of the PrSMs of decoy
  proteins, written by FLASHTnT with -keep_decoy).

  Fragment charges are capped by the precursor charge: the charge of the raw precursor if the precursor mass of the PrSM
  (PrecursorMass) puts an isotope peak in the isolation window at that charge; otherwise the charge implied by the mass
  and the isolation target m/z if that charge does, else the larger of the two. Without a raw charge, the implied charge
  is used; without a precursor mass either, the fragment charges are capped only by algorithm:min_mass_per_charge. The
  proteoform mass is not used: it is not the precursor mass when the precursor is a complex.

  The input mzML must be the raw file FLASHDeconv was run on (FLASHTnT scan numbers refer to it). Only the MS2 spectra of
  the PrSM scans are kept in memory.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FLASHAnnotator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FLASHAnnotator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFLASHAnnotator : public TOPPBase
{
public:
  TOPPFLASHAnnotator():
      TOPPBase("FLASHAnnotator",
               "Annotates FLASHTnT PrSMs on the raw MS2 spectra with fragment isotope envelopes and scores them against the chance "
               "match rate.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "",
                       "Input raw mzML file (profile or centroid) the PrSMs were identified from, or the Thermo .raw file (recognised by its "
                       "extension; read with Thermo's peak picking, so centroided, as ThermoRawFileParser converts it by default; needs OpenMS "
                       "built with WITH_THERMO_RAW and the .NET 8 runtime).", true);
    setValidFormats_("in", ListUtils::create<std::string>("mzML,raw"));

    registerInputFile_("in_prsm", "<file>", "", "Input PrSM tsv file: -out_prsm of FLASHTnT, or any tsv with the columns Scan and ProForma (optionally PrecursorMass).", true);
    setValidFormats_("in_prsm", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out", "<file>", "", "Output PrSM tsv file: the input rows followed by the raw annotation columns.", true);
    setValidFormats_("out", ListUtils::create<std::string>("tsv"));

    registerOutputFile_("out_annotation", "<file>", "", "Output tsv file with one row per accepted fragment isotope envelope.", false);
    setValidFormats_("out_annotation", ListUtils::create<std::string>("tsv"));

    registerStringList_("ion_type", "<types>", ListUtils::toStringList(std::vector<std::string> {"b", "y"}),
                        "Ion types to annotate (as for FLASHTnT).", false);
    setValidStrings_("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

    registerIntOption_("precursor_charge", "<charge>", 0,
                       "Precursor charge used for all PrSMs (caps the fragment charges). 0: the charge of the raw precursor if "
                       "the PrecursorMass of the PrSM has an isotope peak in the isolation window at it; else "
                       "round(PrecursorMass / (precursor m/z - proton)) if that charge does, or if the "
                       "raw charge is not annotated; else the larger of the two; without a mass or precursor m/z, no cap other "
                       "than algorithm:min_mass_per_charge.",
                       false);
    setMinInt_("precursor_charge", 0);

    registerSubsection_("algorithm", "Annotation algorithm parameters");
  }

  Param getSubsectionDefaults_(const std::string& prefix) const override
  {
    if (prefix == "algorithm")
    {
      auto param = FLASHAnnotatorAlgorithm().getDefaults();
      param.remove("ion_type");
      return param;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix);
  }

  static double toDouble_(std::string s, double default_value)
  {
    try
    {
      StringUtils::trim(s);
      return s.empty() ? default_value : StringUtils::toDouble(s);
    }
    catch (Exception::BaseException&)
    {
      return default_value;
    }
  }

  static int toInt_(std::string s, int default_value)
  {
    try
    {
      StringUtils::trim(s);
      return s.empty() ? default_value : (int)std::lround(StringUtils::toDouble(s));
    }
    catch (Exception::BaseException&)
    {
      return default_value;
    }
  }

  /// ';'-separated values, empty ones skipped
  static std::vector<std::string> splitList_(const std::string& s)
  {
    std::vector<std::string> tokens, values;
    StringUtils::split(s, ";", tokens);
    for (auto& t : tokens)
    {
      StringUtils::trim(t);
      if (! t.empty()) values.push_back(t);
    }
    return values;
  }

  static double seconds_(const std::chrono::steady_clock::time_point& since)
  {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
  }

  /// @p value as an ostream writes it, or "NaN"
  static std::string format_(double value)
  {
    if (std::isnan(value)) return "NaN";
    std::ostringstream ss;
    ss << value;
    return ss.str();
  }

  /// the appended PrSM columns and the envelope header
  static constexpr const char* prsm_columns_ = "\tFAScore\tRawScore\tRawLLR\tRawMatchedFragments\tRawFragmentHypotheses\tRawChanceRate\tRawCoveredBonds"
                                                "\tRawCoverage(%)\tRawCalibrationPPM\tRawInternalScore\tRawInternalMatched\tRawInternalHypotheses\tRawInternalChanceRate"
                                                "\tRawAuthenticCleavages\tRawCoverageWithInternal(%)\tRawInternalChanceBonds\tRawInternalTagPieces\tRawInternalTol";
  static constexpr const char* envelope_header_ = "PrSMIndex\tScan\tIonType\tLength\tBond\tCharge\tMonoMass\tMostAbundantMz\tMatchedPeaks"
                                                   "\tPredictedPeaks\tCosine\tPearson\tIntensity\tPPMError\tIntensityRank\tLocalDensity\tDecoyShift\tInternalFirst";

  ExitCodes main_(int, const char**) override
  {
    const auto start_time = std::chrono::steady_clock::now();
    std::string in_file = getStringOption_("in");
    std::string in_prsm_file = getStringOption_("in_prsm");
    std::string out_file = getStringOption_("out");
    std::string out_annotation_file = getStringOption_("out_annotation");
    const int fixed_precursor_charge = getIntOption_("precursor_charge");

    Param algorithm_param = getParam_().copy("algorithm:", true);
    algorithm_param.setValue("ion_type", getParam_().getValue("ion_type"));
    FLASHAnnotatorAlgorithm annotator;
    annotator.setParameters(algorithm_param);
    const bool detection = algorithm_param.getValue("ion_type_detection").toString() == "true";

    //-------------------------------------------------------------
    // PrSMs: columns are found by their header names
    //-------------------------------------------------------------
    auto time = std::chrono::steady_clock::now();
    std::ifstream prsm_stream(in_prsm_file);
    if (! prsm_stream) { throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in_prsm_file); }
    std::string header;
    std::vector<std::string> lines;
    std::vector<std::vector<std::string>> rows;
    {
      std::string line;
      bool first = true;
      while (std::getline(prsm_stream, line))
      {
        if (! line.empty() && line.back() == '\r') line.pop_back();
        if (first)
        {
          header = line;
          first = false;
          continue;
        }
        if (StringUtils::trimmed(line).empty()) continue;
        lines.emplace_back(line);
        std::vector<std::string> fields;
        StringUtils::split(lines.back(), '\t', fields);
        rows.push_back(fields);
      }
    }
    std::vector<std::string> columns;
    StringUtils::split(header, '\t', columns);
    auto column = [&columns](const std::string& name) -> int {
      for (Size i = 0; i < columns.size(); i++)
      {
        if (StringUtils::trimmed(columns[i]) == name) return (int)i;
      }
      return -1;
    };
    int col_index = column("PrSMIndex");
    if (col_index < 0) col_index = 0;
    const int col_scan = column("Scan");
    const int col_sequence = column("ProteinSequence");
    const int col_start = column("StartPosition");
    const int col_end = column("EndPosition");
    const int col_precursor_mass = column("PrecursorMass");
    const int col_mod_mass = column("ModMass");
    const int col_mod_start = column("ModStart");
    const int col_mod_end = column("ModEnd");
    const int col_database = column("DatabaseSequence");
    // a PrSM file of another tool: one ProForma string per PrSM (column ProForma), read when there is no ProteinSequence column
    int col_proforma = column("ProForma");
    if (col_proforma < 0) col_proforma = column("Proforma");
    const bool use_proforma = col_sequence < 0 && col_proforma >= 0;
    if (col_scan < 0 || (col_sequence < 0 && ! use_proforma))
    {
      OPENMS_LOG_ERROR << "The PrSM file needs a Scan column and either a ProteinSequence (FLASHTnT) or a ProForma column: " << in_prsm_file << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    auto field = [](const std::vector<std::string>& row, int col) -> std::string { return (col >= 0 && col < (int)row.size()) ? row[col] : std::string(); };

    std::set<int> needed_scans;
    for (const auto& row : rows)
    {
      needed_scans.insert(toInt_(field(row, col_scan), -1));
    }
    const double read_time = seconds_(time);

    if (rows.empty())
    {
      // nothing to annotate: the raw file is not read
      std::ofstream out(out_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file); }
      out << header << prsm_columns_ << (detection ? "\tDetectedIonTypes" : "") << "\n";
      if (! out_annotation_file.empty())
      {
        std::ofstream out_annotation(out_annotation_file);
        if (! out_annotation) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_annotation_file); }
        out_annotation << envelope_header_ << "\n";
      }
      OPENMS_LOG_INFO << "PrSMs: 0" << endl;
      return EXECUTION_OK;
    }

    //-------------------------------------------------------------
    // raw spectra: only the MS2 spectra of the PrSM scans are kept
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    std::map<int, MSSpectrum> raw_map;
    const Size read_count = FLASHAnnotatorAlgorithm::loadSpectra(in_file, needed_scans, raw_map);
    std::vector<MSSpectrum> raw_spectra;
    std::map<int, int> scan_to_slot;
    raw_spectra.reserve(raw_map.size());
    for (auto& [scan, spec] : raw_map)
    {
      scan_to_slot[scan] = (int)raw_spectra.size();
      raw_spectra.push_back(std::move(spec));
    }
    raw_map.clear();
    const double load_time = seconds_(time);

    //-------------------------------------------------------------
    // the spectra of the PrSMs: centroided (if profile) once per spectrum, in parallel; the raw peaks are then released
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(raw_spectra.size());
    int profile_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : profile_count)
    for (int s = 0; s < (int)raw_spectra.size(); s++)
    {
      MSSpectrum& spec = raw_spectra[s];
      profile_count += FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(spec, spectra[s]) ? 1 : 0;
      spec.clear(false);
      spec.shrink_to_fit();
    }

    // PrSM inputs
    std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(rows.size());
    std::vector<double> prsm_masses(rows.size(), -1.0);
    double max_prsm_mass = 0;
    int missing_scan = 0, invalid_proforma = 0;
    for (Size r = 0; r < rows.size(); r++)
    {
      const auto& row = rows[r];
      auto& prsm = prsms[r];
      if (use_proforma)
      {
        // a ProForma proteoform: both termini known; a terminal modification sits on the terminal residue, one at an unknown position
        // spans the sequence, and a labile one is not on the fragments
        try
        {
          using Site = FLASHHelperClasses::ProFormaModification::Site;
          const ProForma::PeptidoformIon proforma = FLASHHelperClasses::parseProForma(StringUtils::trimmed(field(row, col_proforma)));
          std::string sequence;
          const auto modifications = FLASHHelperClasses::getProFormaModifications(proforma, sequence);
          prsm.sequence = sequence;
          const int last = (int)prsm.sequence.size() - 1;
          for (const auto& m : modifications)
          {
            FLASHAnnotatorAlgorithm::Modification mod;
            mod.mass = m.mass;
            if (m.site == Site::LABILE) continue;
            if (m.site == Site::RESIDUE)
            {
              mod.start = m.start;
              mod.end = m.end;
            }
            else
            {
              mod.start = m.site == Site::C_TERM ? last : 0;
              mod.end = m.site == Site::N_TERM ? 0 : last;
            }
            prsm.modifications.push_back(mod);
          }
        }
        catch (const Exception::BaseException& e)
        {
          invalid_proforma++;
          OPENMS_LOG_WARN << "PrSM " << (r + 1) << ": " << e.what() << endl;
          continue;
        }
      }
      else
      {
      prsm.sequence = StringUtils::trimmed(field(row, col_sequence));
      const int start = toInt_(field(row, col_start), 1);
      const int end = toInt_(field(row, col_end), 1);
      // FLASHTnT writes -1 for a start or end it could not determine: those fragments have no defined mass
      prsm.known_n_term = start > 0;
      prsm.known_c_term = end > 0;
      // modification positions are 1-based protein positions. With an open N-terminus the sequence ends at the known end, or (both
      // termini open) lies in the database sequence: it starts after residue 1 when FLASHTnT bounds it by a PEFF entry's processed form
      int start_in_seq = start < 0 ? 0 : start - 1;
      if (start < 0 && end > 0) start_in_seq = std::max(0, end - (int)prsm.sequence.size());
      else if (start < 0 && col_database >= 0)
      {
        const Size at = StringUtils::trimmed(field(row, col_database)).find(prsm.sequence);
        if (at != std::string::npos) start_in_seq = (int)at;
      }
      const auto mod_masses = splitList_(field(row, col_mod_mass));
      const auto mod_starts = splitList_(field(row, col_mod_start));
      const auto mod_ends = splitList_(field(row, col_mod_end));
      for (Size m = 0; m < mod_masses.size() && m < mod_starts.size() && m < mod_ends.size(); m++)
      {
        FLASHAnnotatorAlgorithm::Modification mod;
        mod.mass = toDouble_(mod_masses[m], 0);
        mod.start = toInt_(mod_starts[m], 1) - 1 - start_in_seq;
        mod.end = toInt_(mod_ends[m], 1) - 1 - start_in_seq;
        prsm.modifications.push_back(mod);
      }
      }

      // the precursor mass only: the proteoform mass is not the precursor mass when the precursor is a complex or a
      // co-isolated species, and a charge derived from it would cap the fragment charges too low
      const double mass = toDouble_(field(row, col_precursor_mass), -1);
      prsm_masses[r] = mass;
      max_prsm_mass = std::max({max_prsm_mass, mass, FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});

      const auto it = scan_to_slot.find(toInt_(field(row, col_scan), -1));
      if (it == scan_to_slot.end())
      {
        missing_scan++;
        continue;
      }
      prsm.spectrum_index = it->second;
    }

    // precursor charges. The isotope table is computed once, beyond the heaviest PrSM (the annotation then reuses it).
    int charge_reasons[5] = {0, 0, 0, 0, 0};
    std::vector<double> iso_mzs, iso_intensities;
    if (fixed_precursor_charge == 0) annotator.getIsotopeEnvelope(max_prsm_mass + 200.0, 1, iso_mzs, iso_intensities);
    for (Size r = 0; r < rows.size(); r++)
    {
      auto& prsm = prsms[r];
      if (prsm.spectrum_index < 0) continue;
      // the isolation target m/z: its charge density selects the internal-ion thresholds
      const auto& precursors = raw_spectra[prsm.spectrum_index].getPrecursors();
      if (! precursors.empty()) prsm.precursor_mz = precursors[0].getMZ();
      if (fixed_precursor_charge > 0)
      {
        prsm.precursor_charge = fixed_precursor_charge;
        continue;
      }
      int reason = 4;
      prsm.precursor_charge = annotator.precursorChargeCap(raw_spectra[prsm.spectrum_index], prsm_masses[r], reason);
      charge_reasons[reason]++;
    }
    const double centroid_time = seconds_(time);

    //-------------------------------------------------------------
    // annotation
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    std::vector<FLASHAnnotatorAlgorithm::Result> results;
    annotator.annotate(spectra, prsms, results);
    const double annotation_time = seconds_(time);

    //-------------------------------------------------------------
    // output, in the input order
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    {
      std::ofstream out(out_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file); }
      out << header << prsm_columns_ << (detection ? "\tDetectedIonTypes" : "") << "\n";
      for (Size r = 0; r < rows.size(); r++)
      {
        const auto& res = results[r];
        out << lines[r] << "\t" << res.fa_score << "\t" << res.score << "\t" << format_(res.llr) << "\t" << res.matched_fragments << "\t" << res.fragment_hypotheses
            << "\t" << res.chance_rate << "\t" << res.covered_bonds << "\t" << res.coverage << "\t" << format_(res.calibration_ppm) << "\t"
            << res.internal_score << "\t" << res.internal_matched << "\t" << res.internal_hypotheses << "\t" << res.internal_chance_rate << "\t"
            << res.authentic_cleavages << "\t" << res.coverage_with_internal << "\t" << res.internal_chance_bonds << "\t" << res.internal_tag_pieces_asp << "/"
            << res.internal_tag_pieces_pro << "/" << res.internal_tag_pieces_other << "\t" << res.internal_tol_ppm;
        if (detection) out << "\t" << res.detected_ion_types;
        out << "\n";
      }
    }
    if (! out_annotation_file.empty())
    {
      std::ofstream out(out_annotation_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_annotation_file); }
      out << envelope_header_ << "\n";
      const auto& ion_types = annotator.getIonTypes();
      const auto& internal_ion_types = annotator.getInternalIonTypes();
      for (Size r = 0; r < rows.size(); r++)
      {
        const std::string prsm_index = field(rows[r], col_index);
        const std::string scan = field(rows[r], col_scan);
        std::ostringstream ss;
        for (const auto& e : results[r].envelopes)
        {
          ss << prsm_index << "\t" << scan << "\t" << (e.first >= 0 ? std::string("int(") + internal_ion_types[e.ion_type] + ")" : ion_types[e.ion_type]) << "\t" << e.length << "\t" << e.bond << "\t" << e.charge << "\t"
             << std::fixed << std::setprecision(5) << e.mono_mass << "\t" << e.most_abundant_mz << "\t" << e.matched_peaks << "\t"
             << e.predicted_peaks << "\t" << std::setprecision(4) << e.cosine << "\t" << e.pearson << "\t" << std::setprecision(1) << e.intensity
             << "\t" << std::setprecision(2) << e.ppm_error << "\t" << std::setprecision(4) << e.intensity_rank << "\t"
             << std::setprecision(1) << e.local_density << "\t" << std::setprecision(3) << e.decoy_shift << "\t" << (e.first + 1) << "\n";
        }
        out << ss.str();
      }
    }
    const double write_time = seconds_(time);

    Size annotated = 0, envelopes = 0;
    for (const auto& res : results)
    {
      annotated += res.annotated ? 1 : 0;
      envelopes += res.envelopes.size();
    }
    OPENMS_LOG_INFO << "PrSMs: " << rows.size() << " (" << annotated << " annotated, " << missing_scan << " without their scan in the mzML, "
                    << invalid_proforma << " with an invalid or unsupported ProForma string, " << (rows.size() - annotated - missing_scan - invalid_proforma)
                    << " with an unknown residue); accepted envelopes: " << envelopes << endl;
    OPENMS_LOG_INFO << "Spectra: " << read_count << " read, " << spectra.size() << " kept and annotated (" << profile_count << " centroided)" << endl;
    if (fixed_precursor_charge == 0)
    {
      OPENMS_LOG_INFO << "Precursor charge: " << charge_reasons[0] << " from the raw precursor, " << charge_reasons[1]
                      << " from the precursor mass and the isolation window (raw charge does not fit), " << charge_reasons[2]
                      << " the larger of both (neither fits), " << charge_reasons[3] << " from the precursor mass (no raw charge), "
                      << charge_reasons[4] << " without a precursor cap" << endl;
    }
    OPENMS_LOG_INFO << std::fixed << std::setprecision(2) << "Timing (s): PrSM file " << read_time << ", mzML load " << load_time
                    << ", centroiding " << centroid_time << ", annotation " << annotation_time << ", writing " << write_time
                    << ", total " << seconds_(start_time) << endl;
    return EXECUTION_OK;
  }
};

/// @endcond

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHAnnotator tool;
  return tool.main(argc, argv);
}
