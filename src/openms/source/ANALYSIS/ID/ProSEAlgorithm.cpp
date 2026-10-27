// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#endif
// Not guarded: header-only constants with no ONNX dependency, needed to register the same
// 'peptdeep:instrument' values in either build.
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

#include <algorithm>
#if defined(_OPENMP) && defined(__GLIBCXX__) // libstdc++'s parallel mode (not libc++)
#include <parallel/algorithm> // __gnu_parallel::nth_element (medianAndMADBySelection)
#define PROSE_GNU_PARALLEL
#endif
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <limits>
#include <numeric>
#include <iomanip>
#include <iterator>
#include <locale>
#include <ostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

#if defined(__linux__)
  #include <cstdint>
  #include <sys/mman.h> // madvise (resizeWithParallelFaults)
  #include <unistd.h>   // sysconf
#endif

using namespace std;

namespace OpenMS
{
  namespace
  {
    // Destroys the elements of a large container with the OpenMP threads, then the container itself. Destroying the
    // spectra, hits and identifications of a search one after the other is a serial step of 10 - 100 ms each at the
    // end of a search; the elements own their memory, so they can be destroyed in any order.
    // In equal static parts per thread: glibc returns a block to the arena of the thread that allocated it, under that
    // arena's lock. The peak lists were made by the static loop of preprocessSpectra_ (the same parts of the same
    // team), so each thread frees into its own arena; in dynamic chunks most frees waited for the lock of another
    // thread's arena (futex system calls: HeLa 25pc_50_25 at 64 threads 0.44 s for the query spectra instead of 0.06).
    // (Only where the memory goes back depends on this: with another team size or a container changed in between,
    // some frees cross arenas, with the same result.)
    template <typename Container>
    void freeInParallel(Container& container)
    {
      const SignedSize n = static_cast<SignedSize>(container.size());
#pragma omp parallel for schedule(static) if (n >= 256)
      for (SignedSize i = 0; i < n; ++i)
      {
        [[maybe_unused]] const auto released = std::move(container[i]);
      }
      Container().swap(container);
    }

    // Destroys the elements of a large container on the OpenMP threads that allocated them (owner(i): the thread
    // number of element i in a team of the same size; other numbers are spread over the team), then the container:
    // as freeInParallel(), for elements that a dynamically scheduled loop made. Each thread looks at all owners and
    // frees its own elements, into its own arena.
    template <typename Container, typename Owner>
    void freeOnOwnerThreads(Container& container, const Owner& owner)
    {
      const SignedSize n = static_cast<SignedSize>(container.size());
#pragma omp parallel if (n >= 256)
      {
#ifdef _OPENMP
        const int team = omp_get_num_threads();
        const int me = omp_get_thread_num();
#else
        const int team = 1;
        const int me = 0;
#endif
        for (SignedSize i = 0; i < n; ++i)
        {
          const int o = owner(i);
          if ((o >= 0 ? o % team : 0) == me) { [[maybe_unused]] const auto released = std::move(container[i]); }
        }
      }
      Container().swap(container);
    }

    // The peak list and native ID of @p from into the empty spectrum @p to (and, if @p with_precursors, its MS level and
    // precursors): a full MSSpectrum copy also copies the meta values, the data arrays and the other settings, heap blocks
    // per spectrum that the readers of these copies never look at (the native ID identifies the copy).
    void copyPeakList_(const MSSpectrum& from, MSSpectrum& to, bool with_precursors)
    {
      to.reserve(from.size());
      to.insert(to.end(), from.begin(), from.end());
      to.setNativeID(from.getNativeID());
      if (with_precursors)
      {
        to.setMSLevel(from.getMSLevel());
        to.setPrecursors(from.getPrecursors());
      }
    }

    void freeInParallel(PeakMap& spectra)
    {
      freeInParallel(spectra.getSpectra());
      spectra.clear(true);
    }

    // Resizes @p v to @p n elements as resize() does; on Linux, the pages of the storage that it allocates for them are
    // first faulted in by the OpenMP threads (madvise(MADV_POPULATE_WRITE) on parts of the range, Linux 5.14 or later; no
    // access to the storage itself), so that the serial construction of the elements does not take the page faults of
    // the whole block one after the other (619,000 empty spectra, 465 MB: 0.21 s, then 0.015 s for the faults in parallel
    // and 0.07 s for the construction at 64 threads). Elsewhere, or if the call fails, as resize().
    template <typename T>
    void resizeWithParallelFaults(std::vector<T>& v, const Size n)
    {
#if defined(__linux__)
      if (n > v.capacity())
      {
        v.reserve(n);
        const uintptr_t page = static_cast<uintptr_t>(sysconf(_SC_PAGESIZE));
        const uintptr_t begin = (reinterpret_cast<uintptr_t>(v.data() + v.size()) + page - 1) / page * page;
        const uintptr_t end = reinterpret_cast<uintptr_t>(v.data() + n) / page * page; // whole pages of the new elements
        const SignedSize pages = end > begin ? static_cast<SignedSize>((end - begin) / page) : 0;
        if (pages >= 1024)
        {
          constexpr int populate_write = 23; // MADV_POPULATE_WRITE (EINVAL before Linux 5.14: then the faults come later)
          const SignedSize parts = 256;
#pragma omp parallel for schedule(dynamic, 1)
          for (SignedSize part = 0; part < parts; ++part)
          {
            const uintptr_t from = begin + static_cast<uintptr_t>(pages * part / parts) * page;
            const uintptr_t to = begin + static_cast<uintptr_t>(pages * (part + 1) / parts) * page;
            if (to > from) { madvise(reinterpret_cast<void*>(from), to - from, populate_write); }
          }
        }
      }
#endif
      v.resize(n);
    }

    // The positions 0 .. n - 1 of @p keys (no NaN) ordered by (key, position) ascending, as MSExperiment::sortSpectra()
    // orders spectra by retention time: a strict total order, so the order does not depend on how it is sorted. Parts
    // sorted by the OpenMP threads, then merged pairwise.
    std::vector<Size> orderByKeyInParallel(const std::vector<double>& keys)
    {
      const Size n = keys.size();
      std::vector<Size> order(n);
      std::iota(order.begin(), order.end(), Size(0));
      const auto less = [&keys](Size a, Size b) { return keys[a] < keys[b] || (!(keys[b] < keys[a]) && a < b); };
      int parts = 1;
#ifdef _OPENMP
      parts = static_cast<int>(std::max<Size>(1, std::min<Size>(static_cast<Size>(omp_get_max_threads()), n / 4096)));
#endif
      std::vector<Size> bounds(static_cast<Size>(parts) + 1);
      for (int p = 0; p <= parts; ++p) { bounds[p] = n / static_cast<Size>(parts) * static_cast<Size>(p) + std::min(n % static_cast<Size>(parts), static_cast<Size>(p)); }
#pragma omp parallel for schedule(static, 1)
      for (int p = 0; p < parts; ++p) { std::sort(order.begin() + bounds[p], order.begin() + bounds[p + 1], less); }
      std::vector<Size> merged(n);
      for (int width = 1; width < parts; width *= 2)
      {
        const int pairs = (parts + 2 * width - 1) / (2 * width);
#pragma omp parallel for schedule(dynamic, 1)
        for (int q = 0; q < pairs; ++q)
        {
          const int p = 2 * width * q;
          const Size lo = bounds[p], mid = bounds[std::min(p + width, parts)], hi = bounds[std::min(p + 2 * width, parts)];
          std::merge(order.begin() + lo, order.begin() + mid, order.begin() + mid, order.begin() + hi, merged.begin() + lo, less);
        }
        order.swap(merged);
      }
      return order;
    }

    // Whether MSExperiment::sortSpectra(false) would leave @p spectra as they are: their retention times ascend (equal ones
    // keep their order in its sort) and none is NaN (which would make the order depend on the sort). An O(n) check that
    // saves the O(n log n) sort (0.1-0.2 s for the 250,000 spectra of a timsTOF run).
    bool inRetentionTimeOrder(const PeakMap& spectra)
    {
      const MSSpectrum::RTLess rt_less;
      return std::none_of(spectra.begin(), spectra.end(), [](const MSSpectrum& spectrum) { return std::isnan(spectrum.getRT()); })
             && std::adjacent_find(spectra.begin(), spectra.end(),
                                   [&rt_less](const MSSpectrum& a, const MSSpectrum& b) { return rt_less(b, a); }) == spectra.end();
    }
  }

  namespace
  {
    /// The m/z of the monoisotopic precursor peak of an observed precursor m/z that FragmentIndex matched at the isotope
    /// offset @p isotope_offset (13C spacings, theoretical minus observed: it searches shifted_mass = precursor_mass +
    /// isotope_offset * C13C12, so M_theo ~ M_obs + isotope_offset * C13C12): observed + offset * C13C12 / charge. The
    /// precursor errors of the PSMs (PRECURSOR_ERROR_PPM_USERPARAM), of the calibration and of the run statistics are
    /// those of this m/z.
    double monoisotopicPrecursorMZ(double observed_mz, int isotope_offset, double charge)
    {
      return observed_mz + static_cast<double>(isotope_offset) * Constants::C13C12_MASSDIFF_U / charge;
    }
  }

  ProSEAlgorithm::ProSEAlgorithm() :
    DefaultParamHandler("ProSEAlgorithm"),
    ProgressLogger()
  {
    defaults_.setValue("precursor:mass_tolerance_lower", 10.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 10.0,
                       "Upper-side precursor-mass tolerance (positive magnitude).");
    defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", {"ppm", "Da"});

    defaults_.setValue("precursor:min_charge", 2,
                       "Minimum precursor charge to be considered. Spectra without a precursor charge are searched with every "
                       "charge of the range; spectra whose recorded charge lies outside of it are searched at that charge, or "
                       "not at all with precursor:recorded_charge_outside_range 'skip'.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered (see precursor:min_charge).");
    defaults_.setValue("precursor:recorded_charge_outside_range", "search",
                       "A spectrum whose recorded precursor charge lies outside of precursor:min_charge..max_charge: 'search' "
                       "(ProSEAlgorithm's default, as ProSE int15) searches it at its recorded charge; 'skip' (the ProSE tool's "
                       "default; OpenMS develop since #10713) "
                       "does not search it. Spectra without a recorded charge are searched with every charge of the range.",
                       {"advanced"});
    defaults_.setValidStrings("precursor:recorded_charge_outside_range", {"search", "skip"});

    defaults_.setSectionDescription("precursor",
      "Precursor (Parent Ion) Options. mass_tolerance_lower/_upper are positive magnitudes "
      "applied as [-lower, +upper] around the precursor mass.");

    defaults_.setValue("fragment:mass_tolerance", 20.0, "Fragment mass tolerance");

    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    defaults_.setValue("fragment:deisotope", "auto", "MS2 deisotoping (single-charge deconvolution) before searching. 'auto' deisotopes only when the fragment tolerance is within the high-resolution deisotoper range (<= 0.1 Da / <= 100 ppm) and skips it for low-resolution (e.g. ion-trap CID) data; 'true' always deisotopes (requires a high-resolution fragment tolerance); 'false' never deisotopes.");
    defaults_.setValidStrings("fragment:deisotope", {"auto", "true", "false"});
    defaults_.setValue("fragment:deisotope_min_peaks", 2,
                       "Minimum number of peaks, the monoisotopic peak included, of an isotope envelope that MS2 deisotoping "
                       "collapses into its monoisotopic peak. With 2 (Sage-like), an ion whose M+2 peak is lost in the noise loses its "
                       "M+1 peak too, so the M+1 peak no longer takes a place in the per-window peak quota, and a multiply charged ion "
                       "seen with two isotope peaks is converted to charge 1; in dense spectra (e.g. timsTOF) more of the two-peak "
                       "multiply charged envelopes are chance pairs. 3 is the earlier behaviour.",
                       {"advanced"});
    defaults_.setMinInt("fragment:deisotope_min_peaks", 2);
    defaults_.setMaxInt("fragment:deisotope_min_peaks", 10);
    defaults_.setValue("fragment:deisotope_charge_cap", "precursor",
                       "Fragment charges that MS2 deisotoping tries. 'precursor' tries charges up to the precursor charge, but at most 3 "
                       "(1-3 when the precursor charge is unknown): a fragment cannot carry more charges than its precursor, so a chance "
                       "peak 1/3 Th above an ion of a 2+ precursor cannot turn the ion into a '3+' envelope at a wrong m/z. This is "
                       "Sage-like; Sage caps at the precursor charge without the limit of 3 and counts an unknown charge as 3. 'none' "
                       "tries charges 1-3 in every spectrum (the earlier behaviour).",
                       {"advanced"});
    defaults_.setValidStrings("fragment:deisotope_charge_cap", {"precursor", "none"});
    defaults_.setValue("fragment:deisotope_sum_intensity", "true",
                       "Give the monoisotopic peak of each isotope envelope the summed intensity of the envelope (Sage-like). "
                       "'false' keeps the monoisotopic peak's own intensity (the earlier behaviour).",
                       {"advanced"});
    defaults_.setValidStrings("fragment:deisotope_sum_intensity", {"true", "false"});


    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("fragment:min_ion_index", 2, "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.");
    defaults_.setMinInt("fragment:min_ion_index", 0);

    defaults_.setValue("peaks:keep_n", 0, "Maximum number of MS2 peaks kept per spectrum (NLargest) before scoring. 0 = auto: a resolution-aware cap (~400 for high-res, ~80 at 0.5 Da, floor 60) that avoids the spurious low-resolution fragment matching which otherwise inflates HyperScore for targets and decoys alike. Set an explicit value to override (e.g. 400 for the legacy behavior).", {"advanced"});
    defaults_.setMinInt("peaks:keep_n", 0);
    defaults_.setValue("peaks:window_top", 20, "Maximum number of MS2 peaks kept per 100 Da window (WindowMower) before scoring.", {"advanced"});
    defaults_.setMinInt("peaks:window_top", 1);
    defaults_.setValue(
      "fragment:query_spectrum", "auto",
      "Spectrum used for candidate retrieval. 'raw' queries the fragment index with all positive-intensity peaks, before "
      "deisotoping and local/top-N filtering, so that spectra whose matching fragments are not among the retained peaks still get "
      "candidates; it costs an additional peak list per spectrum and more fragment matches. 'processed' queries with the scoring peak "
      "list. 'auto' uses 'raw' for high-resolution fragments (<= 0.1 Da / <= 100 ppm) and 'processed' otherwise: at low resolution, "
      "the unfiltered peaks match many random fragments and crowd out the correct candidates. Scoring always uses the processed spectrum.",
      {"advanced"});
    defaults_.setValidStrings("fragment:query_spectrum", {"auto", "raw", "processed"});
    defaults_.setValue("peaks:window_type", "auto",
                       "Local peak filtering. 'auto' uses jump_full for high-resolution fragments (<= 0.1 Da / <= 100 ppm), "
                       "preserving the strongest peaks even in a short final 100 Da window, and legacy jump filtering otherwise. "
                       "'jump' scales the last window's peak quota by its observed width; 'jump_full' keeps the full quota in every window.",
                       {"advanced"});
    defaults_.setValidStrings("peaks:window_type", {"auto", "jump", "jump_full"});

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<std::string> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'. A terminus and a residue carry one modification each: a variable terminal modification (e.g. 'Acetyl (Protein N-term)', 'Gln->pyro-Glu (N-term Q)') is not searched where a fixed one sits on that terminus (e.g. 'TMT6plex (N-term)'), and a variable residue modification (e.g. 'Glutathione (C)') is not searched where a fixed one sits on that residue (e.g. 'Carbamidomethyl (C)').");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    vector<std::string> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));

    defaults_.setValue("decoys", "auto",
      "Decoy handling for target-decoy FDR. "
      "'auto': ensure decoys are available — reuse decoys already present in the database "
      "(marker auto-detected, prefix or suffix) or generate them if none are found. Supplied decoys are reused only "
      "at 0.5 to 2 decoys per target protein: target-decoy FDR needs about one decoy per target, so a database with "
      "far fewer or more is refused (use 'generate'). "
      "'generate': always (re)generate decoys from the target proteins (any pre-existing "
      "decoys are removed first to avoid decoy-of-decoy entries). "
      "'ignore': search the target proteins only — any decoys present in the database are removed.");
    defaults_.setValidStrings("decoys", {"auto","generate","ignore"} );

    defaults_.setValue("decoy_prefix", "DECOY_", "Accession prefix prepended to decoy proteins when ProSE generates them (modes 'auto' without existing decoys, and 'generate'). Pre-existing decoys in the database are detected automatically (any common marker, used as prefix or suffix) and do not require this setting.", {"advanced"});
    defaults_.setValue("decoy_method", "reverse",
      "How ProSE generates decoys (decoys 'generate', or 'auto' with a database without decoys). 'reverse' (ProSEAlgorithm's "
      "default; the ProSE tool's default is 'keep_termini'): "
      "every enzymatic segment of the protein is reversed except its last residue (the last segment entirely), an "
      "initial Met stays first with peptide:clip_nterm_methionine. 'keep_termini': the construction of OpenMS "
      "DecoyDatabase (method reverse, Decoy:keepPeptideNTerm/keepPeptideCTerm true) with the search enzyme and "
      "Decoy:non_shuffle_pattern 'KRP': in every segment of a digest without missed cleavages, the first and the last "
      "residue and every K, R and P stay in place and the rest is reversed. For enzymes whose cleavage sites are fixed by "
      "these residues (Trypsin, Trypsin/P, Lys-C, Arg-C and the like) every target peptide then has a decoy with the same "
      "start, length, mass and cleavage sites: also every protein N-terminal and Met-clipped one, so that "
      "peptide:clip_nterm_methionine and protein N-terminal modifications add the same number of target and decoy "
      "candidates (except peptides that equal their decoy, e.g. made of short segments or of K, R and P only; and "
      "modifications specific to a peptide's first residue, e.g. pyro-Glu on an N-terminal Q, apply to the target and "
      "its decoy alike only where that residue is kept in place). Semi- or "
      "non-specific digestion: as 'reverse'. Other enzymes: a warning names the decoys whose cleavage sites differ from the "
      "target's. No effect on decoys taken from the database (decoys 'auto' with decoys in it). The results differ from "
      "'reverse'.", {"advanced"});
    defaults_.setValidStrings("decoy_method", {"reverse", "keep_termini"});
    defaults_.setValue("decoy_parity_report", "false",
      "Diagnostic: log the candidate counts of the fragment index per label (target, decoy) and subclass (internal, "
      "protein N-terminal with and without initial M, Met-clipped; and their protein N-terminally acetylated occurrences, "
      "none with a fixed N-terminal modification such as TMT, which leaves no room for a variable one), and a "
      "checksum of generated decoys. Output unchanged.",
      {"advanced"});
    defaults_.setValidStrings("decoy_parity_report", {"true", "false"});

    defaults_.setValue("annotate:PSM",  std::vector<std::string>{"ALL"}, "Annotations added to each PSM.");
    defaults_.setValidStrings("annotate:PSM",
      std::vector<std::string>{
        "ALL",
        Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
        Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
        Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
        Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION,
        Constants::UserParam::NUM_MATCHED_PEAKS,
        Constants::UserParam::MATCHED_PREFIX_IONS,
        Constants::UserParam::MATCHED_SUFFIX_IONS,
        Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE,
        Constants::UserParam::MATCHED_ION_CURRENT,
        Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM,
        Constants::UserParam::HYPERSCORE_ZSCORE,
        Constants::UserParam::LN_NUM_CANDIDATES,
        Constants::UserParam::MATCHED_ION_CURRENT_FRACTION,
        Constants::UserParam::COMPLEMENTARY_IONS_FRACTION}
      );

    defaults_.setValue("annotate:self_trained_ion_priors", "true",
      "Learn per-run fragment ion likelihoods from the confident PSMs of each file and add the Percolator features of "
      "annotate:ion_prior_features to every hit. The spectra are split into two halves "
      "by scan parity. Each half trains a model on its rank-one target hits at target-decoy competition q <= "
      "annotate:ion_prior_train_fdr of the native score (their reversed sequences on the same spectra are the noise "
      "model), and every hit is scored by the model of the other half, so no PSM is scored by a model that saw its "
      "spectrum or its label. Needs annotate:ion_prior_min_psms training PSMs in each half, otherwise the features are 0. "
      "Not applied to target-only searches (decoys=ignore), which get no features. Native scores and candidate selection "
      "are unchanged. Independent of annotate:PSM.", {"advanced"});
    defaults_.setValidStrings("annotate:self_trained_ion_priors", {"true", "false"});
    defaults_.setValue("annotate:ion_prior_model", "rich",
      "Contexts and outcomes of the ion priors. 'rich': ion series, precursor charge, fragment charge, relative position, "
      "the residues at the cleavage site (bond before P, bond after D or E) and whether the complementary singly charged "
      "ion matched; outcome: intensity rank x mass error (|error| / tolerance below 0.25, below 0.5, up to 1). 'basic': "
      "ion series, precursor charge, fragment charge and relative position; outcome: intensity rank.", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_model", {"rich", "basic"});
    defaults_.setValue("annotate:ion_prior_peaks", "all",
      "Peaks the ion priors rank and match. 'all': every peak after deisotoping, before the peaks:window_top and "
      "peaks:keep_n filters, so weak fragment ions count as evidence for the model without entering HyperScore. "
      "'scored': the peaks that are scored.", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_peaks", {"all", "scored"});
    defaults_.setValue("annotate:ion_prior_max_fragment_charge", 2,
      "Highest fragment charge of the ions of the ion priors; capped at the precursor charge - 1 (at least 1).", {"advanced"});
    defaults_.setMinInt("annotate:ion_prior_max_fragment_charge", 1);
    defaults_.setValue("annotate:ion_prior_train_fdr", 0.01,
      "Target-decoy competition q-value threshold selecting the training PSMs of annotate:self_trained_ion_priors.", {"advanced"});
    defaults_.setMinFloat("annotate:ion_prior_train_fdr", 0.0);
    defaults_.setMaxFloat("annotate:ion_prior_train_fdr", 1.0);
    defaults_.setValue("annotate:ion_prior_min_psms", 100,
      "Minimum number of training PSMs of annotate:self_trained_ion_priors in each half of the spectra; with fewer, the "
      "features are 0.", {"advanced"});
    defaults_.setMinInt("annotate:ion_prior_min_psms", 1);
    defaults_.setValue("annotate:ion_prior_features",
      std::vector<std::string>{Constants::UserParam::ION_PRIOR_LLR, Constants::UserParam::ION_PRIOR_EXPLAINED},
      "Features of annotate:self_trained_ion_priors written to every hit and listed for Percolator: "
      "ion_prior_llr (summed log-likelihood ratio of the ion outcomes, signal vs reversed noise), ion_prior_explained "
      "(share of the predicted ion presence that was observed), ion_prior_topk_observed (share of the 6 ions most likely "
      "present that were observed).", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_features",
      std::vector<std::string>{Constants::UserParam::ION_PRIOR_LLR, Constants::UserParam::ION_PRIOR_EXPLAINED,
                               Constants::UserParam::ION_PRIOR_TOPK_OBSERVED});

    defaults_.setValue("annotate:local_fragment_evidence", "true",
      "Add chance_match_surprise and mass_competition_evidence annotations and Percolator features. "
      "Uses local peak density and alternative intact/neutral-loss fragment assignments. "
      "Retains an additional peak list before window/top-N filtering, increasing memory use. "
      "Does not change native scores or candidate selection. Independent of annotate:PSM.", {"advanced"});
    defaults_.setValidStrings("annotate:local_fragment_evidence", {"true", "false"});
    defaults_.setSectionDescription("annotate", "Annotation Options");

    defaults_.setValue("peptdeep:enable", "false", "Add PeptDeep prediction-based rescoring features (ms2_cosine, ms2_spectral_angle, ms2_pearson, ms2_frac_pred_found, rt_abs_error) to every PSM. Uses the models shipped in OpenMS' 'share/OpenMS/models' unless 'peptdeep:ms2_model'/'peptdeep:rt_model' name others. Requires an OpenMS built with ONNX support.");
    defaults_.setValidStrings("peptdeep:enable", {"true", "false"});
    defaults_.setValue("peptdeep:ms2_model", "models/peptdeep_ms2_dynamic.onnx", "PeptDeep MS2 fragment-intensity ONNX model, used when 'peptdeep:enable' is true. A relative name is resolved against OpenMS' shared-data directory; an absolute path is used as given.");
    defaults_.setValue("peptdeep:rt_model", "models/peptdeep_rt_dynamic.onnx", "PeptDeep retention-time ONNX model. See 'peptdeep:ms2_model'.");
    defaults_.setValue("peptdeep:instrument", "QE", "Instrument class passed to the PeptDeep MS2 model. An instrument that is not in this list is best approximated by the closest one that is.");
    // peptdeep's own names, in peptdeep's own order, and registered whether or not this build has
    // ONNX so that an ini file means the same thing either way.
    std::vector<std::string> peptdeep_instruments;
    peptdeep_instruments.reserve(ML::ALPHAPEPTDEEP_INSTRUMENTS.size());
    for (const std::string_view name : ML::ALPHAPEPTDEEP_INSTRUMENTS) { peptdeep_instruments.emplace_back(name); }
    defaults_.setValidStrings("peptdeep:instrument", peptdeep_instruments);
    defaults_.setValue("peptdeep:nce", -1.0, "Normalised collision energy for the PeptDeep MS2 model. Negative selects it automatically from the collision energy recorded in the spectra, refined by scoring a small grid on confident PSMs.");
    defaults_.setValue("peptdeep:rt_model_type", "b_spline", "Model mapping predicted onto observed retention time.", {"advanced"});
    defaults_.setValidStrings("peptdeep:rt_model_type", {"b_spline", "lowess", "linear"});
    defaults_.setValue("peptdeep:batch_size", 64, "Peptides per ONNX inference call. The calls run side by side, one per thread: small ones keep the threads evenly busy and the memory low (on iPRG2015 with 16 threads, 64 took 40% less time than 500 and 1 GB less memory).", {"advanced"});
    defaults_.setMinInt("peptdeep:batch_size", 1);
    defaults_.setSectionDescription("peptdeep", "PeptDeep (ONNX) prediction-based rescoring features");

    defaults_.setValue("peptide:min_size", 7, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");
    defaults_.setValue("peptide:clip_nterm_methionine", "true",
                       "Also search protein N-terminal peptides after removal of their initial methionine. "
                       "The retained-M form is still searched; length, mass and missed-cleavage limits apply to each form. "
                       "Set false for the previous search space.");
    defaults_.setValidStrings("peptide:clip_nterm_methionine", {"true", "false"});
    defaults_.setValue("peptide:deduplicate", "true",
                       "Index each exact peptidoform once before candidate selection. Protein mappings still list every protein occurrence. "
                       "Across database chunks, count each peptide/charge/isotope hypothesis once; retaining queried keys adds memory per spectrum. "
                       "With low_memory, the repeated peptidoforms among the candidates of a spectrum are dropped instead. "
                       "Set false for occurrence-based legacy candidates.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:deduplicate", {"true", "false"});
    defaults_.setValue("peptide:protein_mapping", "index",
                       "How the hits are mapped to their proteins (peptide evidences, target_decoy, protein_references, protein hits). "
                       "'index': from the digest of the fragment index, which holds every protein occurrence of a candidate (with a "
                       "protein-terminal modification: from a pass over the database for the hits' sequences); the result is "
                       "that of PeptideIndexing, which runs instead wherever the index cannot reproduce it exactly (e.g. low_memory, a specificity "
                       "other than full, a chunked database, symbols other than letters or long stretches "
                       "of ambiguous residues in the database). 'PeptideIndexing': always search every hit in the whole database "
                       "(Aho-Corasick), as before.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:protein_mapping", {"index", "PeptideIndexing"});
    defaults_.setValue("peptide:enzyme_specificity", "full",
      "Enzyme cleavage specificity required for both peptide termini.\n"
      "  'full' : both termini must be enzyme-specific (canonical, e.g. tryptic).\n"
      "  'semi' : only one terminus needs to be enzyme-specific (semi-tryptic).\n"
      "  'none' : no enzyme constraint at either terminus; every substring of length\n"
      "           [min_size, max_size] is enumerated. Use this for immunopeptidomics\n"
      "           (e.g. HLA peptides 8..12mers). For very large search spaces consider\n"
      "           tightening 'peptide:min_size'/'peptide:max_size'.");
    defaults_.setValidStrings("peptide:enzyme_specificity", {"full", "semi", "none"});
    defaults_.setValue("peptide:motif", "", ""); // description: FragmentIndex::restrictSearchParameters()
    defaults_.setSectionDescription("peptide", "Peptide Options");

    // forwarded to FragmentIndex
    defaults_.setValue("low_memory", "false",
      "Non-specific searches (peptide:enzyme_specificity=none; ignored otherwise): index the peptide masses, as the "
      "prefixes of one mother peptide per protein position, instead of the fragments of every peptidoform, and "
      "generate the fragments of the peptidoforms in the precursor window at query time. Finds the same candidates "
      "with a fraction of the memory (e.g. for immunopeptidomics against a human database), but takes more time per "
      "spectrum, the more the wider the precursor window.");
    defaults_.setValidStrings("low_memory", {"true", "false"});

    defaults_.setValue("report:top_hits", 1,
                       "Maximum number of top scoring hits per spectrum that are reported. With 'rescore', only the best "
                       "hit per spectrum is reported after the rescoring, and the hits that Percolator rescores are set by "
                       "'rescoring:top_hits' (report:top_hits then applies to a run that is not rescored and to -out_pin).");
    defaults_.setValue("report:isotope_error_convention", "observed_minus_theoretical",
                       "Sign of the precursor isotope error reported for each PSM (meta value 'isotope_error', in 13C "
                       "spacings of 1.00336 Da). 'observed_minus_theoretical': +1 means the observed precursor lies one "
                       "13C spacing above the peptide (its first 13C isotope peak was selected), as in MS-GF+, pepXML and "
                       "Sage; PercolatorInfile and PercolatorAdapter (-out_pin, -percolator_executable) remove the offset "
                       "from the precursor mass difference with this sign. The search parameters of the result record it "
                       "as 'isotope_error_convention'. 'theoretical_minus_observed': the opposite sign, as written by "
                       "earlier ProSE versions (no record in the search parameters); with it the Percolator input "
                       "doubles the offset in dm/absdm instead of removing it.",
                       {"advanced"});
    defaults_.setValidStrings("report:isotope_error_convention", {"observed_minus_theoretical", "theoretical_minus_observed"});
    defaults_.setSectionDescription("report", "Reporting Options");

    defaults_.setValue("rescore", "true",
                       "Rescore the PSMs of every search run (input file) in-process with Percolator (OpenMS' built-in "
                       "implementation, as PercolatorAdapter uses by default) before the FDR filters: Percolator q-values become "
                       "the main score (the SVM score and PEP are kept as meta values), Percolator rescores the top "
                       "'rescoring:top_hits' hits of each spectrum and only its best one per spectrum is kept, "
                       "and FDR:PSM and FDR:protein are applied to the rescored PSMs. Needs decoys (see 'decoys') and at least "
                       "100 PSMs: a run without them, or for which Percolator fails, is not rescored and keeps its search engine "
                       "scores, which the FDR filters then use (the ProSE tool warns). If that happens to some files of a "
                       "multi-file search, the files rescored on their own keep their results, and the others get the "
                       "scores of a model trained on the PSMs of those (on the scale of its training PSMs, with their PEP "
                       "at the same score; q-values by target-decoy competition of each file's PSMs on that score); if "
                       "Percolator can rescore no file on its own, all files are rescored together with one model. Only "
                       "if that fails too, the files keep their search engine scores, and the pooled evaluation of the "
                       "files (the ProSE tool's -out_merged) scores the proteins by ProSE's score 'ln(hyperscore)' of all "
                       "PSMs instead of Percolator's PEP. See percolator:max_train_psms for the training PSMs. Set to "
                       "'false' to report the search engine scores.");
    defaults_.setValidStrings("rescore", {"true", "false"});
    defaults_.setValue("percolator:max_train_psms", 500000,
                       "With 'rescore': the most PSMs Percolator trains a model on (0: all). Every PSM is scored, also "
                       "beyond the cap. For a run with more PSMs, Percolator samples the training PSMs by spectrum (see "
                       "'rescoring:subset_mode'); for the model of the files of a multi-file search too small for a model of "
                       "their own (trained on the files rescored on their own), each of these files contributes an equal "
                       "share, so that none dominates. The run log says when the training PSMs were sampled. Searches below "
                       "the cap are not affected.", {"advanced"});
    defaults_.setMinInt("percolator:max_train_psms", 0);
    defaults_.setSectionDescription("percolator", "Percolator rescoring ('rescore')");
    defaults_.setValue("rescore_solver", "newton",
                       "Solver of Percolator's SVMs with 'rescore' (Percolator's parameter 'solver'): 'newton' (default; "
                       "explicit-Hessian Newton method with exact line search and warm starts: the same SVM optimum for a "
                       "given training set as TRON, within the solver tolerance, much faster on large runs; as the earlier "
                       "iterations of the cross-validation stop at a looser tolerance, they may choose slightly different "
                       "training sets, so the models can differ from TRON's), 'tron' (Percolator's trust region Newton "
                       "method, the solver of the percolator executable and of earlier ProSE versions) or "
                       "'lda' (Fisher's linear discriminant instead of the SVM: fastest, an approximation).", {"advanced"});
    defaults_.setValidStrings("rescore_solver", {"newton", "tron", "lda"});
    defaults_.setValue("rescoring:top_hits", 2,
                       "With 'rescore': the top scoring hits per spectrum (by the search engine score) that Percolator "
                       "rescores: its training and scoring rows, of which its target-decoy competition keeps the best per "
                       "spectrum, which is reported. 2 (default) gives Percolator the runner-up of each spectrum as well; "
                       "more hits cost time and memory. 0 = 'report:top_hits' (as earlier "
                       "ProSE versions, whose Percolator rescored the reported hits). The search keeps these hits for "
                       "Percolator (with -out_pin also the reported ones): a run that Percolator does not rescore (fewer "
                       "than 100 PSMs, no decoy PSM, a failure) reports at most rescoring:top_hits (if above 0) of its "
                       "report:top_hits hits per spectrum; with decoys 'ignore' (no rescoring) report:top_hits.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:top_hits", 0);
    defaults_.setValue("rescoring:subset_mode", "out_of_fold",
                       "With a training cap ('percolator:max_train_psms' > 0) on a run above it (Percolator's parameter 'subset_mode'): 'out_of_fold' "
                       "(default) assigns every spectrum to one of Percolator's three cross-validation folds first, draws "
                       "the training sample (whole spectra, at random from the entire run, seeded), trains each fold's "
                       "model on the sampled spectra of the other folds only, and scores every PSM with the model of its "
                       "own fold, which never saw its spectrum. 'out_of_fold_avg': the same, but the PSMs of the spectra "
                       "outside the sample (which no model saw) are scored by the average of the fold models. "
                       "'in_sample' = the percolator executable's --subset-max-train: cross-validation on the sample and "
                       "all PSMs scored by the average of the fold models (the sampled PSMs by two models that trained on "
                       "them). 'learning_curve' (also without a cap): out-of-fold training on nested samples of "
                       "'rescoring:curve_start', twice, four times ... as many PSMs, up to 'percolator:max_train_psms' (0: "
                       "1,000,000) or the whole run but the validation spectra (the end of the random order, at most "
                       "250,000 PSMs or 20 %); each fold's model stops at the first doubling that adds less than "
                       "'rescoring:curve_min_gain' or 2 standard errors of targets at 1 % among the validation spectra of the "
                       "other folds (no label of its own fold is used), and is not trained further; every PSM is scored by its fold's model. The sizes and the curve are "
                       "recorded as Percolator:subset_curve.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:subset_mode", {"out_of_fold", "out_of_fold_avg", "in_sample", "learning_curve"});
    defaults_.setValue("rescoring:curve_start", 100000, "'rescoring:subset_mode' learning_curve: PSMs of the first sample.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:curve_start", 1);
    defaults_.setValue("rescoring:curve_min_gain", 0.002,
                       "'rescoring:subset_mode' learning_curve: a fold stops at the first doubling that adds less than "
                       "this fraction of validation targets.", {"advanced"});
    defaults_.setMinFloat("rescoring:curve_min_gain", 0.0);
    defaults_.setValue("rescoring:iterations", 10,
                       "Iterations of Percolator's cross-validation training (the percolator executable's --maxiter). "
                       "Fewer are faster; the results differ.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:iterations", 1);
    defaults_.setValue("rescoring:tron_epsilon", 1e-7,
                       "Stopping tolerance of the SVM solver (TRON) of Percolator's training (Percolator: 1e-7). Larger "
                       "values (e.g. 1e-3) train much faster on large runs; the results differ.",
                       {"advanced"});
    defaults_.setMinFloat("rescoring:tron_epsilon", 1e-15);
    defaults_.setValue("rescoring:warm_start", "false",
                       "Each SVM training of an iteration starts from the weights of the previous iteration instead of "
                       "0 (faster; the results differ).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:warm_start", {"true", "false"});
    defaults_.setValue("rescoring:single_cneg", "false",
                       "Percolator trains with one Cneg (3 x the target/decoy ratio x Cpos) instead of selecting among "
                       "three by cross-validation (faster; the results differ).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:single_cneg", {"true", "false"});
    defaults_.setValue("rescoring:use_pi0", "true",
                       "Percolator's q-values and PEPs with its pi0 estimate (default) or with pi0 = 1, as the percolator "
                       "executable computes them with target-decoy competition (-Y). 'false' gives slightly fewer PSMs and "
                       "a lower measured false discovery proportion; the results differ. 'auto' decides as the executable (OpenMS "
                       "develop's in-process default since #10721): with ProSE's target-decoy competition that is 'false'.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:use_pi0", {"true", "false", "auto"});
    defaults_.setValue("rescoring:calc_mass_feature", "true",
                       "Percolator's training uses the CalcMass column (the theoretical m/z) as a feature. 'true' (ProSEAlgorithm's "
                       "default, as ProSE int15); 'false' (the ProSE tool's default): as the percolator executable, which reads it as metadata (OpenMS develop since "
                       "#10717). The results differ.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:calc_mass_feature", {"true", "false"});
    defaults_.setValue("rescoring:fixed_enzyme_termini", "false",
                       "A protein terminus counts as an enzymatic site for the enzN/enzC features of the rescoring and of "
                       "'-out_pin' (Percolator's convention). Default 'false' keeps the former features, in which a peptide "
                       "at a protein N-terminus has enzN = 0 and one at a protein C-terminus that does not end in a "
                       "cleavage residue has enzC = 0; the results differ. The peptide after an initial Met (the Met-clipped protein "
                       "N-terminal peptide of peptide:clip_nterm_methionine) counts as N-terminal too. Experimental: tested "
                       "only with decoy_method 'keep_termini'; without the Met-clip rule, with the decoys of a database of "
                       "reversed proteins (which rarely start with M) enzN = 0 was left almost only to Met-clipped peptides "
                       "of targets and the FDR was no longer controlled (entrapment FDP 13-24 % at 1 % q).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:fixed_enzyme_termini", {"true", "false"});
    defaults_.setValue("rescoring:terminal_features", "false",
                       "Two more rescoring features: 'protein_nterm' (1 if the peptide occurs at a protein N-terminus) and "
                       "'met_clipped' (1 if it occurs only as the Met-clipped protein N-terminal peptide of "
                       "peptide:clip_nterm_methionine). They let the rescoring learn the different rate of correct "
                       "matches of these small peptide classes. Use them only with decoys that populate both classes like "
                       "the targets (decoy_method 'keep_termini'): with decoys that rarely start with M (reversed "
                       "proteins) they tell targets from decoys. The results differ.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:terminal_features", {"true", "false"});
    defaults_.setValue("rescoring:stratum_q_values", "false",
                       "Also q-values per stratum (meta values 'stratum' and 'stratum_q_value'): the rescored PSMs of "
                       "protein N-terminal peptides (incl. Met-clipped ones) and the others each get target-decoy "
                       "competition of their own ((decoys + 1) / targets, as mokapot's grouped confidence). For a list of "
                       "N-terminal peptides at an FDR, use the stratum q-value; the main score stays the pooled q-value. "
                       "Needs decoys that populate the N-terminal classes like the targets (decoy_method 'keep_termini'). With "
                       "rescore_method 'reset' the values are computed on RESET's scores (both decoy sets count as decoys "
                       "with reset:decoys 'two'), not as RESET's estimate.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:stratum_q_values", {"true", "false"});
    defaults_.setValue("rescoring:peptide_level", "true",
                       "Also Percolator's peptide-level q-values (default; 'false' leaves them out), as the percolator executable reports peptides (without -U): "
                       "after the target-decoy competition of the spectra, the best PSM of each peptide (modified sequence) "
                       "and label competes, and q-values and PEPs are estimated on these peptides. Every rescored hit gets "
                       "the values of its peptide as meta values 'percolator_peptide_q_value' and 'percolator_peptide_pep'; "
                       "a list of peptides at an FDR is the target peptides with percolator_peptide_q_value at most that FDR "
                       "(the distinct peptides of the PSMs at an FDR are not FDR-controlled at the peptide level). The main "
                       "score stays the PSM-level q-value. Default Percolator rescoring only (not with rescore_method 'reset', "
                       "which reports RESET:peptide_q_value, nor with rescore_reuse models).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:peptide_level", {"true", "false"});
    defaults_.setSectionDescription("rescoring", "Options of the Percolator training of the rescoring (see 'rescore')");
    defaults_.setValue("rescore_method", "percolator",
                       "With 'rescore': 'percolator' = Percolator's PSM-level FDR (default); 'reset' = Percolator-RESET "
                       "(Freestone et al., J Proteome Res 2025): half of the decoy peptides (stems = unmodified sequences, by a "
                       "random coin per stem) train the model, the other half estimates the FDR, which is controlled at the "
                       "peptide (stem) level: every rescored hit gets the q-value of its stem as meta value "
                       "'RESET:peptide_q_value'. The main score (q-value) stays PSM-level, but under RESET it is a heuristic "
                       "(not FDR-controlled): 2 (estimating decoy PSMs + 1) / target PSMs. Experimental.", {"advanced"});
    defaults_.setValidStrings("rescore_method", {"percolator", "reset"});
    defaults_.setValue("reset:salt", "",
                       "rescore_method 'reset': salt of the decoy coins (unsigned 64-bit integer); empty = random per run. "
                       "The salt is recorded in the search parameters ('RESET:salt'). Never pick the best of several salts.",
                       {"advanced"});
    defaults_.setValue("reset:learner", "svm",
                       "rescore_method 'reset': 'svm' = linear SVM; 'mlp' = small MLP ensemble (experimental).", {"advanced"});
    defaults_.setValidStrings("reset:learner", {"svm", "mlp"});
    defaults_.setValue("reset:mlp_members", 8, "rescore_method 'reset', learner 'mlp': networks per ensemble.", {"advanced"});
    defaults_.setMinInt("reset:mlp_members", 1);
    defaults_.setValue("reset:mlp_hidden", 16, "rescore_method 'reset', learner 'mlp': units of each of the two hidden layers.", {"advanced"});
    defaults_.setMinInt("reset:mlp_hidden", 1);
    defaults_.setValue("reset:mlp_epochs", 30, "rescore_method 'reset', learner 'mlp': training epochs of the first round.", {"advanced"});
    defaults_.setMinInt("reset:mlp_epochs", 1);
    defaults_.setValue("reset:mlp_warm_epochs", 3,
                       "rescore_method 'reset', learner 'mlp': > 0 = the networks continue from the previous round for this many "
                       "epochs (warm starts); 0 = trained afresh (mlp_epochs) every round.", {"advanced"});
    defaults_.setMinInt("reset:mlp_warm_epochs", 0);
    defaults_.setValue("reset:mlp_crossfit", "true",
                       "rescore_method 'reset', learner 'mlp': score each peptide by the networks that did not train on it (3 folds "
                       "by sequence hash); 'false' = one final ensemble (memorizes its training decoys: fewer peptides).", {"advanced"});
    defaults_.setValidStrings("reset:mlp_crossfit", {"true", "false"});
    defaults_.setValue("reset:mlp_norm_min_decoys", 1000,
                       "rescore_method 'reset', learner 'mlp', cross-fitted: normalize the folds' scores (as Percolator) only if every "
                       "fold has at least this many training decoys (0 = always).", {"advanced"});
    defaults_.setMinInt("reset:mlp_norm_min_decoys", 0);
    defaults_.setValue("reset:representative", "model",
                       "rescore_method 'reset': the PSM that represents a peptide (stem): 'model' (default) = after each training "
                       "round, best score under the new model (+2 % peptides on DIA; the FDR then also assumes that null targets "
                       "and decoys are exchangeable over their whole PSM collections, including their sizes); 'search' = best "
                       "search engine score (the paper).", {"advanced"});
    defaults_.setValidStrings("reset:representative", {"search", "model"});
    defaults_.setValue("reset:iterations", 25,
                       "rescore_method 'reset': training rounds (the paper: 5; on DIA pseudo-spectra 5 rounds stop far from "
                       "convergence: 82 % of the peptides of 20 rounds on HeLa).", {"advanced"});
    defaults_.setMinInt("reset:iterations", 1);
    defaults_.setValue("reset:audit", "false",
                       "rescore_method 'reset': feature-symmetry audit before training (log).", {"advanced"});
    defaults_.setValidStrings("reset:audit", {"true", "false"});
    defaults_.setValue("reset:decoys", "single",
                       "rescore_method 'reset': 'single' = one decoy database (default); 'two' = the two-decoy mode of "
                       "Percolator-RESET: the database holds a second decoy set of the same construction (accessions starting "
                       "with reset:second_decoy_prefix, also recognized as decoys by the decoy marker), the spectra are "
                       "searched against targets and both decoy sets, and RESET estimates the FDR with the factors of two decoy "
                       "databases (fewer peptides lost to the decoy split, less variation between salts; the search is larger).",
                       {"advanced"});
    defaults_.setValidStrings("reset:decoys", {"single", "two"});
    defaults_.setValue("reset:second_decoy_prefix", "DECOY_X_",
                       "rescore_method 'reset', reset:decoys 'two': accession prefix of the proteins of the second decoy set.",
                       {"advanced"});
    defaults_.setSectionDescription("reset", "Percolator-RESET options (rescore_method 'reset')");
    defaults_.setValue("reset:stem_table", "",
                       "rescore_method 'reset': write the stems with their q-values (tab-separated) to this file, or, if it "
                       "ends with '/' or is a directory, to <dir>/<spectra file base name>.reset_stems.tsv per run.",
                       {"advanced"});
    defaults_.setValue("rescore_reuse:mode", "per_file",
                       "Which PSMs train the Percolator model of a run ('rescore'). 'per_file': the run itself, with Percolator's "
                       "three-fold cross-validation. The other modes cross-fit: the spectra are split into two halves at random "
                       "(by Percolator's seed), one model is trained per half on about 'rescore_reuse:max_train' PSMs, and every PSM is "
                       "scored by the model of the other half, so that no PSM is scored by a model that was trained on it (5 % of the "
                       "spectra of each half are not trained on: their decoys put the scores of their own half's model on one scale; they are scored by the other model). 'joint' (a "
                       "multi-file search): the training PSMs are drawn from all files, the same number from each, and the files are "
                       "rescored once all are searched (mokapot's joint model). 'static' (a multi-file search): they are drawn from the first "
                       "'rescore_reuse:static_files' files, and every later file is rescored with these models after its search, unless "
                       "the features of its PSMs differ from those of the training PSMs by more than 'rescore_reuse:max_shift': "
                       "then it is rescored as with 'per_file'. 'own': every file trains its two models on its own PSMs. In every mode the "
                       "q-values and PEPs of a file come from its own target-decoy competition (Percolator's, as with 'per_file'). 'joint' and "
                       "'static' need an unchunked closed search of several files and keep all files' PSMs in memory until they are "
                       "rescored; otherwise 'per_file' applies.",
                       {"advanced"});
    defaults_.setValidStrings("rescore_reuse:mode", {"per_file", "joint", "static", "own"});
    defaults_.setValue("rescore_reuse:max_train", 250000,
                       "Training PSMs per model of 'rescore_reuse:mode' (0: all PSMs of the half but the calibration spectra). The spectra are "
                       "drawn at random, all PSMs of a drawn spectrum together, so the number is met on average; the files share it "
                       "equally (a file with fewer PSMs gives all).",
                       {"advanced"});
    defaults_.setMinInt("rescore_reuse:max_train", 0);
    defaults_.setValue("rescore_reuse:static_files", 1, "With 'rescore_reuse:mode' static: the number of (first) input files that train the models.",
                       {"advanced"});
    defaults_.setMinInt("rescore_reuse:static_files", 1);
    defaults_.setValue("rescore_reuse:max_shift", 0.5,
                       "With 'rescore_reuse:mode' static: the largest domain shift for which a file is rescored with the static models. "
                       "The shift is max |2 AUC - 1| over the features the models use, AUC comparing the file's PSMs with the "
                       "training PSMs, targets and decoys alike (0: same distribution, 1: separated).",
                       {"advanced"});
    defaults_.setMinFloat("rescore_reuse:max_shift", 0.0);
    defaults_.setMaxFloat("rescore_reuse:max_shift", 1.0);
    defaults_.setSectionDescription("rescore_reuse", "Percolator model reuse across the files of a search ('rescore')");

    defaults_.setValue("FDR:PSM", 0.0,
                       "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR; 0 disables the filter). Target and decoy PSMs are filtered "
                       "alike by the threshold; decoys that pass are kept (no decoy-specific stripping here — decoys are removed only at "
                       "protein-FDR finalization). With 'rescore' (the default), the threshold applies to the Percolator q-values, the main "
                       "score of the rescored PSMs. Without rescoring ('rescore' false, or a run that Percolator did not rescore), 0 computes "
                       "no q-values: every PSM keeps its search engine score ('ln(hyperscore)') as main score, e.g. as input for Percolator "
                       "or a downstream FDR step. Above 0, the target-decoy q-value (see FDR:PSM_groups) replaces the main score (the search "
                       "engine score is kept as meta value 'ln(hyperscore)_score') and only the best hit of each spectrum is kept. If no PSM "
                       "matched a decoy, the q-value of the T PSMs is the conservative 1/T (target-decoy competition with +1), "
                       "so a lower threshold reports no PSM instead of lifting the filter. Needs decoys ('decoys' auto or "
                       "generate); with 'ignore', a threshold above 0 is skipped with a warning (the QPX run view of the ProSE tool's "
                       "-qpx_extended records per input whether the filter was applied: 'psm_fdr_applied'). Afterwards, the "
                       "proteins that no remaining PSM references are removed.");
    defaults_.setMinFloat("FDR:PSM", 0.0);
    defaults_.setMaxFloat("FDR:PSM", 1.0);
    defaults_.setValue("FDR:PSM_groups", "scored_charges",
                       "Target-decoy competitions of FDR:PSM. 'pooled' ranks all PSMs together. 'scored_charges' computes q-values "
                       "separately for PSMs scored with different numbers of fragment charges (scoring:fragment_charges): HyperScore "
                       "grows with the number of theoretical ions, so when multiply charged fragments are scored for precursors of "
                       "charge 3 and above but not for charge 2, a pooled competition is dominated by the higher precursor charges. "
                       "With singly charged fragments only (deisotoped, high-resolution spectra) there is one group, as with 'pooled'. "
                       "With fragment:max_charge above 2, precursor charges of 4 and above form further groups, which accept few PSMs "
                       "per run (tens), so their q-values (D/T, as FalseDiscoveryRate) are coarse. "
                       "A group without decoy or without target PSMs falls back to 'pooled'.",
                       {"advanced"});
    defaults_.setValidStrings("FDR:PSM_groups", {"scored_charges", "pooled"});
    defaults_.setValue("FDR:protein", 0.0, "Filter proteins based on picked-protein FDR q-value (e.g., 0.01 = 1% protein FDR, set to 0 to disable). Applied after PSM-level FDR on a complete protein set (single file, or the -out_merged aggregate). Setting this > 0 finalizes the result: identified decoys are removed. With 0, decoys are retained for downstream/merged FDR. Uses the picked-protein approach (Savitski et al. 2015) which pairs target and decoy proteins by accession. The protein inference and the picked-protein FDR use the PSMs up to FDR:protein_inference_PSM, and a protein scores with the best posterior error probability of its PSMs (Percolator's, with 'rescore'; without it the best search engine score). The proteins and the protein groups (which have q-values of their own) are filtered at the threshold, and only the proteins of the groups that pass are kept. Proteins without a reported PSM (one that passes FDR:PSM) are removed afterwards, and so are the PSMs all of whose proteins failed the protein FDR. Needs decoys ('decoys' auto or generate).");
    defaults_.setMinFloat("FDR:protein", 0.0);
    defaults_.setMaxFloat("FDR:protein", 1.0);
    defaults_.setValue("FDR:protein_inference_PSM", 0.1,
                       "q-value cutoff of the PSMs that protein inference and the protein FDR (FDR:protein) use; FDR:PSM still "
                       "filters the reported PSMs. The q-values are those FDR:PSM filters on (Percolator's with 'rescore', else "
                       "see FDR:PSM_groups). A looser cutoff than FDR:PSM leaves more decoy proteins for the picked-protein FDR, "
                       "as ProteomicsLFQ infers proteins on PSMs at a looser q-value than its PSM FDR. Values at or below FDR:PSM "
                       "use the PSMs of FDR:PSM; with FDR:PSM 0 (no PSM filter) all PSMs are used.",
                       {"advanced"});
    defaults_.setMinFloat("FDR:protein_inference_PSM", 0.0);
    defaults_.setMaxFloat("FDR:protein_inference_PSM", 1.0);
    defaults_.setSectionDescription("FDR", "False Discovery Rate control (requires decoys)");

    // Add parameters which are only used by FragmentIndex
    defaults_.setValue("peptide:min_mass", 100, "Minimal peptide mass for database");
    defaults_.setValue("peptide:max_mass", 9000, "Maximal peptide mass for database");

    // Fragment-level filtering
    defaults_.setValue("fragment:min_matched_ions", 5, "Minimal number of matched ions to report a PSM");

    // Precursor isotope error handling. The search convention (theoretical minus observed) is kept for these two
    // parameters; the reported isotope_error follows report:isotope_error_convention.
    defaults_.setValue("precursor:isotope_error_min", -1,
                       "Minimum precursor isotope error searched, in 13C spacings (1.00336 Da) added to the observed "
                       "precursor mass: -1 finds a peptide whose first 13C isotope peak was selected as the precursor. "
                       "The PSM meta value 'isotope_error' reports the opposite sign by default (see "
                       "report:isotope_error_convention).");
    defaults_.setValue("precursor:isotope_error_max", 1,
                       "Maximum precursor isotope error searched, with the sign of precursor:isotope_error_min: +1 finds "
                       "a peptide whose precursor was selected one 13C spacing below its monoisotopic peak.");

    // Fragment and scoring limits
    defaults_.setValue("fragment:max_charge", 2, "max fragment charge");
    defaults_.setValue("scoring:max_candidates_per_spectrum", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue(
      "scoring:fragment_charges", "auto",
      "Fragment charges of the theoretical spectra that score each candidate. 'single' scores singly charged fragments only. "
      "'multiple' adds charges up to min(precursor charge - 1, fragment:max_charge). 'auto' uses 'multiple' when the MS2 spectra "
      "are not deisotoped (fragment:deisotope=false, or a fragment tolerance above the deisotoper's limit of 0.1 Da / 100 ppm, "
      "as for ion-trap CID) and 'single' otherwise, because deisotoping converts fragments to charge 1.",
      {"advanced"});
    defaults_.setValidStrings("scoring:fragment_charges", {"auto", "single", "multiple"});
    defaults_.setValue("scoring:method", "auto",
                       "Native score. 'hyperscore' is the X!Tandem HyperScore. 'mass_accuracy' weights each matched fragment's "
                       "contribution to HyperScore's ion counts and intensity product by exp(-0.5 * (error / scoring:mass_error_sd)^2), "
                       "with the error in ppm, so that matches near the theoretical m/z count fully and matches at the edge of the "
                       "tolerance count little; exact matches score as HyperScore. The kernel is centred at 0 and suits high-resolution "
                       "fragments only. 'auto' (default) uses 'mass_accuracy' for fragment tolerances in ppm within the high-resolution "
                       "range (<= 100 ppm) and 'hyperscore' otherwise (Da tolerances, e.g. ion-trap CID). Matching tolerances, candidates "
                       "and the unweighted match annotations are unchanged. The score type stays 'ln(hyperscore)'; the search parameters "
                       "record scoring:method_resolved when the weighted score is used.",
                       {"advanced"});
    defaults_.setValidStrings("scoring:method", {"hyperscore", "mass_accuracy", "auto"});
    defaults_.setValue("scoring:mass_error_sd", 7.0,
                       "Standard deviation in ppm of the Gaussian fragment mass-error kernel of scoring:method 'mass_accuracy'. "
                       "Fixed (not fitted) and independent of the matching tolerance and of calibration.",
                       {"advanced"});
    defaults_.setMinFloat("scoring:mass_error_sd", 1e-6);
    defaults_.setSectionDescription("scoring", "Search/Scoring Limits");

    // Ion series toggles
    defaults_.setValue("ions:add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("ions:add_y_ions", {"true","false"});
    defaults_.setValue("ions:add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("ions:add_b_ions", {"true","false"});
    defaults_.setValue("ions:add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("ions:add_a_ions", {"true","false"});
    defaults_.setValue("ions:add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("ions:add_c_ions", {"true","false"});
    defaults_.setValue("ions:add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("ions:add_x_ions", {"true","false"});
    defaults_.setValue("ions:add_z_ions", "false", "Add peaks of z-ions (y - NH3) to the spectrum. These are not the z+1 ions of ETD-type spectra; see ions:by_activation and ions:add_zp1_ions.");
    defaults_.setValidStrings("ions:add_z_ions", {"true","false"});
    defaults_.setValue("ions:add_zp1_ions", "false", "Add peaks of z+1 ions (z-dot, y - NH2) to the spectrum, the main C-terminal fragments of ETD, EThcD and ETciD spectra, for all spectra. ions:by_activation adds them (with c ions) for spectra recorded as electron-activated; set this, typically with ions:add_c_ions, for ETD-type data without activation information.");
    defaults_.setValidStrings("ions:add_zp1_ions", {"true","false"});
    defaults_.setValue("ions:by_activation", "true", "Score spectra whose precursor was activated by electrons (ETD, ECD, EThcD or ETciD, as recorded in the input file) with c and z+1 ions in addition to the ion series above. Other spectra, and spectra without activation information, use the ion series above.");
    defaults_.setValidStrings("ions:by_activation", {"true","false"});
    defaults_.setSectionDescription("ions", "Theoretical ion series toggles");

    defaults_.setValue("calibration:enabled", "false",
      "If enabled, run a fast calibration pass on a subset of spectra before the main search. "
      "Estimates tighter precursor and fragment tolerances from confident PSMs. "
      "The fragment index is NOT rebuilt — only query-time tolerances are tightened. "
      "Inspired by MSFragger's calibrate_mass and OpenNuXL's autotune.");
    defaults_.setValidStrings("calibration:enabled", {"true", "false"});
    defaults_.setValue("calibration:subset_ratio", 0.1,
      "Fraction of spectra (by TIC, highest first) used for the calibration pass (0.0-1.0).");
    defaults_.setMinFloat("calibration:subset_ratio", 0.01);
    defaults_.setMaxFloat("calibration:subset_ratio", 1.0);
    defaults_.setValue("calibration:min_psms", 50,
      "Minimum number of confident PSMs required for calibration. If fewer are found, "
      "calibration is skipped and the user-configured tolerances are used.");
    defaults_.setMinInt("calibration:min_psms", 1);

    defaults_.setValue("database:chunk_size", 0,
      "Split the protein database into chunks of at most this many proteins for fragment index "
      "building. 0 = disabled (load entire database at once). Enable for very large databases "
      "(e.g. MHC-II immunopeptidomics with variants) that exceed available memory. Each chunk "
      "builds its own fragment index; results are merged across chunks before post-processing. "
      "Calibration (when enabled) runs on a strided, size-bounded sample of the full "
      "decoy-augmented database — sample size is tied to chunk_size so calibration memory "
      "respects the same budget as main-search chunks. In multi-file mode (-in a.mzML b.mzML), "
      "the chunk-major path builds each chunk's fragment index once and scores all files "
      "against it before moving to the next chunk. Note that chunk_size bounds only the "
      "fragment-index memory: the chunk-major schedule holds ALL input files' preprocessed "
      "MS2 spectra in memory for the whole search, so budget for their sum as well, or split "
      "very large cohorts across separate ProSE invocations.");
    defaults_.setMinInt("database:chunk_size", 0);

    defaults_.setValue("chimeric:enabled", "false",
      "Also identify precursors co-isolated with the selected one (after the andes chimeric cascade). The survey scan of "
      "each MS2 spectrum is searched for isotope envelopes in its isolation window that fit averagine (needs centroided "
      "MS1 spectra in the input and isolation windows on the precursors). After the search, the peaks the best hit of the "
      "spectrum explains are removed, and the rest is searched for each co-isolated precursor in turn, at its m/z and "
      "charge; the best hit of each removes the peaks it explains before the next is searched. The best hit's peptide is "
      "not a hit of a co-isolated precursor, and a co-isolated precursor whose best hit repeats a PSM of the spectrum "
      "(same peptide and charge) is dropped. These PSMs get the spectrum reference of their spectrum and the PSM meta "
      "value 'chimeric_rank' 2 (the others 1). Not with database:chunk_size > 0.");
    defaults_.setValidStrings("chimeric:enabled", {"true", "false"});
    defaults_.setValue("chimeric:max_coisolated", 4, "At most this many co-isolated precursors per spectrum, the most intense first.");
    defaults_.setMinInt("chimeric:max_coisolated", 1);
    defaults_.setValue("chimeric:max_kl", 0.3,
      "Largest Kullback-Leibler divergence of an isotope envelope in the survey scan from averagine (first four isotope "
      "peaks) to accept it as a co-isolated precursor.");
    defaults_.setMinFloat("chimeric:max_kl", 0.0);
    defaults_.setValue("chimeric:top_hits", 3,
      "At most this many hits of a co-isolated precursor are reported (and at most report:top_hits). Each takes about as "
      "much memory as a hit of a spectrum.");
    defaults_.setMinInt("chimeric:top_hits", 1);
    defaults_.setValue("chimeric:mz_tolerance", 10.0, "Tolerance (ppm) for the m/z of the isotope peaks in the survey scan.");
    defaults_.setMinFloat("chimeric:mz_tolerance", 0.0);
    defaults_.setSectionDescription("chimeric", "Search of precursors co-isolated with the selected one (chimeric spectra).");

    defaults_.setSectionDescription("calibration",
      "Automatic mass accuracy calibration (two-pass search). A fast first pass on a subset of "
      "spectra estimates instrument-specific mass accuracy, then the main search uses the "
      "calibrated (typically tighter) tolerances for better discrimination.");

    FragmentIndex::restrictSearchParameters(defaults_); // the bounds of the parameters it shares with the FragmentIndex

    defaultsToParam_();
  }

  void ProSEAlgorithm::updateMembers_()
  {
    // inverted ranges and an invalid peptide:motif, before anything is read (the FragmentIndex checks them again)
    FragmentIndex::checkSearchParameters(param_);

    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();

    precursor_min_charge_ = param_.getValue("precursor:min_charge");
    precursor_max_charge_ = param_.getValue("precursor:max_charge");

    peptdeep_enable_ = param_.getValue("peptdeep:enable").toString() == "true";
    peptdeep_ms2_model_ = param_.getValue("peptdeep:ms2_model").toString();
    peptdeep_rt_model_ = param_.getValue("peptdeep:rt_model").toString();
    peptdeep_instrument_ = param_.getValue("peptdeep:instrument").toString();
    peptdeep_nce_ = param_.getValue("peptdeep:nce");
    peptdeep_rt_model_type_ = param_.getValue("peptdeep:rt_model_type").toString();
    peptdeep_batch_size_ = param_.getValue("peptdeep:batch_size");

    peaks_keep_n_ = (Size)(int)param_.getValue("peaks:keep_n");
    peaks_window_top_ = (Int)param_.getValue("peaks:window_top");
    peaks_window_type_ = param_.getValue("peaks:window_type").toString();

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");
    // Every scoring method and the fragment index match within this tolerance (0, negative and NaN tolerances match
    // nothing, infinity everything), and HyperScore::computeMassAccuracy() and local fragment evidence require it.
    // Checked here, before the parallel scoring loops.
    if (! std::isfinite(fragment_mass_tolerance_) || fragment_mass_tolerance_ <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "fragment:mass_tolerance must be finite and positive (got " + std::to_string(fragment_mass_tolerance_) + ").");
    }

    fragment_mass_tolerance_unit_ = param_.getValue("fragment:mass_tolerance_unit").toString();

    // Resolve the MS2 deisotoping mode (fragment:deisotope = auto|true|false) ONCE.
    // The Deisotoper supports a fragment tolerance <= 100 ppm / <= 0.1 Da and throws
    // otherwise; consult its non-throwing predicate so the limit lives in one place
    // (OpenMS#9619). 'false' never deisotopes; 'true'/'auto' deisotope only when the
    // tolerance is supported -- 'true' fails fast here with a clear parameter error
    // rather than aborting later inside preprocessSpectra_'s OpenMP region, while
    // 'auto' silently skips deisotoping for low-resolution data.
    const std::string deisotope_mode = param_.getValue("fragment:deisotope").toString();
    const bool deisotope_supported =
      Deisotoper::isToleranceSupported(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm");
    // Resolved from the configured tolerance, with the deisotoper's resolution boundary; calibration does not switch it.
    const std::string query_spectrum = param_.getValue("fragment:query_spectrum").toString();
    query_raw_spectrum_ = query_spectrum == "raw" || (query_spectrum == "auto" && deisotope_supported);
    deisotope_requested_ = (deisotope_mode != "false");
    if (deisotope_mode == "true" && !deisotope_supported)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "fragment:deisotope=true requires a high-resolution fragment tolerance (<= 0.1 Da or <= 100 ppm). "
        "Use 'auto' (deisotopes only for high-resolution data) or 'false' for low-resolution data.");
    }
    if (deisotope_mode == "auto" && !deisotope_supported)
    {
      OPENMS_LOG_WARN << "[ProSE] Fragment tolerance " << fragment_mass_tolerance_
                      << (fragment_mass_tolerance_unit_ == "ppm" ? " ppm" : " Da")
                      << " exceeds the deisotoping limit (100 ppm / 0.1 Da); skipping MS2 "
                      << "deisotoping (expected for low-resolution data)." << endl;
    }
    deisotoping_.min_peaks = static_cast<unsigned int>(static_cast<int>(param_.getValue("fragment:deisotope_min_peaks")));
    deisotoping_.charge_cap_precursor = param_.getValue("fragment:deisotope_charge_cap").toString() == "precursor";
    deisotoping_.sum_intensity = param_.getValue("fragment:deisotope_sum_intensity").toBool();

    // Spectra that are not deisotoped keep multiply charged fragments at their own m/z;
    // deisotoped spectra hold charge-1 fragments only.
    const std::string fragment_charges = param_.getValue("scoring:fragment_charges").toString();
    const bool deisotoped = deisotope_requested_ && deisotope_supported;
    scoring_multiple_charges_ = fragment_charges == "multiple" || (fragment_charges == "auto" && ! deisotoped);
    scoring_max_charge_ = static_cast<int>(param_.getValue("fragment:max_charge"));

    // Mass-accuracy weighting needs fragment errors of a few ppm; 'auto' uses it for high-resolution ppm tolerances only.
    const std::string scoring_method = param_.getValue("scoring:method").toString();
    const bool high_resolution_ppm = fragment_mass_tolerance_unit_ == "ppm" && deisotope_supported;
    mass_accuracy_score_ = scoring_method == "mass_accuracy" || (scoring_method == "auto" && high_resolution_ppm);
    mass_error_sd_ppm_ = param_.getValue("scoring:mass_error_sd");
    if (mass_accuracy_score_ && (! std::isfinite(mass_error_sd_ppm_) || mass_error_sd_ppm_ <= 0.0))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "scoring:mass_error_sd must be finite and positive.");
    }
    if (scoring_method == "mass_accuracy" && ! high_resolution_ppm)
    {
      OPENMS_LOG_WARN << "[ProSE] scoring:method=mass_accuracy weights fragment matches by a " << mass_error_sd_ppm_
                      << " ppm kernel, but the fragment tolerance is " << fragment_mass_tolerance_
                      << (fragment_mass_tolerance_unit_ == "ppm" ? " ppm" : " Da")
                      << "; matches far from the theoretical m/z score little. 'auto' uses HyperScore for such tolerances." << endl;
    }

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    set<std::string> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());
    if (fixed_unique.size() != modifications_fixed_.size())
    {
      OPENMS_LOG_WARN << "Duplicate fixed modification provided. Making them unique." << endl;
      modifications_fixed_.assign(fixed_unique.begin(), fixed_unique.end());
    }
    // fixed terminal modifications the fragment index cannot restrict: fail here, not after reading the database
    FragmentIndex::checkFixedModifications(modifications_fixed_);

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<std::string> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      OPENMS_LOG_WARN << "Duplicate variable modification provided. Making them unique." << endl;
      modifications_variable_.assign(var_unique.begin(), var_unique.end());
    }
    // A modification that is fixed as well applies everywhere already: the remedy is to drop the fixed one.
    const auto also_fixed = [this](const std::string& mod)
    {
      return std::find(modifications_fixed_.begin(), modifications_fixed_.end(), mod) != modifications_fixed_.end();
    };
    for (const std::string& mod : FragmentIndex::shadowedVariableTerminalModifications(modifications_fixed_, modifications_variable_))
    {
      if (also_fixed(mod))
      {
        OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: it is also a fixed modification, which "
                        << "applies to every such terminus. To search the terminus with and without it, remove it from "
                        << "modifications:fixed." << endl;
        continue;
      }
      OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: a fixed modification already sits on "
                      << "that terminus, which carries one modification. To search it, specify the fixed terminal "
                      << "modification as a variable one as well." << endl;
    }
    for (const std::string& mod : FragmentIndex::shadowedVariableResidueModifications(modifications_fixed_, modifications_variable_))
    {
      if (also_fixed(mod))
      {
        OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: it is also a fixed modification, which "
                        << "applies to every such residue. To search the residue with and without it, remove it from "
                        << "modifications:fixed." << endl;
        continue;
      }
      OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: a fixed modification already sits on "
                      << "that residue, which carries one modification. To search it, specify the fixed residue "
                      << "modification as a variable one instead." << endl;
    }

    modifications_max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    enzyme_ = param_.getValue("enzyme").toString();

    peptide_min_size_ = param_.getValue("peptide:min_size");
    peptide_max_size_ = param_.getValue("peptide:max_size");
    peptide_missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_enzyme_specificity_ = EnzymaticDigestion::getSpecificityByName(
      param_.getValue("peptide:enzyme_specificity").toString());

    report_top_hits_ = param_.getValue("report:top_hits");
    rescore_top_hits_ = static_cast<Size>(static_cast<int>(param_.getValue("rescoring:top_hits")));
    isotope_error_observed_minus_theoretical_ =
      param_.getValue("report:isotope_error_convention").toString() == "observed_minus_theoretical";

    const std::string decoy_mode_str = param_.getValue("decoys").toString();
    if (decoy_mode_str == "generate")   { decoy_mode_ = DecoyMode_::GENERATE; }
    else if (decoy_mode_str == "ignore") { decoy_mode_ = DecoyMode_::IGNORE; }
    else                                 { decoy_mode_ = DecoyMode_::AUTO; }
    decoy_prefix_ = param_.getValue("decoy_prefix").toString();
    decoy_keep_termini_ = param_.getValue("decoy_method").toString() == "keep_termini";
    decoy_parity_report_ = param_.getValue("decoy_parity_report").toBool();
    terminal_features_ = param_.getValue("rescoring:terminal_features").toString() == "true";
    stratum_q_values_ = param_.getValue("rescoring:stratum_q_values").toString() == "true";
    annotate_psm_ = ListUtils::toStringList<std::string>(param_.getValue("annotate:PSM"));
    if (peptdeep_enable_
        && std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") == annotate_psm_.end()
        && std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) == annotate_psm_.end())
    {
      // PeptDeepRescoring compares the predicted with the annotated fragment peaks of each PSM
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "peptdeep:enable needs the fragment annotation of the PSMs: annotate:PSM must contain 'ALL' or '"
        + std::string(Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) + "'.");
    }
    self_trained_ion_priors_ = param_.getValue("annotate:self_trained_ion_priors").toBool();
    ion_prior_model_ = param_.getValue("annotate:ion_prior_model").toString() == "basic"
                         ? FragmentIonLikelihoodModel::ContextSet::BASIC : FragmentIonLikelihoodModel::ContextSet::RICH;
    ion_prior_scored_peaks_ = param_.getValue("annotate:ion_prior_peaks").toString() == "scored";
    ion_prior_max_fragment_charge_ = static_cast<int>(param_.getValue("annotate:ion_prior_max_fragment_charge"));
    ion_prior_train_fdr_ = param_.getValue("annotate:ion_prior_train_fdr");
    ion_prior_min_psms_ = static_cast<Size>(static_cast<int>(param_.getValue("annotate:ion_prior_min_psms")));
    {
      const StringList features = ListUtils::toStringList<std::string>(param_.getValue("annotate:ion_prior_features"));
      auto selected = [&features](const std::string& name) { return std::find(features.begin(), features.end(), name) != features.end(); };
      ion_prior_feature_llr_ = selected(Constants::UserParam::ION_PRIOR_LLR);
      ion_prior_feature_explained_ = selected(Constants::UserParam::ION_PRIOR_EXPLAINED);
      ion_prior_feature_topk_ = selected(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED);
      if (self_trained_ion_priors_ && !(ion_prior_feature_llr_ || ion_prior_feature_explained_ || ion_prior_feature_topk_))
      {
        OPENMS_LOG_WARN << "[ProSE] annotate:self_trained_ion_priors is set, but annotate:ion_prior_features is empty: "
                        << "no ion priors are learned." << endl;
        self_trained_ion_priors_ = false;
      }
    }
    if (self_trained_ion_priors_ && decoy_mode_ == DecoyMode_::IGNORE)
    {
      // A target-only search has no target-decoy competition to select training PSMs, and Percolator cannot use its
      // output: nothing is learned or written (the output equals the one with the ion priors switched off).
      OPENMS_LOG_INFO << "[ProSE] annotate:self_trained_ion_priors: not applied to a target-only search (decoys=ignore)." << endl;
      self_trained_ion_priors_ = false;
    }
    rescore_ = param_.getValue("rescore").toString() == "true";
    max_train_psms_ = static_cast<Size>(static_cast<int>(param_.getValue("percolator:max_train_psms")));
    percolator_training_.subset_max_train = static_cast<int>(std::min<Size>(max_train_psms_, std::numeric_limits<int>::max()));
    percolator_training_.subset_mode = param_.getValue("rescoring:subset_mode").toString();
    percolator_training_.curve_start = static_cast<int>(param_.getValue("rescoring:curve_start"));
    percolator_training_.curve_min_gain = static_cast<double>(param_.getValue("rescoring:curve_min_gain"));
    percolator_training_.iterations = static_cast<int>(param_.getValue("rescoring:iterations"));
    percolator_training_.tron_epsilon = static_cast<double>(param_.getValue("rescoring:tron_epsilon"));
    percolator_training_.warm_start = param_.getValue("rescoring:warm_start").toString() == "true";
    percolator_training_.single_cneg = param_.getValue("rescoring:single_cneg").toString() == "true";
    percolator_training_.solver = param_.getValue("rescore_solver").toString();
    percolator_training_.use_pi0 = param_.getValue("rescoring:use_pi0").toString();
    percolator_training_.calc_mass_feature = param_.getValue("rescoring:calc_mass_feature").toString() == "true";
    percolator_training_.fixed_enzyme_termini = param_.getValue("rescoring:fixed_enzyme_termini").toString() == "true";
    percolator_training_.peptide_level = param_.getValue("rescoring:peptide_level").toString() == "true";
    percolator_extra_ = Param();
    if (param_.getValue("rescore_method").toString() == "reset")
    {
      percolator_extra_.setValue("algorithm", "reset");
      for (const std::string key : {"salt", "learner", "audit", "stem_table", "representative", "mlp_crossfit", "decoys",
                                    "second_decoy_prefix"})
      {
        percolator_extra_.setValue("reset:" + key, param_.getValue("reset:" + key).toString());
      }
      for (const std::string key : {"iterations", "mlp_members", "mlp_hidden", "mlp_epochs", "mlp_warm_epochs", "mlp_norm_min_decoys"})
      {
        percolator_extra_.setValue("reset:" + key, static_cast<int>(param_.getValue("reset:" + key)));
      }
      if (static_cast<double>(param_.getValue("FDR:PSM")) > 0.0)
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_method 'reset' with FDR:PSM > 0: the PSM q-values that FDR:PSM filters on are "
                           "heuristic under RESET (FDR is controlled at the peptide level: RESET:peptide_q_value)." << endl;
      }
    }
    reuse_mode_ = param_.getValue("rescore_reuse:mode").toString();
    reuse_max_train_ = static_cast<Size>(static_cast<int>(param_.getValue("rescore_reuse:max_train")));
    reuse_static_files_ = static_cast<Size>(static_cast<int>(param_.getValue("rescore_reuse:static_files")));
    reuse_max_shift_ = param_.getValue("rescore_reuse:max_shift");
    if (percolator_extra_.exists("algorithm") && reuse_mode_ != "per_file")
    {
      OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode '" << reuse_mode_ << "' does not apply with rescore_method 'reset' (the "
                         "reused models are Percolator models): every file is rescored by RESET ('per_file')." << endl;
      reuse_mode_ = "per_file";
    }
    fdr_psm_ = param_.getValue("FDR:PSM");
    fdr_psm_by_scored_charges_ = param_.getValue("FDR:PSM_groups").toString() == "scored_charges";
    fdr_protein_ = param_.getValue("FDR:protein");
    fdr_protein_inference_psm_ = param_.getValue("FDR:protein_inference_PSM");

    // Open search mode is automatically determined based on precursor tolerance in isOpenSearchMode_()

    add_a_ions_ = param_.getValue("ions:add_a_ions").toBool();
    add_b_ions_ = param_.getValue("ions:add_b_ions").toBool();
    add_c_ions_ = param_.getValue("ions:add_c_ions").toBool();
    add_x_ions_ = param_.getValue("ions:add_x_ions").toBool();
    add_y_ions_ = param_.getValue("ions:add_y_ions").toBool();
    add_z_ions_ = param_.getValue("ions:add_z_ions").toBool();
    add_zp1_ions_ = param_.getValue("ions:add_zp1_ions").toBool();
    ions_by_activation_ = param_.getValue("ions:by_activation").toBool();
    protein_mapping_from_index_ = param_.getValue("peptide:protein_mapping").toString() == "index";

    database_chunk_size_ = param_.getValue("database:chunk_size");
    chimeric_enabled_ = param_.getValue("chimeric:enabled").toBool();

    calibration_enabled_ = param_.getValue("calibration:enabled") == "true";
    calibration_subset_ratio_ = param_.getValue("calibration:subset_ratio");
    calibration_min_psms_ = param_.getValue("calibration:min_psms");
  }

  // static
  struct ProSEAlgorithm::SpectrumGenerators_
  {
    TheoreticalSpectrumGenerator standard; ///< configured ion series
    TheoreticalSpectrumGenerator electron; ///< configured ion series plus c and z+1 ions
    bool by_activation = false;

    /// True if @p spectrum is matched against, and scored with, c and z+1 ions as well: its
    /// precursor was activated by electrons
    bool electronIons(const MSSpectrum& spectrum) const
    {
      return by_activation && isElectronActivated_(spectrum);
    }

    /// The generator for @p spectrum, chosen by the activation of its precursor
    const TheoreticalSpectrumGenerator& forSpectrum(const MSSpectrum& spectrum) const
    {
      return electronIons(spectrum) ? electron : standard;
    }
  };

  ProSEAlgorithm::SpectrumGenerators_ ProSEAlgorithm::spectrumGenerators_() const
  {
    // Scoring, annotation and calibration use the ion series of the fragment index, so that
    // a candidate is scored against the ions it was retrieved by.
    auto configure = [this](TheoreticalSpectrumGenerator& tsg, bool electron)
    {
      Param p(tsg.getParameters());
      p.setValue("add_first_prefix_ion", "true");
      p.setValue("add_metainfo", "true");
      p.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
      p.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
      p.setValue("add_c_ions", (add_c_ions_ || electron) ? "true" : "false");
      p.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
      p.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
      p.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
      p.setValue("add_zp1_ions", (add_zp1_ions_ || electron) ? "true" : "false");
      tsg.setParameters(p);
    };
    SpectrumGenerators_ generators;
    configure(generators.standard, false);
    configure(generators.electron, true);
    generators.by_activation = ions_by_activation_;
    return generators;
  }

  std::vector<double> ProSEAlgorithm::localPeakDensities_(const MSSpectrum& spectrum)
  {
    // Both bounds only move forward: O(number of peaks), once per spectrum.
    constexpr double half_window = 50.0;
    std::vector<double> densities(spectrum.size());
    Size left = 0, right = 0;
    for (Size i = 0; i < spectrum.size(); ++i)
    {
      const double mz = spectrum[i].getMZ();
      while (left < spectrum.size() && spectrum[left].getMZ() < mz - half_window) ++left;
      while (right < spectrum.size() && spectrum[right].getMZ() <= mz + half_window) ++right;
      densities[i] = static_cast<double>(right - left) / (2.0 * half_window);
    }
    return densities;
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical,
      const std::vector<double>& densities, double tolerance, bool ppm)
  {
    std::vector<double> alternatives;
    return localFragmentEvidence_(spectrum, theoretical, std::numeric_limits<int>::max(), densities, tolerance, ppm, alternatives);
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical, int max_charge,
      const std::vector<double>& densities, double tolerance, bool ppm, std::vector<double>& alternatives)
  {
    if (spectrum.empty() || theoretical.empty()) return LocalFragmentEvidence_();
    OPENMS_PRECONDITION(!theoretical.getIntegerDataArrays().empty()
      && theoretical.getIntegerDataArrays()[0].size() == theoretical.size(), "Aligned fragment charges required")
    return localFragmentEvidence_(spectrum, theoretical, theoretical.getIntegerDataArrays()[0], max_charge, densities, tolerance,
                                  ppm, alternatives);
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical, const std::vector<Int>& charges, int max_charge,
      const std::vector<double>& densities, double tolerance, bool ppm, std::vector<double>& alternatives)
  {
    // The ions of charge <= max_charge of a TSG spectrum of charges 1..Z are the TSG spectrum of charges
    // 1..max_charge, in the same order: the generator adds each charge independently and sorts stably.
    // Annotation therefore passes its own theoretical spectrum instead of generating a second one.
    LocalFragmentEvidence_ evidence;
    if (spectrum.empty() || theoretical.empty()) return evidence;
    OPENMS_PRECONDITION(densities.size() == spectrum.size(), "One local density per experimental peak required")
    OPENMS_PRECONDITION(charges.size() == theoretical.size(), "Aligned fragment charges required")

    static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
    static const double ammonia_mass = EmpiricalFormula("NH3").getMonoWeight();
    alternatives.clear();
    alternatives.reserve(3 * theoretical.size());
    for (Size i = 0; i < theoretical.size(); ++i)
    {
      OPENMS_PRECONDITION(charges[i] > 0, "Positive fragment charge required")
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      alternatives.push_back(mz);
      for (const double loss : {water_mass, ammonia_mass})
      {
        const double loss_mz = mz - loss / charges[i];
        if (loss_mz > 0.0) alternatives.push_back(loss_mz);
      }
    }
    std::sort(alternatives.begin(), alternatives.end());

    for (Size i = 0; i < theoretical.size(); ++i)
    {
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      const double width = ppm ? mz * tolerance * 1e-6 : tolerance;
      const Size nearest = spectrum.findNearest(mz);
      const double observed = spectrum[nearest].getMZ();
      if (std::abs(mz - observed) > width) continue;

      // The matched peak itself is in its density window, hence density > 0.
      const double density = densities[nearest];
      evidence.chance_match_surprise += std::max(0.0, -std::log(2.0 * width * density));
      const auto first = std::lower_bound(alternatives.begin(), alternatives.end(), observed - width);
      const auto last = std::upper_bound(first, alternatives.end(), observed + width);
      // One intact hypothesis is the ion currently being scored. All other
      // hypotheses, including coincident ions, represent alternative assignments.
      const Size count = static_cast<Size>(last - first);
      const Size competing = count > 0 ? count - 1 : 0;
      evidence.mass_competition_evidence += 1.0 / (1.0 + competing + density);
    }
    return evidence;
  }

  bool ProSEAlgorithm::isElectronActivated_(const MSSpectrum& spectrum)
  {
    if (spectrum.getPrecursors().empty()) return false;
    for (const Precursor::ActivationMethod method : spectrum.getPrecursors()[0].getActivationMethods())
    {
      if (method == Precursor::ActivationMethod::ETD || method == Precursor::ActivationMethod::ECD
          || method == Precursor::ActivationMethod::EThcD || method == Precursor::ActivationMethod::ETciD)
      {
        return true;
      }
    }
    return false;
  }

  Size ProSEAlgorithm::countElectronActivated_(const PeakMap& spectra) const
  {
    if (!ions_by_activation_) return 0;
    return static_cast<Size>(std::count_if(spectra.begin(), spectra.end(),
                                           [](const MSSpectrum& spectrum) { return isElectronActivated_(spectrum); }));
  }

  namespace
  {
    // Calls f when the scope ends, also when an exception leaves it (f must not throw then).
    template <typename F>
    class OnScopeExit
    {
    public:
      explicit OnScopeExit(F f) : f_(std::move(f)) {}
      OnScopeExit(const OnScopeExit&) = delete;
      OnScopeExit& operator=(const OnScopeExit&) = delete;
      ~OnScopeExit() { f_(); }

    private:
      F f_;
    };

    // Background reads of spectra that lend OpenMP threads of this thread (BackgroundSpectra, below):
    // takeBackThreadsOfEndedRead() gives this thread its threads back if the read has ended (called between the phases
    // of a search, outside of parallel regions); startDeferredRead() starts a read that waits for the rescoring of the
    // search before it; readLendsThreads(): whether a read lends threads of this thread now.
    void takeBackThreadsOfEndedRead();
    void startDeferredRead();
    bool readLendsThreads();

    // Keeps the peaks WindowMower::filterPeakSpectrumForTopNInJumpingWindow() keeps: the peak_count most intense peaks of
    // each m/z window of window_size that starts at a peak (of the last window a share of peak_count by its width), and
    // every peak equal to a kept one. It selects them through indices instead of copies of the peaks, and finds the
    // equal peaks among their neighbours instead of searching all kept peaks for each peak. The same std::partial_sort
    // calls on the same intensities select the same peaks.
    void filterTopNInJumpingWindows(MSSpectrum& spectrum, double window_size, UInt peak_count)
    {
      if (spectrum.empty()) { return; }
      spectrum.sortByPosition();
      const Size n = spectrum.size();
      std::vector<char> kept(n, 0);
      std::vector<Size> window;
      const auto more_intense = [&spectrum](Size a, Size b) { return spectrum[b].getIntensity() < spectrum[a].getIntensity(); };
      const auto keep_most_intense = [&](Size begin, Size end, Size count)
      {
        if (end - begin > count)
        {
          window.resize(end - begin);
          std::iota(window.begin(), window.end(), begin);
          std::partial_sort(window.begin(), window.begin() + count, window.end(), more_intense);
          for (Size k = 0; k < count; ++k) { kept[window[k]] = 1; }
        }
        else
        {
          std::fill(kept.begin() + begin, kept.begin() + end, 1);
        }
      };
      double window_start = spectrum[0].getMZ();
      Size begin = 0;
      for (Size i = 0; i != n; ++i)
      {
        if (spectrum[i].getMZ() - window_start < window_size) { continue; }
        // a gap may leave windows empty: the next window starts at the next peak
        window_start = spectrum[i].getMZ();
        keep_most_intense(begin, i, peak_count);
        begin = i;
      }
      const double last_window_fraction = (spectrum[n - 1].getMZ() - window_start) / window_size;
      keep_most_intense(begin, n, static_cast<Size>(std::round(last_window_fraction * peak_count)));

      // peaks of equal m/z are neighbours (the spectrum is sorted by m/z)
      std::vector<Size> selected;
      selected.reserve(n);
      for (Size i = 0; i < n;)
      {
        Size j = i + 1;
        while (j < n && spectrum[j].getMZ() == spectrum[i].getMZ()) { ++j; }
        for (Size a = i; a < j; ++a)
        {
          for (Size b = i; b < j; ++b)
          {
            if (kept[b] && spectrum[b] == spectrum[a])
            {
              selected.push_back(a);
              break;
            }
          }
        }
        i = j;
      }
      spectrum.select(selected);
    }
  }

  void ProSEAlgorithm::filterLocalPeaks_(MSSpectrum& spectrum, Size peaks_per_window)
  {
    // Work on indices so all peak-associated data arrays survive the selection.
    // An isolated high-m/z ion is still evidence; its window's observed width
    // must not round the retention quota down to zero.
    std::vector<Size> indices(spectrum.size()), selected;
    std::iota(indices.begin(), indices.end(), 0);
    selected.reserve(spectrum.size());
    for (Size begin = 0; begin < spectrum.size();)
    {
      Size end = begin + 1;
      while (end < spectrum.size() && spectrum[end].getMZ() - spectrum[begin].getMZ() < 100.0)
      {
        ++end;
      }
      const Size keep = std::min(peaks_per_window, end - begin);
      std::partial_sort(indices.begin() + begin, indices.begin() + begin + keep, indices.begin() + end, [&spectrum](Size a, Size b) {
        if (spectrum[a].getIntensity() != spectrum[b].getIntensity()) { return spectrum[a].getIntensity() > spectrum[b].getIntensity(); }
        return a < b;
      });
      selected.insert(selected.end(), indices.begin() + begin, indices.begin() + begin + keep);
      begin = end;
    }
    std::sort(selected.begin(), selected.end());
    spectrum.select(selected);
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type)
  {
    preprocessSpectra_(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, deisotope_requested, peaks_keep_n,
                       peaks_window_top, window_type, DeisotopingSettings_{});
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type,
                                          const DeisotopingSettings_& deisotoping,
                                          FragmentIonLikelihoodModel::PeakLists* ion_evidence,
                                          bool ion_evidence_scored_peaks,
                                          PeakMap* evidence_spectra,
                                          PeakMap* query_spectra)
  {
    // Intensity threshold + normalization used to run here as two extra SERIAL full-map
    // passes. Both are strictly per-spectrum: ThresholdMower::filterPeakMap and
    // Normalizer::filterPeakMap are literally "for (auto& s : exp) filterSpectrum(s);" and
    // neither iterates chromatograms. They are therefore applied at the top of the parallel
    // loop below instead, which is per-spectrum equivalent and removes two full sweeps over
    // the peak data. Both objects are configured once here; like nlargest_filter below, each
    // OpenMP thread works on its own copy (made before the loop): ThresholdMower stores its 'threshold'
    // Param in a member on every call, and concurrent writes are a data race even when every
    // thread writes the same value. One copy per thread costs a few Param copies per search.
    // Peaks without intensity (zero or negative, e.g. empty centroids) would still count as
    // matched ions, so they are removed. Nothing else is (every positive float intensity,
    // subnormal ones included, passes): the ThresholdMower default of 0.05 is an absolute
    // cutoff on the raw intensities before normalization and deleted real peaks from
    // intensity-scaled input (e.g. spectra normalized to a base peak of 1).
    ThresholdMower threshold_mower_filter;
    Param threshold_param = threshold_mower_filter.getParameters();
    threshold_param.setValue("threshold", static_cast<double>(std::numeric_limits<float>::denorm_min()));
    threshold_mower_filter.setParameters(threshold_param);
    Normalizer normalizer;

    // sort by rt; done before the loop because sortSpectra(false) only permutes whole
    // spectra (it does not touch peak order, see MSExperiment::sortSpectra) and therefore
    // commutes with the per-spectrum filters that now run inside the loop.
    if (!inRetentionTimeOrder(exp)) { exp.sortSpectra(false); } // the spectra of loadMS2Spectra() are sorted already
    if (ion_evidence != nullptr) { ion_evidence->reset(exp.size()); }
    if (evidence_spectra != nullptr)
    {
      evidence_spectra->clear(true);
      resizeWithParallelFaults(evidence_spectra->getSpectra(), exp.size());
    }
    if (query_spectra != nullptr)
    {
      query_spectra->clear(true);
      resizeWithParallelFaults(query_spectra->getSpectra(), exp.size());
    }

    // filter settings: the most intense peaks_window_top peaks per 100 Th window (jumping windows as WindowMower's,
    // unless full_window_quota)
    const bool full_window_quota
      = window_type == "jump_full"
        || (window_type == "auto" && Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm));

    // Resolution-aware peak retention. peaks_keep_n == 0 => auto. For HIGH-resolution fragments
    // (within the deisotoper range, <= 0.1 Da / <= 100 ppm) keep the legacy cap of 400. For
    // LOW-resolution data (e.g. 0.5 Da ion-trap CID — the same regime where deisotoping is
    // skipped) a wide match window admits many spurious low-intensity matches that inflate the
    // count-sensitive HyperScore for targets AND decoys alike, collapsing the target/decoy
    // margin; retain far fewer peaks. Random matches scale with peak density x tolerance width,
    // so sqrt-dampen around the 0.02 Da reference and clamp to [60, 400]:
    //   0.2 Da -> 126, 0.5 Da -> 80, 1 Da -> 60.
    Size effective_keep_n = peaks_keep_n;
    if (effective_keep_n == 0)
    {
      if (Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm))
      {
        effective_keep_n = 400; // high-resolution: unchanged legacy behavior
      }
      else
      {
        const double eff_tol_da = fragment_mass_tolerance_unit_ppm
          ? fragment_mass_tolerance * 500.0 * 1e-6   // ppm -> Da at a 500 m/z reference
          : fragment_mass_tolerance;
        const double n = std::round(400.0 * std::sqrt(0.02 / std::max(eff_tol_da, 1e-6)));
        effective_keep_n = static_cast<Size>(std::min(400.0, std::max(60.0, n)));
      }
    }
    NLargest nlargest_filter = NLargest(effective_keep_n);

    // Deisotope only when requested (param fragment:deisotope, resolved in
    // updateMembers_) AND the tolerance is in the Deisotoper's supported range.
    // The range guard keeps the Deisotoper call below from ever throwing, so no
    // exception can escape the OpenMP region (OpenMS#9619). Mode resolution and the
    // low-resolution warning are handled once, in updateMembers_.
    const bool do_deisotope = deisotope_requested &&
      Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);

    // The per-thread copies of the filters are made here: copied on entry to the parallel region (firstprivate), they
    // were outside any exception handling, and an exception from a copy (std::bad_alloc) terminated the process.
#ifdef _OPENMP
    const Size num_threads = static_cast<Size>(omp_get_max_threads());
#else
    const Size num_threads = 1;
#endif
    std::vector<ThresholdMower> threshold_mower_filters(num_threads, threshold_mower_filter);
    std::vector<Normalizer> normalizers(num_threads, normalizer);
    std::vector<NLargest> nlargest_filters(num_threads, nlargest_filter);

    // An exception on a worker thread must not leave the parallel region, which would terminate the process: the
    // first one is rethrown after the loop.
    std::exception_ptr preprocessing_error;
    std::atomic<bool> preprocessing_failed{false};
    // (static: freeInParallel() frees the copies made here in the same parts, each in its own thread's arena)
#pragma omp parallel for schedule(static) default(none) shared(exp, evidence_spectra, query_spectra, do_deisotope, fragment_mass_tolerance, \
                                                fragment_mass_tolerance_unit_ppm, full_window_quota, peaks_window_top, \
                                                deisotoping, ion_evidence, ion_evidence_scored_peaks, \
                                                threshold_mower_filters, normalizers, nlargest_filters, \
                                                preprocessing_error, preprocessing_failed)
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      if (preprocessing_failed.load(std::memory_order_relaxed)) continue;
      try
      {
#ifdef _OPENMP
      const Size thread = static_cast<Size>(omp_get_thread_num());
#else
      const Size thread = 0;
#endif
      // remove 0 intensities, then normalize (formerly two serial full-map passes)
      threshold_mower_filters[thread].filterPeakSpectrum(exp[exp_index]);
      normalizers[thread].filterPeakSpectrum(exp[exp_index]);

      // sort by mz
      exp[exp_index].sortByPosition();

      // fragment:query_spectrum=raw retrieves candidates with this peak list (before deisotoping and local filtering).
      // (only what the candidate query reads: the peaks, MS level and precursors; not the meta data and data arrays)
      if (query_spectra != nullptr) { copyPeakList_(exp[exp_index], (*query_spectra)[exp_index], true); }

      // deisotope (skipped for low-resolution data; see do_deisotope above)
      // Isotope intensities must fall from the monoisotopic peak on (start_intensity_check = 1).
      // With the library default of 2, a small peak one isotope spacing below a fragment ion
      // became the envelope's monoisotopic peak and the ion itself was removed as its isotope.
      // TMT/TMTpro-labelled fragments carry such a peak (reagent isotope impurity), and dense
      // Orbitrap Astral and timsTOF spectra often hold one by chance.
      // The envelope rule is set by fragment:deisotope_min_peaks, _charge_cap and _sum_intensity; their defaults
      // are Sage-like (two peaks, charges up to the precursor charge but at most 3, summed intensity; the earlier
      // rule was three peaks, charges 1-3, own intensity). The Sage-like rule differs from Sage itself: Sage links isotope
      // pairs with a fixed 10 ppm tolerance, requires each isotope peak to be weaker than its parent and does not
      // cap the charge at 3, while the OpenMS deisotoper extends an envelope from its monoisotopic peak with the
      // search tolerance and stops at the first isotope peak that is more intense than its predecessor.
      // The deisotoper does not exclude peaks that already belong to an earlier envelope, so an isotope peak can be
      // claimed by two monoisotopic peaks and, with summing, counted in both (two-peak envelopes make this more
      // frequent; Sage likewise adds a peak to several parents of the same charge). It does not depend on labels.
      if (do_deisotope)
      {
        int max_charge = 3;
        if (deisotoping.charge_cap_precursor && ! exp[exp_index].getPrecursors().empty())
        {
          const int precursor_charge = exp[exp_index].getPrecursors()[0].getCharge();
          if (precursor_charge > 0) { max_charge = std::min(max_charge, precursor_charge); }
        }
        Deisotoper::deisotopeAndSingleCharge(exp[exp_index],
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
          1, max_charge,  // min / max charge
          false,  // keep only deisotoped
          deisotoping.min_peaks, 10,  // min / max isopeaks
          true,   // convert fragment m/z to mono-charge
          false,  // annotate charge
          false,  // annotate isotopic peak counts
          true,   // decreasing isotope intensities
          1,      // start the intensity check at the monoisotopic peak
          deisotoping.sum_intensity);  // the monoisotopic peak carries the envelope's intensity
      }

      // ion priors: every peak after deisotoping, before the window and top-N filters
      if (ion_evidence != nullptr && !ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
      // Local fragment evidence uses this peak list: after deisotoping, so that densities and
      // matching share the m/z space of the scoring spectrum, but before local and top-N filtering.
      // (only what localPeakDensities_ and localFragmentEvidence_ read: the peaks)
      if (evidence_spectra != nullptr) { copyPeakList_(exp[exp_index], (*evidence_spectra)[exp_index], false); }

      // remove noise
      if (full_window_quota) { filterLocalPeaks_(exp[exp_index], static_cast<Size>(peaks_window_top)); }
      else { filterTopNInJumpingWindows(exp[exp_index], 100.0, static_cast<UInt>(peaks_window_top)); }
      nlargest_filters[thread].filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();

      // ion priors on the scored peaks
      if (ion_evidence != nullptr && ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_preprocessing_error)
        {
          if (!preprocessing_error) preprocessing_error = std::current_exception();
        }
        preprocessing_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (preprocessing_error) std::rethrow_exception(preprocessing_error);
  }

  double ProSEAlgorithm::CandidatePoolStats_::zScore() const
  {
    if (count < 2) return 0.0;
    const double n = static_cast<double>(count);
    const double mean = sum / n;
    const double var = std::max(0.0, sumsq / n - mean * mean);
    if (var <= 0.0) return 0.0; // every candidate scored the same: the best one is not an outlier
    return (best - mean) / std::sqrt(var);
  }

  void ProSEAlgorithm::annotatePeptDeepFeatures_(const PeakMap& spectra,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    if (!peptdeep_enable_) { return; }
#ifdef WITH_ONNX
    startProgress(0, 1, "Adding PeptDeep rescoring features...");
    PeptDeepRescoring rescoring;
    Param p = rescoring.getParameters();
    // A bare name resolves against share/OpenMS, an absolute path is returned unchanged.
    // The hint is needed because the compiled-in data path is the *source* share/OpenMS, which has
    // no models/ -- they are downloaded into the build tree. On an installed OpenMS the hint is a
    // dead path and File::find falls through to the installed data path, where 'make install' puts
    // them. This is the same resolution PeptDeepLibraryPredictor uses, deliberately: not relative
    // to the executable, which would only hold for a tool sitting in bin/.
    const StringList model_dirs = {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"};
    p.setValue("ms2_model", File::find(peptdeep_ms2_model_, model_dirs));
    p.setValue("rt_model", File::find(peptdeep_rt_model_, model_dirs));
    p.setValue("instrument", peptdeep_instrument_);
    p.setValue("nce", peptdeep_nce_);
    p.setValue("rt_model_type", peptdeep_rt_model_type_);
    p.setValue("batch_size", peptdeep_batch_size_);
    // Inference otherwise runs at its own default while the search around it scales with
    // the thread count the user actually asked for.
#ifdef _OPENMP
    p.setValue("threads", std::max(1, omp_get_max_threads()));
#else
    p.setValue("threads", 1);
#endif
    rescoring.setParameters(p);
    rescoring.setLogType(getLogType());
    rescoring.annotate(spectra, protein_ids, peptide_ids);
    endProgress();
    // 'peptdeep:enable' is an explicit request: an input that ends up without a single
    // predicted PSM must not silently continue with search-engine features only. An input
    // without PSMs has nothing to predict: its (empty) result is reported as it is.
    if (rescoring.getNumPredictedPSMs() == 0 && rescoring.getNumPSMs() > 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "[ProSE] 'peptdeep:enable' is set, but PeptDeep prediction added features to 0 / "
        + std::to_string(rescoring.getNumPSMs()) + " PSMs. Check that 'annotate:PSM' includes fragment annotations "
        "('ALL' or 'fragment_annotation').");
    }
#else
    (void)spectra; (void)protein_ids; (void)peptide_ids;
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "[ProSE] 'peptdeep:enable' is set, but this OpenMS was built without ONNX support "
      "(WITH_ONNX=OFF), so the PeptDeep prediction-based rescoring features cannot be computed. "
      "Use an OpenMS build with -DWITH_ONNX=ON or set 'peptdeep:enable' to 'false'.");
#endif
  }

  namespace
  {
    /// Scratch buffers of alignAbsoluteTolerance_(), reused across the hits and spectra of a thread (rewritten by every call).
    struct AlignmentScratch_
    {
      std::vector<double> prev;         ///< scores of the computed cells of row i - 1
      std::vector<double> cur;          ///< scores of the computed cells of row i
      std::vector<unsigned char> dir;   ///< traceback direction of every computed cell, rows concatenated
      std::vector<Size> row_begin;      ///< first computed column of row i
      std::vector<Size> row_offset;     ///< position of row i in dir (entry i + 1 closes row i)
    };

    /**
      @brief Absolute-tolerance (Da) branch of SpectrumAlignment::getSpectrumAlignment on flat storage.

      Returns exactly the alignment of the banded dynamic programme in SpectrumAlignment.h, which keeps
      its cells in std::map<Size, std::map<Size, ...>> (two tree nodes and several lookups per cell).
      Identity: the reference computes, per row i, one contiguous run of columns (from the left border
      it carries along to the column where it leaves the band), in the same order and with the same
      expressions as below. Every other cell it reads is either a border cell (i * tolerance or
      j * tolerance, (0,0) = 0) or absent, in which case it substitutes (i + j) * tolerance: the same
      value, so one formula serves both. A traceback step onto an absent cell reads a default (0,0)
      entry there, which ends the walk.
    */
    void alignAbsoluteTolerance_(std::vector<std::pair<Size, Size>>& alignment,
                                 const MSSpectrum& s1, const MSSpectrum& s2,
                                 const double tolerance, AlignmentScratch_& scratch)
    {
      if (!s1.isSorted() || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
      }
      alignment.clear();

      const Size n1 = s1.size();
      const Size n2 = s2.size();
      // cell (i, j) outside the computed runs: border value or the reference's substitute for an absent cell
      const auto outside = [tolerance](const Size i_plus_j) { return i_plus_j == 0 ? 0.0 : i_plus_j * tolerance; };
      enum : unsigned char { DIAGONAL = 1, UP = 2, LEFT = 3 }; // predecessor (i-1,j-1), (i,j-1), (i-1,j)

      std::vector<double>& prev = scratch.prev;
      std::vector<double>& cur = scratch.cur;
      std::vector<unsigned char>& dir = scratch.dir;
      prev.clear();
      dir.clear();
      scratch.row_begin.assign(n1 + 2, 0);
      scratch.row_offset.assign(n1 + 2, 0);

      Size left_ptr(1);
      Size last_i(0), last_j(0);
      Size prev_begin(1), prev_end(0); // computed columns of row i - 1 (row 0: none)
      for (Size i = 1; i <= n1; ++i)
      {
        const double pos1(s1[i - 1].getMZ());
        const Size cur_begin = left_ptr;
        cur.clear();
        scratch.row_begin[i] = cur_begin;
        scratch.row_offset[i] = dir.size();

        for (Size j = cur_begin; j <= n2; ++j)
        {
          bool off_band(false);
          const double pos2(s2[j - 1].getMZ());
          const double diff_align = fabs(pos1 - pos2);

          // running off the right border of the band?
          if (pos2 > pos1 && diff_align > tolerance)
          {
            if (i < n1 && j < n2 && s1[i].getMZ() < pos2)
            {
              off_band = true;
            }
          }

          // can we tighten the left border of the band?
          if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
          {
            ++left_ptr;
          }

          double score_align = diff_align;
          if (j - 1 >= prev_begin && j - 1 <= prev_end)
          {
            score_align += prev[j - 1 - prev_begin];
          }
          else
          {
            score_align += outside(i - 1 + j - 1);
          }

          double score_up = tolerance;
          if (j > cur_begin)
          {
            score_up += cur.back();
          }
          else
          {
            score_up += outside(i + j - 1);
          }

          double score_left = tolerance;
          if (j >= prev_begin && j <= prev_end)
          {
            score_left += prev[j - prev_begin];
          }
          else
          {
            score_left += outside(i - 1 + j);
          }

          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
            cur.push_back(score_align);
            dir.push_back(DIAGONAL);
            last_i = i;
            last_j = j;
          }
          else if (score_up <= score_left)
          {
            cur.push_back(score_up);
            dir.push_back(UP);
          }
          else
          {
            cur.push_back(score_left);
            dir.push_back(LEFT);
          }

          if (off_band)
          {
            break;
          }
        }
        prev_begin = cur_begin;
        prev_end = cur_begin + cur.size() - 1;
        prev.swap(cur);
      }
      scratch.row_offset[n1 + 1] = dir.size();

      // do traceback
      Size i = last_i;
      Size j = last_j;
      while (i >= 1 && j >= 1)
      {
        const Size begin = scratch.row_begin[i];
        const Size width = scratch.row_offset[i + 1] - scratch.row_offset[i];
        if (j < begin || j >= begin + width)
        {
          // absent cell: the reference reads (0,0), which counts as a diagonal step only at (1,1)
          if (i == 1 && j == 1) alignment.emplace_back(0, 0);
          break;
        }
        const unsigned char d = dir[scratch.row_offset[i] + (j - begin)];
        if (d == DIAGONAL)
        {
          alignment.emplace_back(i - 1, j - 1);
          --i;
          --j;
        }
        else if (d == UP)
        {
          --j;
        }
        else
        {
          --i;
        }
      }
      std::reverse(alignment.begin(), alignment.end());
    }
  }

  void ProSEAlgorithm::postProcessHits_(const PeakMap& exp,
        std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
        const std::vector<CandidatePoolStats_>& pool_stats,
        std::vector<ProteinIdentification>& protein_ids,
        PeptideIdentificationList& peptide_ids,
        Size top_hits,
  //      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
  //      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
  //      Size max_variable_mods_per_peptide, TODO: what about this parameter?
        const StringList& modifications_fixed,
        const StringList& modifications_variable,
        Int peptide_missed_cleavages,
        double precursor_mass_tolerance,
        double fragment_mass_tolerance,
        const std::string& precursor_mass_tolerance_unit_ppm,
        const std::string& fragment_mass_tolerance_unit_ppm,
        const Int precursor_min_charge,
        const Int precursor_max_charge,
        const std::string& enzyme,
        const std::string& database_name,
        const PeakMap* evidence_spectra,
        const std::vector<Size>* evidence_parent,
        bool defer_peak_annotations) const
  {
    // Candidate-pool features (delta score, z-score, candidate count) are derived
    // from @p pool_stats rather than from @p annotated_hits: scoreSpectraAgainstIndex_
    // already pruned each spectrum to max(top_hits, 2) candidates, so the hits left
    // here are the report set, not the search space. pool_stats summarises every
    // candidate that was scored, zero-scoring ones included, and is accumulated
    // across chunks by the chunked search paths.
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }

    // One value per spectrum, kept in side vectors so scoring does not have to carry
    // them on every AnnotatedHit_ candidate.
    std::vector<float> delta_scores(annotated_hits.size(), 0.0f);
    std::vector<float> hyperscore_zscores(annotated_hits.size(), 0.0f);
    std::vector<float> ln_num_candidates(annotated_hits.size(), 0.0f);
#pragma omp parallel for default(none) shared(annotated_hits, top_hits, pool_stats, delta_scores, hyperscore_zscores, ln_num_candidates)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      const CandidatePoolStats_& stats = pool_stats[scan_index];
      delta_scores[scan_index] = static_cast<float>(stats.deltaScore());
      hyperscore_zscores[scan_index] = static_cast<float>(stats.zScore());
      // log1p rather than log: a spectrum can end up with zero candidates, and
      // ln(0) is not representable. The +1 offset is order-preserving, so this
      // stays a monotone proxy for the size of the scored search space.
      ln_num_candidates[scan_index] = static_cast<float>(std::log1p(static_cast<double>(stats.count)));

      // sort and keep n best elements according to score (unchanged)
      auto& hits = annotated_hits[scan_index];
      Size topn = top_hits > hits.size() ? hits.size() : top_hits;
      std::partial_sort(hits.begin(), hits.begin() + topn, hits.end(), AnnotatedHit_::hasBetterScore);
      hits.resize(topn);
      hits.shrink_to_fit();
    }

    bool annotation_precursor_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_fragment_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_prefix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_suffix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_num_matched_peaks = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::NUM_MATCHED_PEAKS) != annotate_psm_.end();
    bool annotation_matched_prefix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS) != annotate_psm_.end();
    bool annotation_matched_suffix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS) != annotate_psm_.end();
    bool annotation_longest_ion_run = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE) != annotate_psm_.end();
    bool annotation_matched_ion_current = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT) != annotate_psm_.end();
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end();
    bool annotation_hyperscore_zscore = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::HYPERSCORE_ZSCORE) != annotate_psm_.end();
    bool annotation_ln_num_candidates = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LN_NUM_CANDIDATES) != annotate_psm_.end();
    bool annotation_matched_ion_current_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT_FRACTION) != annotate_psm_.end();
    bool annotation_complementary_ions_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::COMPLEMENTARY_IONS_FRACTION) != annotate_psm_.end();
    const bool annotation_local_evidence = param_.getValue("annotate:local_fragment_evidence").toBool();
    const bool evidence_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    const Size n_evidence = evidence_spectra != nullptr ? evidence_spectra->size() : 0;
    if (annotation_local_evidence && (evidence_spectra == nullptr
        || n_evidence + (evidence_parent != nullptr ? evidence_parent->size() : 0) != exp.size()))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Local fragment evidence requires the aligned peak lists retained during preprocessing.");
    }

    // "ALL" adds all annotations
    if (std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end())
    {
      annotation_precursor_error_ppm = true;
      annotation_fragment_error_ppm = true;
      annotation_prefix_fraction = true;
      annotation_suffix_fraction = true;
      annotation_num_matched_peaks = true;
      annotation_matched_prefix_ions = true;
      annotation_matched_suffix_ions = true;
      annotation_longest_ion_run = true;
      annotation_matched_ion_current = true;
      annotation_fragment_annotations = true;
      annotation_hyperscore_zscore = true;
      annotation_ln_num_candidates = true;
      annotation_matched_ion_current_fraction = true;
      annotation_complementary_ions_fraction = true;
    }

    // Alignment is needed for fragment error, fragment annotations, longest ion run, MIC,
    // normalized MIC, and complementary ion pairs
    // (deferred fragment annotations are aligned when they are made)
    const bool need_alignment = annotation_fragment_error_ppm || (annotation_fragment_annotations && !defer_peak_annotations)
      || annotation_longest_ion_run || annotation_matched_ion_current || annotation_matched_ion_current_fraction
      || annotation_complementary_ions_fraction;
    // no ion name is written here: the theoretical spectra may come without them
    const bool without_ion_names = !annotation_fragment_annotations || defer_peak_annotations;

    // Both configurations depend only on the fragment tolerance, so they are
    // shared read-only by every thread of the loop below: getSpectrum() and
    // getSpectrumAlignment() are const and hold no mutable state (the scoring
    // loop in scoreSpectraAgainstIndex_ shares its generators the same way).
    // The alignment tolerance mirrors the search's fragment tolerance so the
    // reported FRAGMENT_ERROR_MEDIAN_PPM is not polluted by far-off spurious
    // matches — SpectrumAlignment's default is 0.3 Da absolute, which is ~30×
    // looser than a typical 20 ppm search.
    // The generators mirror the ion series of the scoring, so annotation features (fragment
    // error, peak annotations, longest ion run, MIC) are computed on the same theoretical
    // spectrum the candidate was scored against.
    const SpectrumGenerators_ generators = spectrumGenerators_();
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }
    // Da tolerance: same alignment from a flat-storage copy of SpectrumAlignment's banded DP
    // (alignAbsoluteTolerance_ above); ppm tolerance keeps SpectrumAlignment's cheap matching.
    const bool sa_absolute = !sa.getParameters().getValue("is_relative_tolerance").toBool();
    const double sa_tolerance = (double)sa.getParameters().getValue("tolerance");

    // Meta value keys of the loop below, resolved to registry indices once per call: every
    // string-keyed setMetaValue()/getMetaValue() takes the registry's process-wide lock (also for
    // names that are already registered), which serialised the annotation threads.
    // Identity: MetaInfo keeps (and idXML writes) the values of a hit ordered by registry index, so
    // new names must get their indices in the order in which the loop registered them: the first
    // spectrum with hits registers "scan_index", then IM if it carries a drift time, then its first
    // hit registers the per-hit names in the order they are set below. All hits set the same names,
    // so no thread could register a later name before an earlier one and this order did not depend
    // on the thread count. registerName() only looks up names that are already known, and without
    // any hit nothing was registered, so nothing is registered here either. Names that stay
    // string-keyed (spectrum reference, IM, rank; at most one lookup per spectrum or hit) are
    // still registered inside the loop, after the ones below, as before (checked: the registry
    // contents after a search are the same as before this change, at 1 and at 16 threads).
    const int isotope_error_sign = isotope_error_observed_minus_theoretical_ ? -1 : 1;
    const bool open_search_mode = isOpenSearchMode_();
    UInt mv_scan_index{}, mv_fragment_error{}, mv_precursor_error{}, mv_prefix_fraction{}, mv_suffix_fraction{},
         mv_num_matched_peaks{}, mv_matched_prefix_ions{}, mv_matched_suffix_ions{}, mv_delta_score{},
         mv_hyperscore_zscore{}, mv_ln_num_candidates{}, mv_matched_ion_current{}, mv_matched_ion_current_fraction{},
         mv_longest_ion_run{}, mv_complementary_ions_fraction{}, mv_chance_match_surprise{},
         mv_mass_competition_evidence{}, mv_isotope_error{}, mv_delta_mass{}, mv_rank{};
    {
      const auto first_with_hits = std::find_if(annotated_hits.begin(), annotated_hits.end(),
                                                [](const std::vector<AnnotatedHit_>& hits) { return !hits.empty(); });
      if (first_with_hits != annotated_hits.end())
      {
        MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
        mv_scan_index = registry.registerName("scan_index");
        const MSSpectrum& first_spectrum = exp[first_with_hits - annotated_hits.begin()];
        if (IMTypes::determineIMFormat(first_spectrum) == IMFormat::IM_SPECTRUM)
        {
          registry.registerName(first_spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE
                                  ? Constants::UserParam::FAIMS_CV : Constants::UserParam::IM);
        }
        if (annotation_fragment_error_ppm) mv_fragment_error = registry.registerName(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
        if (annotation_precursor_error_ppm) mv_precursor_error = registry.registerName(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
        if (annotation_prefix_fraction) mv_prefix_fraction = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
        if (annotation_suffix_fraction) mv_suffix_fraction = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
        if (annotation_num_matched_peaks) mv_num_matched_peaks = registry.registerName(Constants::UserParam::NUM_MATCHED_PEAKS);
        if (annotation_matched_prefix_ions) mv_matched_prefix_ions = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS);
        if (annotation_matched_suffix_ions) mv_matched_suffix_ions = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS);
        mv_delta_score = registry.registerName(Constants::UserParam::DELTA_SCORE);
        if (annotation_hyperscore_zscore) mv_hyperscore_zscore = registry.registerName(Constants::UserParam::HYPERSCORE_ZSCORE);
        if (annotation_ln_num_candidates) mv_ln_num_candidates = registry.registerName(Constants::UserParam::LN_NUM_CANDIDATES);
        if (annotation_matched_ion_current) mv_matched_ion_current = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT);
        if (annotation_matched_ion_current_fraction) mv_matched_ion_current_fraction = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
        if (annotation_longest_ion_run) mv_longest_ion_run = registry.registerName(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
        if (annotation_complementary_ions_fraction) mv_complementary_ions_fraction = registry.registerName(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
        if (annotation_local_evidence)
        {
          mv_chance_match_surprise = registry.registerName(Constants::UserParam::CHANCE_MATCH_SURPRISE);
          mv_mass_competition_evidence = registry.registerName(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
        }
        mv_isotope_error = registry.registerName(Constants::UserParam::ISOTOPE_ERROR);
        if (open_search_mode) mv_delta_mass = registry.registerName("DeltaMass");
        // The rank of the hits below the top one (PeptideHit::setRank() by name: one lookup under the registry's
        // lock per hit), registered where the loop below registered it (by its first second-best hit), after the names
        // above: no other name of a hit is registered in between.
        if (std::any_of(annotated_hits.begin(), annotated_hits.end(), [](const std::vector<AnnotatedHit_>& hits) { return hits.size() > 1; }))
        {
          mv_rank = registry.registerName(Constants::UserParam::RANK);
        }
      }
    }

    // The meta values a hit gets below, and later in the search (protein mapping: 2, ion priors: 2): reserved at once.
    const Size hit_meta_values = Size(annotation_fragment_error_ppm) + annotation_precursor_error_ppm + annotation_prefix_fraction
      + annotation_suffix_fraction + annotation_num_matched_peaks + annotation_matched_prefix_ions + annotation_matched_suffix_ions
      + 1 /* delta score */ + annotation_hyperscore_zscore + annotation_ln_num_candidates + annotation_matched_ion_current
      + annotation_matched_ion_current_fraction + annotation_longest_ion_run + annotation_complementary_ions_fraction
      + 2 * Size(annotation_local_evidence) + 1 /* isotope error */ + open_search_mode + 1 /* rank */ + 4;

    // With an empty list to start from, every spectrum with hits gets a slot in scan order: the identifications
    // end up in scan order without a lock and without the sort below (which orders them by scan_index, i.e. the
    // same order). Otherwise they are appended under a lock and sorted, as before.
    const bool slots = peptide_ids.empty();
    std::vector<Size> slot_of_scan;
    if (slots)
    {
      slot_of_scan.resize(annotated_hits.size());
      Size used_slots = 0;
      for (Size scan_index = 0; scan_index < annotated_hits.size(); ++scan_index)
      {
        slot_of_scan[scan_index] = used_slots;
        if (!annotated_hits[scan_index].empty()) ++used_slots;
      }
      peptide_ids.resize(used_slots);
    }

    // The work per spectrum varies with its number of hits: hand out small blocks
#pragma omp parallel
    {
    // scratch of the theoretical spectra, alignments and local evidence, reused by all hits of all spectra of a thread
    // (as before by all hits of a spectrum: each use overwrites what it reads)
    AlignmentScratch_ alignment_scratch;
    std::vector<double> evidence_alternatives;
    // the theoretical spectrum without ion names, see below
    MSSpectrum coded_spec;
    std::vector<uint8_t> coded_series;
    std::vector<UInt> coded_ordinals;
    std::vector<Int> coded_charges;
#pragma omp for schedule(dynamic, 16)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spec = exp[scan_index];
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setSpectrumReference( spec.getNativeID());
        pi.setMetaValue(mv_scan_index, static_cast<unsigned int>(scan_index));
        pi.setScoreType(SEARCH_ENGINE_SCORE); // also for scoring:method mass_accuracy, a log-space HyperScore (recorded in the search parameters)
        pi.setHigherScoreBetter(true);
        double mz = spec.getPrecursors()[0].getMZ();
        pi.setRT(spec.getRT());
        pi.setMZ(mz);

        // Annotate the spectrum's single drift time: a FAIMS compensation voltage as FAIMS_CV (which the
        // FAIMS-aware tools read), an ion mobility (e.g. DDA-PASEF, 1/K0) as IM
        if (IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM)
        {
          pi.setMetaValue(spec.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE
                            ? Constants::UserParam::FAIMS_CV : Constants::UserParam::IM,
                          spec.getDriftTime());
        }

        Size charge = spec.getPrecursors()[0].getCharge();

        // Spectrum-level quantity, identical for every candidate of this spectrum, so it is
        // computed once here rather than per hit.
        const double spectrum_tic =
          annotation_matched_ion_current_fraction ? spec.calculateTIC() : 0.0;
        const MSSpectrum& evidence_spec = !annotation_local_evidence ? spec
          : (*evidence_spectra)[static_cast<Size>(scan_index) < n_evidence ? static_cast<Size>(scan_index)
                                                                           : (*evidence_parent)[scan_index - n_evidence]];
        const std::vector<double> local_densities = annotation_local_evidence
          ? localPeakDensities_(evidence_spec) : std::vector<double>{};

        // create full peptide hit structure from annotated hits (their sequences are moved: the hits are freed after
        // this function)
        vector<PeptideHit> phs;
        phs.reserve(annotated_hits[scan_index].size());
        for (auto& ah : annotated_hits[scan_index])
        {
          PeptideHit ph;
          ph.reserveMetaValues(hit_meta_values); // (one block instead of growing by doubling)
          // Prefer spectrum charge; if absent (0), fall back to the charge actually used by FI for this candidate
          const Size used_charge = (charge > 0) ? charge : static_cast<Size>(ah.applied_charge);
          ph.setCharge(used_charge);
          ph.setScore(ah.score);
          ph.setSequence(std::move(ah.sequence));
          const AASequence& sequence = ph.getSequence();

          // Generate theoretical spectrum + alignment for annotations that need it. Unless the ion names are written
          // (fragment annotations here), the peaks of getSpectrum() come without them, with the series, ordinal and
          // charge of each ion as numbers (getIonSeriesSpectrum(): the same peaks in the same order; about a fifth of
          // the cost), if the generator supports that.
          std::vector<std::pair<Size, Size>> alignment;
          MSSpectrum theoretical_spec;
          bool coded = false;
          // Annotate the charges actually scored when higher charges are scored.
          const int max_frag_z = scoring_multiple_charges_ ? scoringMaxCharge_(static_cast<int>(used_charge))
                                                           : ((charge >= 2) ? std::min<int>(charge - 1, 2) : 1);
          if (need_alignment)
          {
            coded = without_ion_names
                    && tsg.getIonSeriesSpectrum(coded_spec, coded_series, coded_ordinals, coded_charges, sequence, 1, max_frag_z);
            if (!coded) tsg.getSpectrum(theoretical_spec, sequence, 1, max_frag_z);
            if (sa_absolute)
            {
              alignAbsoluteTolerance_(alignment, coded ? coded_spec : theoretical_spec, spec, sa_tolerance, alignment_scratch);
            }
            else
            {
              sa.getSpectrumAlignment(alignment, coded ? coded_spec : theoretical_spec, spec);
            }
          }
          const MSSpectrum& theo_spec = coded ? coded_spec : theoretical_spec;

          if (annotation_fragment_error_ppm)
          {
            std::vector<double> err;
            for (const auto& match : alignment)
            {
              double fragment_error = fabs(Math::getPPM(spec[match.second].getMZ(), theo_spec[match.first].getMZ()));
              err.push_back(fragment_error);
            }
            double median_ppm_error(0);
            if (!err.empty()) { median_ppm_error = Math::median(err.begin(), err.end(), false); }
            ph.setMetaValue(mv_fragment_error, median_ppm_error);
          }

          if (annotation_precursor_error_ppm)
          {
            // Subtract out the isotope offset FI matched at — FragmentIndex searches
            // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs
            // + isotope_error * C13C12, and the observed-to-monoiso correction in m/z is
            //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
            // (ah.isotope_error is this search offset, whatever sign the PSM reports).
            // Without this, a ±1 Da FI match reports ~1000 ppm / charge for the Percolator
            // feature, corrupting target/decoy discrimination.
            const double corrected_mz = monoisotopicPrecursorMZ(mz, ah.isotope_error, static_cast<double>(used_charge));
            double theo_mz = sequence.getMZ(used_charge);
            double ppm_difference = Math::getPPM(corrected_mz, theo_mz);
            ph.setMetaValue(mv_precursor_error, ppm_difference);
          }

          if (annotation_prefix_fraction)
          {
            ph.setMetaValue(mv_prefix_fraction, ah.prefix_fraction);
          }

          if (annotation_suffix_fraction)
          {
            ph.setMetaValue(mv_suffix_fraction, ah.suffix_fraction);
          }

          // Matched ion counts (from scoring, no alignment needed)
          if (annotation_num_matched_peaks)
          {
            ph.setMetaValue(mv_num_matched_peaks, static_cast<int>(ah.matched_prefix_ions + ah.matched_suffix_ions));
          }
          if (annotation_matched_prefix_ions)
          {
            ph.setMetaValue(mv_matched_prefix_ions, static_cast<int>(ah.matched_prefix_ions));
          }
          if (annotation_matched_suffix_ions)
          {
            ph.setMetaValue(mv_matched_suffix_ions, static_cast<int>(ah.matched_suffix_ions));
          }

          ph.setMetaValue(mv_delta_score, delta_scores[scan_index]);

          if (annotation_hyperscore_zscore)
          {
            ph.setMetaValue(mv_hyperscore_zscore, hyperscore_zscores[scan_index]);
          }
          if (annotation_ln_num_candidates)
          {
            ph.setMetaValue(mv_ln_num_candidates, ln_num_candidates[scan_index]);
          }

          // Fragment annotations, longest ion run, MIC, normalized MIC, and complementary
          // ion pairs all iterate the alignment + ion names
          if ((annotation_fragment_annotations && !defer_peak_annotations) || annotation_longest_ion_run || annotation_matched_ion_current
            || annotation_matched_ion_current_fraction || annotation_complementary_ions_fraction)
          {
            // (the names and charges of getSpectrum(); none for a coded spectrum)
            const auto* ion_names = theoretical_spec.getStringDataArrays().empty() ? nullptr : &theoretical_spec.getStringDataArrays()[0];
            const auto* ion_charges = theoretical_spec.getIntegerDataArrays().empty() ? nullptr : &theoretical_spec.getIntegerDataArrays()[0];

            // Build PeakAnnotation vector + collect ion ordinals for longest run.
            // Prefix = a/b/c (N-terminal), suffix = x/y/z (C-terminal). Ordinals
            // for different ion types at the same cleavage position (e.g. a3 + b3)
            // are merged via std::unique below — each ordinal is a backbone
            // position, not an ion-type-specific identifier.
            std::vector<PeptideHit::PeakAnnotation> peak_annotations;
            std::vector<int> prefix_ordinals, suffix_ordinals;
            double matched_ion_current = 0.0;
            const bool need_mic = annotation_matched_ion_current || annotation_matched_ion_current_fraction;
            const bool need_ordinals = annotation_longest_ion_run || annotation_complementary_ions_fraction;
            // Dedup guard for MIC: in ppm-alignment mode a single experimental
            // peak can match multiple theoretical peaks (e.g. b-ion and near
            // isotope), so we must sum each exp_idx at most once. Sized only
            // when MIC is actually requested.
            std::vector<char> counted_exp_peaks(need_mic ? spec.size() : 0, 0);
            if (annotation_fragment_annotations && !defer_peak_annotations) { peak_annotations.reserve(alignment.size()); }

            for (const auto& [theo_idx, exp_idx] : alignment)
            {
              if (annotation_fragment_annotations && !defer_peak_annotations)
              {
                PeptideHit::PeakAnnotation pa;
                pa.mz = spec[exp_idx].getMZ();
                pa.intensity = spec[exp_idx].getIntensity();
                pa.annotation = (*ion_names)[theo_idx];
                pa.charge = (*ion_charges)[theo_idx];
                peak_annotations.push_back(pa);
              }

              if (need_mic && !counted_exp_peaks[exp_idx])
              {
                matched_ion_current += spec[exp_idx].getIntensity();
                counted_exp_peaks[exp_idx] = 1;
              }

              if (need_ordinals && coded)
              {
                // the number of the ion name ("b5+", "z.4+": 5, 4), never 0
                (coded_series[theo_idx] == 0 ? prefix_ordinals : suffix_ordinals).push_back(static_cast<int>(coded_ordinals[theo_idx]));
              }
              else if (need_ordinals && (*ion_names)[theo_idx].size() >= 2)
              {
                const std::string& name = (*ion_names)[theo_idx];
                const char c = name[0];
                const bool is_prefix = (c == 'a' || c == 'b' || c == 'c');
                const bool is_suffix = (c == 'x' || c == 'y' || c == 'z');
                if (is_prefix || is_suffix)
                {
                  // Extract ordinal: "b5", "y3-H2O1+", "c12++", "z.4+" (z+1) -> 5, 3, 12, 4
                  Size pos = 1;
                  while (pos < name.size() && (name[pos] == '.' || name[pos] == '\'')) ++pos; // z. (z+1), z' (z+2)
                  const Size ordinal_begin = pos;
                  while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
                  if (pos > ordinal_begin)
                  {
                    int ordinal = StringUtils::toInt32(StringUtils::substr(name, ordinal_begin, pos - ordinal_begin));
                    (is_prefix ? prefix_ordinals : suffix_ordinals).push_back(ordinal);
                  }
                }
              }
            }

            if (annotation_fragment_annotations && !defer_peak_annotations)
            {
              ph.setPeakAnnotations(std::move(peak_annotations));
            }

            if (annotation_matched_ion_current)
            {
              ph.setMetaValue(mv_matched_ion_current, matched_ion_current);
            }

            if (annotation_matched_ion_current_fraction)
            {
              ph.setMetaValue(mv_matched_ion_current_fraction,
                              spectrum_tic > 0 ? matched_ion_current / spectrum_tic : 0.0);
            }

            if (need_ordinals)
            {
              // Compute longest consecutive run across prefix and suffix series.
              // Sorts + deduplicates each ordinal vector in place (a backbone position
              // matched by multiple ion types, e.g. a3 and b3, counts once).
              auto longestRun = [](std::vector<int>& v) -> int {
                if (v.empty()) return 0;
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                int best = 1, run = 1;
                for (Size i = 1; i < v.size(); ++i)
                {
                  if (v[i] == v[i - 1] + 1) { ++run; if (run > best) best = run; }
                  else run = 1;
                }
                return best;
              };
              int longest_prefix = longestRun(prefix_ordinals);
              int longest_suffix = longestRun(suffix_ordinals);

              if (annotation_longest_ion_run)
              {
                ph.setMetaValue(mv_longest_ion_run, std::max(longest_prefix, longest_suffix));
              }

              if (annotation_complementary_ions_fraction)
              {
                // A prefix ion at backbone position i (a_i/b_i/c_i) is complementary to a
                // suffix ion at position (peptide_length - i) (x/y/z at the same cleavage
                // site). Andromeda-style structural signal, distinct from the independently
                // computed prefix/suffix fractions above.
                const int pep_len = static_cast<int>(sequence.size());
                double complementary_fraction = 0.0;
                if (pep_len > 1)
                {
                  Size n_complementary = 0;
                  for (int p : prefix_ordinals)
                  {
                    if (std::binary_search(suffix_ordinals.begin(), suffix_ordinals.end(), pep_len - p)) { ++n_complementary; }
                  }
                  complementary_fraction = static_cast<double>(n_complementary) / static_cast<double>(pep_len - 1);
                }
                ph.setMetaValue(mv_complementary_ions_fraction, complementary_fraction);
              }
            }
          }


          if (annotation_local_evidence)
          {
            // Use the fragment charges that were scored. The annotation spectrum (charges 1..max_frag_z)
            // contains them whenever it was generated, so it is reused instead of generating a second one.
            const int evidence_z = scoringMaxCharge_(static_cast<int>(used_charge));
            LocalFragmentEvidence_ evidence;
            if (need_alignment && evidence_z <= max_frag_z)
            {
              evidence = coded ? localFragmentEvidence_(evidence_spec, coded_spec, coded_charges, evidence_z, local_densities,
                                                        fragment_mass_tolerance, evidence_ppm, evidence_alternatives)
                               : localFragmentEvidence_(evidence_spec, theoretical_spec, evidence_z, local_densities,
                                                        fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            else if (tsg.getIonSeriesSpectrum(coded_spec, coded_series, coded_ordinals, coded_charges, sequence, 1, evidence_z))
            {
              // (the peaks and charges of getSpectrum(), without the names)
              evidence = localFragmentEvidence_(evidence_spec, coded_spec, coded_charges, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            else
            {
              MSSpectrum evidence_theory;
              tsg.getSpectrum(evidence_theory, sequence, 1, evidence_z);
              evidence = localFragmentEvidence_(evidence_spec, evidence_theory, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            ph.setMetaValue(mv_chance_match_surprise, evidence.chance_match_surprise);
            ph.setMetaValue(mv_mass_competition_evidence, evidence.mass_competition_evidence);
          }

          // Add isotope error metavalue (always; exposed as Percolator feature). ah.isotope_error is the offset
          // FragmentIndex added to the observed mass (theoretical minus observed); report:isotope_error_convention
          // selects the reported sign.
          ph.setMetaValue(mv_isotope_error, isotope_error_sign * ah.isotope_error);

          // Add delta mass metavalue for open search
          if (open_search_mode)
          {
            ph.setMetaValue(mv_delta_mass, ah.delta_mass);
          }

          // store PSM
          phs.push_back(std::move(ph));
        }
        pi.setHits(std::move(phs));
        // Ensure hits are sorted by score (best first), then assign ranks explicitly (0 = top hit)
        pi.sort();
        {
          std::vector<PeptideHit>& hits = pi.getHits();
          for (Size r = 0; r < hits.size(); ++r)
          {
            // as setRank(r): the rank of the top hit (0) is the default, which a new hit does not carry
            if (r != 0) hits[r].setMetaValue(mv_rank, static_cast<UInt>(r));
          }
        }

        // Debug: log spectrum-level top hit details before storing PeptideIdentification.
        // DEBUG (not INFO) because multi-file mode would emit one line per scan per input.
        if (!pi.getHits().empty())
        {
          const PeptideHit& top_hit = pi.getHits().front();
          OPENMS_LOG_DEBUG << "[ProSE] scan_index=" << scan_index
                           << " top_ln(hyperscore)=" << top_hit.getScore()
                           << " top_charge=" << top_hit.getCharge()
                           << " top_isotope_error=" << (int)top_hit.getMetaValue(mv_isotope_error)
                           << std::endl;
        }
        if (slots)
        {
          peptide_ids[slot_of_scan[scan_index]] = std::move(pi);
        }
        else
        {
#pragma omp critical (peptide_ids_access)
          {
            //clang-tidy: seems to be a false-positive in combination with omp
            peptide_ids.push_back(std::move(pi));
          }
        }
      }
    }
    } // omp parallel

#ifdef _OPENMP
    // we need to sort the peptide_ids by scan_index in order to have the same output in the idXML-file
    if (omp_get_max_threads() > 1 && !slots)
    {
      // one registry lookup for the whole sort instead of two (locked) lookups per comparison;
      // getMetaValue(name) is getMetaValue(getIndex(name)), so the comparisons are unchanged
      const UInt scan_index_key = MetaInfoInterface::metaRegistry().getIndex("scan_index");
      std::sort(peptide_ids.begin(), peptide_ids.end(), [scan_index_key](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return a.getMetaValue(scan_index_key) < b.getMetaValue(scan_index_key);
      });
    }
#endif

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("ProSE");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    DateTime now = DateTime::now();
    std::string identifier("ProSE_" + now.get());
    protein_ids[0].setIdentifier(identifier);
    for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges =StringUtils::toStr(precursor_min_charge) + ":" + StringUtils::toStr(precursor_max_charge);

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = modifications_fixed;
    search_parameters.variable_modifications = modifications_variable;
    search_parameters.missed_cleavages = peptide_missed_cleavages;
    search_parameters.fragment_mass_tolerance = fragment_mass_tolerance;
    search_parameters.precursor_mass_tolerance = precursor_mass_tolerance;
    search_parameters.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);

    // add additional percolator features or post-processing
    StringList feature_set{"score"};
    if (annotation_fragment_error_ppm) feature_set.push_back(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
    if (annotation_prefix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
    if (annotation_suffix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
    if (annotation_longest_ion_run) feature_set.push_back(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
    if (annotation_matched_prefix_ions) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS);
    if (annotation_matched_suffix_ions) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS);
    if (annotation_matched_ion_current) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT);
    if (annotation_matched_ion_current_fraction) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
    if (annotation_complementary_ions_fraction) feature_set.push_back(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
    if (annotation_hyperscore_zscore) feature_set.push_back(Constants::UserParam::HYPERSCORE_ZSCORE);
    if (annotation_ln_num_candidates) feature_set.push_back(Constants::UserParam::LN_NUM_CANDIDATES);
    if (annotation_local_evidence)
    {
      feature_set.push_back(Constants::UserParam::CHANCE_MATCH_SURPRISE);
      feature_set.push_back(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
    }
    feature_set.push_back(Constants::UserParam::DELTA_SCORE);
    feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
    if (chimericSearch_()) feature_set.push_back(CHIMERIC_RANK); // set by search() after this
    // note: precursor error is calculated by percolator itself
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    // Readers tell the sign of the PSMs' isotope_error by this record; without it (earlier ProSE versions, or
    // report:isotope_error_convention=theoretical_minus_observed) the sign is theoretical minus observed.
    if (isotope_error_observed_minus_theoretical_)
    {
      search_parameters.setMetaValue("isotope_error_convention", "observed_minus_theoretical");
    }
    // record whether open-search mode was used
    search_parameters.setMetaValue("open_search", isOpenSearchMode_() ? "true" : "false");

    search_parameters.setMetaValue("peptide:clip_nterm_methionine", param_.getValue("peptide:clip_nterm_methionine"));
    search_parameters.setMetaValue("peptide:deduplicate", param_.getValue("peptide:deduplicate"));
    // recorded only when not the default (the default output stays as before)
    if (decoy_keep_termini_) search_parameters.setMetaValue("decoy_method", "keep_termini");
    search_parameters.setMetaValue("fragment:query_spectrum", param_.getValue("fragment:query_spectrum"));
    search_parameters.setMetaValue("fragment:query_spectrum_resolved", query_raw_spectrum_ ? "raw" : "processed");
    search_parameters.setMetaValue("annotate:local_fragment_evidence", param_.getValue("annotate:local_fragment_evidence"));
    search_parameters.setMetaValue("scoring:fragment_charges", param_.getValue("scoring:fragment_charges"));
    search_parameters.setMetaValue("scoring:fragment_charges_resolved", scoring_multiple_charges_ ? "multiple" : "single");
    if (mass_accuracy_score_) // recorded when it changes the score, so HyperScore searches keep their parameter list
    {
      search_parameters.setMetaValue("scoring:method", param_.getValue("scoring:method"));
      search_parameters.setMetaValue("scoring:method_resolved", "mass_accuracy");
      search_parameters.setMetaValue("scoring:mass_error_sd", mass_error_sd_ppm_);
    }
    search_parameters.setMetaValue("peaks:window_type", peaks_window_type_);
    search_parameters.setMetaValue(
      "peaks:window_type_resolved",
      peaks_window_type_ == "auto"
        ? (Deisotoper::isToleranceSupported(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm") ? "jump_full" : "jump")
        : peaks_window_type_);
    // MS2 deisotoping rule (used only when the spectra are deisotoped, see fragment:deisotope)
    search_parameters.setMetaValue("fragment:deisotope_min_peaks", static_cast<int>(deisotoping_.min_peaks));
    search_parameters.setMetaValue("fragment:deisotope_charge_cap", deisotoping_.charge_cap_precursor ? "precursor" : "none");
    search_parameters.setMetaValue("fragment:deisotope_sum_intensity", deisotoping_.sum_intensity ? "true" : "false");

    search_parameters.enzyme_term_specificity = peptide_enzyme_specificity_;
    protein_ids[0].setSearchParameters(std::move(search_parameters));

    // Annotate the IM unit on the ProteinIdentification if all PeptideIdentifications have an IM value. The
    // unit is that of the spectra the values came from (FAIMS compensation voltages are written as FAIMS_CV,
    // so they never count here)
    if (!peptide_ids.empty())
    {
      bool all_have_im = std::all_of(
        peptide_ids.begin(), peptide_ids.end(),
        [](const PeptideIdentification& pid) { return pid.metaValueExists(Constants::UserParam::IM); });

      if (all_have_im)
      {
        const auto with_im = std::find_if(exp.begin(), exp.end(), [](const MSSpectrum& spec) {
          return IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM
                 && spec.getDriftTimeUnit() != DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
        });
        protein_ids[0].setMetaValue(Constants::UserParam::IM,
          (with_im != exp.end() ? *with_im : exp[0]).getDriftTimeUnitAsString());
      }
    }
  }

  // =====================================================================
  // Build a SearchContext (decoys + FragmentIndex) from a FASTA database.
  // Hoisted out of the in-memory search() body so callers can build the
  // index once and reuse it across many spectrum files.
  // =====================================================================
  Param ProSEAlgorithm::fragmentIndexParameters_(bool electron_ions) const
  {
    Param p = getParameters();
    // FragmentIndex has its own boolean 'decoys' flag {true,false}; ProSE's enum
    // value (auto/generate/ignore) would fail its validation. ProSE builds the
    // decoy-augmented database itself, so the FragmentIndex must never generate.
    p.remove("decoys");
    p.setValue("decoys", "false");
    p.remove("decoy_method"); // ProSE's own (the decoys of the database it builds)
    p.remove("decoy_parity_report");
    // ions:by_activation is resolved by the caller, which knows the spectra. The c and z+1 ions
    // for electron-activated spectra go into a set of their own, which the other spectra are not
    // matched against (see scoreSpectraAgainstIndex_).
    p.remove("ions:by_activation");
    p.setValue("ions:electron_ions", electron_ions ? "true" : "false");
    // annotation-only parameters of ProSE (see annotateIonPriors_)
    p.remove("annotate:self_trained_ion_priors");
    p.removeAll("annotate:ion_prior_");
    p.remove("peptide:protein_mapping"); // ProSE's own (FragmentIndex::getProteinOccurrences() serves it)
    return p;
  }

  namespace
  {
    /**
      @brief DecoyHelper::findDecoyString(@p db, quiet = true), with the OpenMP threads for large databases.

      DecoyHelper::countDecoys() sums, over the proteins, the matches of the prefix and the suffix pattern and keeps,
      per match, the spelling of its last occurrence (within a protein the suffix after the prefix). Chunks of
      consecutive proteins counted on their own and joined in their order give the same. The decision is the one of
      findDecoyString() for these statistics; it does not depend on the order in which it looks at the matches,
      since at most one prefix and one suffix can reach 80% of all prefixes or suffixes.

      At most 8 threads: the step (about 60 ms serially for 50,000 entries) gains little beyond 4 threads, and a
      larger team costs CPU time (measured at 64 threads: 2-4 CPU-s per run for no shorter wall time).

      Most accessions cannot match a pattern, and searching for the suffix pattern tries every position: an accession
      is searched only if mayMatchDecoyPattern() says it can match.
    */
    DecoyHelper::Result findDecoyStringInParallel(const std::vector<FASTAFile::FASTAEntry>& db)
    {
      // Whether DecoyHelper's prefix (or suffix) pattern can match the lower-case accession: "^(decoy_*|dec_*|...)"
      // matches exactly the accessions that start with an affix, "(_decoy*|_dec*|...)$" exactly those that end with
      // '_', an affix without its last character and this character any number of times. Boost's ^ and $ also match
      // at the line breaks inside a text: such an accession is always searched.
      const auto mayMatchDecoyPattern = [](const std::string& lower, bool suffix)
      {
        if (lower.find_first_of("\n\r\f") != std::string::npos) { return true; }
        for (const std::string& affix : DecoyHelper::affixes)
        {
          if (!suffix)
          {
            if (StringUtils::hasPrefix(lower, affix)) { return true; }
            continue;
          }
          const std::string_view stem(affix.data(), affix.size() - 1);
          for (std::string_view head(lower); ; head.remove_suffix(1))
          {
            if (head.size() > stem.size() && head.substr(head.size() - stem.size()) == stem && head[head.size() - stem.size() - 1] == '_')
            {
              return true;
            }
            if (head.empty() || head.back() != affix.back()) { break; }
          }
        }
        return false;
      };
#ifdef _OPENMP
      constexpr int max_threads = 8;
      const int threads = omp_in_parallel() ? 1 : std::min(omp_get_max_threads(), max_threads);
#else
      const int threads = 1;
#endif
      if (threads < 2 || db.size() < 8192)
      {
        FASTAContainer<TFI_Vector> container(db);
        return DecoyHelper::findDecoyString(container, /*quiet=*/true);
      }
      struct Counts
      {
        std::map<std::string, std::pair<Size, Size>> count; ///< match -> occurrences as prefix, as suffix
        std::map<std::string, std::string> spelling;        ///< match -> spelling of its last occurrence
        Size prefixes = 0, suffixes = 0;
      };
      const SignedSize n = static_cast<SignedSize>(db.size());
      const SignedSize num_chunks = 4 * static_cast<SignedSize>(threads);
      std::vector<Counts> chunks(num_chunks);
#pragma omp parallel num_threads(threads)
      {
        const RegularExpression prefix_pattern(DecoyHelper::regexstr_prefix);
        const RegularExpression suffix_pattern(DecoyHelper::regexstr_suffix);
        std::string lower, match;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize c = 0; c < num_chunks; ++c)
        {
          Counts& chunk = chunks[c];
          for (SignedSize i = n * c / num_chunks; i < n * (c + 1) / num_chunks; ++i)
          {
            const std::string& identifier = db[i].identifier;
            lower = identifier;
            StringUtils::toLower(lower);
            if (mayMatchDecoyPattern(lower, false) && prefix_pattern.search(lower, &match))
            {
              ++chunk.prefixes;
              ++chunk.count[match].first;
              chunk.spelling[match] = StringUtils::prefix(identifier, match.length());
            }
            if (mayMatchDecoyPattern(lower, true) && suffix_pattern.search(lower, &match))
            {
              ++chunk.suffixes;
              ++chunk.count[match].second;
              chunk.spelling[match] = StringUtils::suffix(identifier, match.length());
            }
          }
        }
      }
      Counts all;
      for (const Counts& chunk : chunks) // in order: later chunks hold later proteins
      {
        all.prefixes += chunk.prefixes;
        all.suffixes += chunk.suffixes;
        for (const auto& [match, count] : chunk.count)
        {
          all.count[match].first += count.first;
          all.count[match].second += count.second;
        }
        for (const auto& [match, spelling] : chunk.spelling) all.spelling[match] = spelling;
      }
      // DecoyHelper::findDecoyString()'s decision
      const double proteins = static_cast<double>(db.size());
      if (static_cast<double>(all.prefixes + all.suffixes) < 0.4 * proteins || all.prefixes == all.suffixes)
      {
        return {false, "?", true};
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.first) / static_cast<double>(all.prefixes) >= 0.8 && static_cast<double>(count.first) / proteins >= 0.4)
        {
          return {true, all.spelling[match], true};
        }
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.second) / static_cast<double>(all.suffixes) >= 0.8 && static_cast<double>(count.second) / proteins >= 0.4)
        {
          return {true, all.spelling[match], false};
        }
      }
      return {false, "?", true};
    }
  }

  // static
  Size ProSEAlgorithm::removeDuplicateProteins(std::vector<FASTAFile::FASTAEntry>& db, const std::string& source)
  {
    // accession -> position of its first entry; the entries are compared by their sequences
    std::unordered_map<std::string_view, Size> first;
    first.reserve(db.size());
    std::vector<bool> repeated(db.size(), false);
    Size n_repeated = 0;
    for (Size i = 0; i < db.size(); ++i)
    {
      const auto [it, inserted] = first.emplace(db[i].identifier, i);
      if (inserted) { continue; }
      if (db[it->second].sequence != db[i].sequence)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source,
          "The protein database '" + source + "' has two entries of the accession '" + db[i].identifier + "' with "
          "different sequences (entries " + StringUtils::toStr(it->second + 1) + " and " + StringUtils::toStr(i + 1)
          + "). An accession must name one protein: rename or remove one of them.");
      }
      repeated[i] = true;
      ++n_repeated;
    }
    if (n_repeated == 0) { return 0; }
    first.clear(); // (views of the identifiers, which are moved below)
    Size kept = 0;
    for (Size i = 0; i < db.size(); ++i)
    {
      if (repeated[i]) { continue; }
      if (kept != i) { db[kept] = std::move(db[i]); }
      ++kept;
    }
    db.resize(kept);
    OPENMS_LOG_WARN << "[ProSE] The protein database '" << source << "' repeats " << n_repeated << " entr"
                    << (n_repeated == 1 ? "y" : "ies") << " (same accession and sequence as an entry before, e.g. the "
                    << "file or the contaminants were added twice); searching each protein once." << std::endl;
    return n_repeated;
  }

  ProSEAlgorithm::DecoyStrategy_
  ProSEAlgorithm::resolveDecoyStrategy_(const std::vector<FASTAFile::FASTAEntry>& db) const
  {
    // Detect pre-existing decoys. First the common-marker heuristic
    // (DecoyHelper: "decoy", "rev", "xxx", ... as prefix or suffix), then a
    // literal fall-back to the configured decoy_prefix so custom markers
    // outside the vocabulary are still recognised.
    // quiet=true: a target-only database is a normal case here (auto/generate
    // then synthesise decoys), so suppress DecoyHelper's "unable to determine
    // decoy string" ERROR/WARN noise — we handle the negative result ourselves.
    const DecoyHelper::Result det = findDecoyStringInParallel(db);

    bool existing = det.success;
    std::string ext_string = det.success ? det.name : decoy_prefix_;
    bool ext_is_prefix = det.success ? det.is_prefix : true;

    Size n_decoys = 0; // the entries with the marker (as PeptideIndexing finds them: this spelling exactly)
    if (!existing && !decoy_prefix_.empty())
    {
      Size n_prefix = 0, n_suffix = 0;
      for (const FASTAFile::FASTAEntry& e : db)
      {
        n_prefix += StringUtils::hasPrefix(e.identifier, decoy_prefix_);
        n_suffix += StringUtils::hasSuffix(e.identifier, decoy_prefix_);
      }
      if (n_prefix + n_suffix > 0)
      {
        existing = true;
        ext_string = decoy_prefix_;
        ext_is_prefix = n_prefix > 0; // prefer prefix when both orientations occur
        n_decoys = ext_is_prefix ? n_prefix : n_suffix;
      }
    }
    else if (existing)
    {
      n_decoys = countDecoyEntries_(db, ext_string, ext_is_prefix);
    }

    DecoyStrategy_ s;
    switch (decoy_mode_)
    {
      case DecoyMode_::AUTO:
        if (existing)
        {
          // Reuse the decoys already in the database; never double-generate. Target-decoy FDR estimates the false
          // targets by the decoys, one decoy per target: with far fewer (e.g. a few entries of the decoy_prefix found
          // literally, which DecoyHelper's 40 % rule does not detect) or far more decoys, every q-value is off by
          // about that factor. Refuse such a database rather than reporting wrong FDRs.
          const Size n_targets = db.size() - n_decoys;
          if (static_cast<double>(n_decoys) < MIN_DECOYS_PER_TARGET * static_cast<double>(n_targets)
              || static_cast<double>(n_decoys) > MAX_DECOYS_PER_TARGET * static_cast<double>(n_targets))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "decoys=auto: the protein database has " + StringUtils::toStr(n_decoys) + " decoy protein(s) (marker '"
              + ext_string + "', " + (ext_is_prefix ? "prefix" : "suffix") + ") for " + StringUtils::toStr(n_targets)
              + " target protein(s). Target-decoy FDR needs about one decoy per target (ProSE reuses decoys at "
              + StringUtils::toStr(MIN_DECOYS_PER_TARGET) + " to " + StringUtils::toStr(MAX_DECOYS_PER_TARGET)
              + " per target); with these its q-values would be wrong. Set 'decoys' to 'generate' to replace them "
              "with decoys generated from the targets, or use a database with one decoy per target.");
          }
          s.generate = false; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = ext_string; s.is_prefix = ext_is_prefix;
        }
        else
        {
          // No decoys present: synthesise them so target-decoy FDR is possible.
          s.generate = true; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = decoy_prefix_; s.is_prefix = true;
        }
        break;

      case DecoyMode_::GENERATE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=generate: removing pre-existing decoys (marker '"
                          << ext_string << "') and regenerating from the target proteins.\n";
        }
        else
        {
          // Nothing detected to strip. If the database still carries some decoy-like accessions
          // with a non-standard or too-sparse marker (below DecoyHelper's threshold, and not the
          // configured decoy_prefix), they are treated as targets and reversed into decoys-of-decoys.
          // Warn so the user can set -Search:decoy_prefix to the real marker or pre-clean the FASTA.
          FASTAContainer<TFI_Vector> stats_container(db);
          const DecoyHelper::DecoyStatistics stats = DecoyHelper::countDecoys(stats_container);
          if (stats.all_prefix_occur + stats.all_suffix_occur > 0)
          {
            OPENMS_LOG_WARN << "[ProSE] decoys=generate: " << (stats.all_prefix_occur + stats.all_suffix_occur)
                            << " accession(s) carry a decoy-like marker but too few to auto-detect, so "
                            << "they are not stripped and will be treated as targets (risking "
                            << "decoys-of-decoys). Set -Search:decoy_prefix to the actual marker or "
                            << "pre-clean the database.\n";
          }
        }
        s.generate = true; s.strip_existing = existing; s.have_decoys = true;
        s.decoy_string = decoy_prefix_; s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;

      case DecoyMode_::IGNORE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=ignore: removing pre-existing decoys (marker '"
                          << ext_string << "'); searching the target proteins only.\n";
        }
        s.generate = false; s.strip_existing = existing; s.have_decoys = false;
        s.decoy_string.clear(); s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;
    }
    return s;
  }

  namespace
  {
    // The first exception that a parallel loop caught (lowest index; others are dropped), rethrown after the loop: an
    // exception must not leave an OpenMP region
    struct FirstError
    {
      std::exception_ptr error;
      SignedSize index = std::numeric_limits<SignedSize>::max();
      void set(SignedSize i)
      {
#pragma omp critical (ProSEAlgorithm_firstError)
        if (i < index)
        {
          index = i;
          error = std::current_exception();
        }
      }
      void rethrow() const
      {
        if (error) std::rethrow_exception(error);
      }
    };

    // A cleavage rule of the form "(?<=[KR])" with an optional "(?!P)" (Trypsin, Trypsin/P, Lys-C, Arg-C, ...; as
    // FragmentIndex's regex-free digest): a site after a residue of the first set unless the next one blocks it
    struct SegmentRule
    {
      std::array<bool, 256> after{};
      std::array<bool, 256> not_before{};
    };

    // True (and @p rule filled) if @p regex is exactly "(?<=[" + letters A-Z + "])", optionally followed by "(?!" + one
    // letter + ")"; every other rule is left to the regex
    bool parseSegmentRule(const std::string& regex, SegmentRule& rule)
    {
      const auto is_letter = [](char c) { return c >= 'A' && c <= 'Z'; };
      if (regex.compare(0, 5, "(?<=[") != 0) return false;
      Size i = 5;
      for (; i < regex.size() && is_letter(regex[i]); ++i) rule.after[static_cast<unsigned char>(regex[i])] = true;
      if (i == 5 || regex.compare(i, 2, "])") != 0) return false;
      i += 2;
      if (i == regex.size()) return true;
      if (regex.size() != i + 5 || regex.compare(i, 3, "(?!") != 0 || !is_letter(regex[i + 3]) || regex[i + 4] != ')') return false;
      rule.not_before[static_cast<unsigned char>(regex[i + 3])] = true;
      return true;
    }

    // The segments {start, length} of @p sequence: digestor.digestUnmodified(sequence, spans, 1, 0) (fully specific, no
    // missed cleavages), with @p rule without the regex (~1 ms per protein with the regex): its tokenize_() yields 0 and
    // every site 0 < p < size, and every segment between consecutive ones (and the last one up to the end) is a product.
    void digestSegments(const ProteaseDigestion& digestor, const SegmentRule* rule, const std::string& sequence,
                        std::vector<std::pair<Size, Size>>& spans)
    {
      if (rule == nullptr)
      {
        digestor.digestUnmodified(sequence, spans, 1, 0);
        return;
      }
      spans.clear();
      if (sequence.empty()) return;
      Size start = 0;
      for (Size p = 1; p < sequence.size(); ++p)
      {
        if (rule->after[static_cast<unsigned char>(sequence[p - 1])] && !rule->not_before[static_cast<unsigned char>(sequence[p])])
        {
          spans.emplace_back(start, p - start);
          start = p;
        }
      }
      spans.emplace_back(start, sequence.size() - start);
    }

    // ProteaseDigestion::isValidProduct(protein, pos, length, true, true, false) of a fully specific digestion with the
    // enzyme of @p rule, without the regex: the start is a site (or the protein start, or within the first three residues
    // of a protein that starts with M), and the end is a site or the protein end. (isValidProduct_() tokenizes
    // [pos - 2, end + 2): every site p inside it with pos - 2 < p < end + 2 sees the residues p - 1 and p.)
    // Requires 0 <= pos, 0 < length, pos + length <= protein.size().
    bool isValidFullySpecific(const SegmentRule& rule, const std::string& protein, Size pos, Size length)
    {
      const auto site = [&rule, &protein](Size p)
      {
        return rule.after[static_cast<unsigned char>(protein[p - 1])] && !rule.not_before[static_cast<unsigned char>(protein[p])];
      };
      const Size end = pos + length;
      const bool spec_n = pos == 0 || (pos <= 2 && protein[0] == 'M') || site(pos);
      const bool spec_c = end == protein.size() || site(end);
      return spec_n && spec_c;
    }

    // decoy_method 'keep_termini': OpenMS DecoyDatabase's construction (method reverse, Decoy:keepPeptideNTerm
    // and keepPeptideCTerm true; Decoy:non_shuffle_pattern "KRP", its default is empty) with @p digestor (the search enzyme, no missed cleavages,
    // full specificity): every segment is MRMDecoy::reversePeptideLight(segment, keepN, keepC, "KRP"). Residues not in
    // any segment (none for a digest that tiles the protein) stay in place.
    std::string keepTerminiDecoy(const std::string& protein, const ProteaseDigestion& digestor, const SegmentRule* rule,
                                 std::vector<std::pair<Size, Size>>& spans)
    {
      digestSegments(digestor, rule, protein, spans);
      std::string decoy = protein;
      for (const auto& [start, length] : spans)
      {
        decoy.replace(start, length, MRMDecoy::reversePeptideLight(protein.substr(start, length), {}, true, true, "KRP").first);
      }
      return decoy;
    }

    // FNV-1a 64 of the generated decoy sequences, each followed by '\n', in target order (decoy_parity_report)
    void fnv1a(UInt64& h, const std::string& s)
    {
      for (const unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
      h ^= static_cast<unsigned char>('\n'); h *= 1099511628211ULL;
    }
  }

  // Build the searched database according to @p strategy: optionally strip
  // pre-existing decoys, optionally generate fresh decoys from the targets.
  // Produces FASTA entries only — does not build a FragmentIndex.
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    return buildDecoyAugmentedDB_(std::vector<FASTAFile::FASTAEntry>(fasta_db), strategy);
  }

  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    // The entries are filtered in place instead of being copied one by one: the same entries
    // in the same order.
    std::vector<FASTAFile::FASTAEntry> db = std::move(fasta_db);
    // For in-memory databases; one read from a file had its repeats removed by loadFASTA() already (and its conflicts
    // reported with the file name), so this finds none.
    removeDuplicateProteins(db, "the protein database");

    // decoys=auto reusing pre-existing decoys logs nothing otherwise; surface the auto-detected
    // marker so a rare DecoyHelper misdetection (a target DB whose accessions start with a decoy
    // affix) is diagnosable. Single emission: buildDecoyAugmentedDB_ runs once per search.
    const bool reuses_decoys = decoy_mode_ == DecoyMode_::AUTO && !strategy.generate && strategy.have_decoys;
    if (reuses_decoys)
    {
      OPENMS_LOG_INFO << "[ProSE] decoys=auto: reusing existing decoys detected in the database "
                      << "(marker '" << strategy.decoy_string << "', "
                      << (strategy.is_prefix ? "prefix" : "suffix") << ")." << std::endl;
    }
    // Initial-Met clipping adds N-terminal peptides of every protein that starts with M. Generated decoys keep the
    // initial Met (below); supplied decoys often do not (a reversed protein ends with it), and then the clipped
    // peptides enlarge the target space only. Count both kinds while the entries are filtered.
    const bool check_met_symmetry = reuses_decoys && param_.getValue("peptide:clip_nterm_methionine").toBool();
    Size targets = 0, targets_with_met = 0, decoys = 0, decoys_with_met = 0;

    // 1. Keep targets, dropping pre-existing decoys when requested. A stop codon ('*') that ends
    //    a sequence, as in databases translated from genomes (e.g. SGD), is not a residue: remove
    //    it so the C-terminal peptide stays searchable and decoys are built from the protein
    //    alone. FragmentIndex skips peptides that contain a stop codon inside the sequence.
    //    An entry left without residues has nothing to search, and decoy generation needs
    //    residues: drop it. Lowercase residues are read as uppercase ones, as the FASTA format
    //    defines them (NCBI); they are neither indexed nor reversed otherwise.
    const auto existing_decoy = [&strategy](const FASTAFile::FASTAEntry& e)
    {
      return strategy.strip_existing &&
        (strategy.strip_is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.strip_string)
                                  : StringUtils::hasSuffix(e.identifier, strategy.strip_string));
    };
    // (the lowercase check reads every residue: in parallel, 10-20 ms on one thread for a human target-decoy database)
    Size lowercase_entries = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+: lowercase_entries) if (db.size() >= 1024)
    for (SignedSize i = 0; i < static_cast<SignedSize>(db.size()); ++i)
    {
      std::string& sequence = db[i].sequence;
      if (!existing_decoy(db[i]) && std::any_of(sequence.begin(), sequence.end(), [](char c) { return c >= 'a' && c <= 'z'; }))
      {
        StringUtils::toUpper(sequence);
        ++lowercase_entries;
      }
    }
    auto kept = db.begin();
    for (FASTAFile::FASTAEntry& e : db)
    {
      const bool is_existing_decoy = existing_decoy(e);
      if (!is_existing_decoy)
      {
        std::string& sequence = e.sequence;
        while (!sequence.empty() && sequence.back() == '*') { sequence.pop_back(); }
        if (!sequence.empty())
        {
          if (check_met_symmetry)
          {
            const bool is_decoy = strategy.is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.decoy_string)
                                                     : StringUtils::hasSuffix(e.identifier, strategy.decoy_string);
            (is_decoy ? decoys : targets) += 1;
            (is_decoy ? decoys_with_met : targets_with_met) += (sequence.size() > 1 && sequence[0] == 'M') ? 1 : 0;
          }
          if (&*kept != &e) { *kept = std::move(e); }
          ++kept;
        }
      }
    }
    db.erase(kept, db.end());
    if (lowercase_entries > 0)
    {
      OPENMS_LOG_INFO << "[ProSE] " << lowercase_entries << " database entries hold lowercase residues; they are searched "
                      << "as uppercase ones." << std::endl;
    }
    // Warn when the share of decoys that start with M is below half the share of targets that do.
    if (check_met_symmetry && targets_with_met > 0 && 2 * decoys_with_met * targets < targets_with_met * decoys)
    {
      OPENMS_LOG_WARN << "[ProSE] peptide:clip_nterm_methionine: " << targets_with_met << " of " << targets
                      << " target proteins but only " << decoys_with_met << " of " << decoys << " decoy proteins "
                      << "in the database start with M. Removing the initial Met adds N-terminal peptides to the "
                      << "targets that have no decoy counterpart, which makes target-decoy FDR estimates slightly "
                      << "optimistic. Use '-Search:decoys generate' (generated decoys keep the initial Met) or "
                      << "'-Search:peptide:clip_nterm_methionine false' for a symmetric search space." << std::endl;
    }
    // The options that rely on symmetric decoys, with decoys taken from the database: say so (non-default options only)
    if (reuses_decoys && decoy_keep_termini_)
    {
      OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini has no effect: the decoys of the database are used "
                      << "(decoys 'auto'). Use '-Search:decoys generate' to replace them." << std::endl;
    }
    if (check_met_symmetry && (terminal_features_ || stratum_q_values_)
        && targets_with_met > 0 && 2 * decoys_with_met * targets < targets_with_met * decoys)
    {
      OPENMS_LOG_WARN << "[ProSE] rescoring:terminal_features / rescoring:stratum_q_values with these decoys: the "
                      << "Met-clipped (and protein N-terminal) peptides have (almost) no decoys, so the features tell "
                      << "targets from decoys and the stratum q-values are too low; the FDR is not controlled. Use '-Search:decoys generate "
                      << "-Search:decoy_method keep_termini'." << std::endl;
    }

    // 2. Generate decoys by reversing the (remaining) target proteins.
    if (strategy.generate)
    {
      // decoy_string is the prefix prepended to generated accessions; an empty prefix would
      // produce decoys with the same accession as their targets, silently breaking FDR.
      if (strategy.decoy_string.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "decoy_prefix must be non-empty to generate decoys (decoys=auto without existing decoys, "
          "or decoys=generate).");
      }
      DecoyGenerator decoy_generator;
      const size_t old_size = db.size();
      db.reserve(old_size * 2);
      // decoy_method 'keep_termini' (fully specific digestion only: the candidates are unions of whole segments; otherwise
      // 'reverse')
      const bool keep_termini = decoy_keep_termini_ && peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_FULL
                                && enzyme_ != EnzymaticDigestion::UnspecificCleavage && enzyme_ != EnzymaticDigestion::NoCleavage;
      if (decoy_keep_termini_ && !keep_termini)
      {
        OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini needs a fully specific digestion with a cleaving enzyme; "
                        << "the decoys are generated with 'reverse'." << std::endl;
      }
      ProteaseDigestion segments;
      Size site_mismatches = 0;
      UInt64 checksum = 14695981039346656037ULL;
      SegmentRule segment_rule;
      bool simple_segments = false;
      if (keep_termini)
      {
        segments.setEnzyme(enzyme_);
        segments.setMissedCleavages(0);
        segments.setSpecificity(EnzymaticDigestion::SPEC_FULL);
        simple_segments = parseSegmentRule(ProteaseDB::getInstance()->getEnzyme(enzyme_)->getRegEx(), segment_rule);
      }
      const SegmentRule* const segment_rule_ptr = simple_segments ? &segment_rule : nullptr;
      // The decoy sequences in parallel (each from its target alone; serially 1.2-1.6 s per run for a human database, as
      // AASequence parsing and the digestion per protein dominate), then appended in target order as before.
      const bool clip_nterm_methionine = param_.getValue("peptide:clip_nterm_methionine").toBool();
      std::vector<std::string> decoy_sequences(old_size);
      std::vector<char> site_mismatch(old_size, 0);
      std::vector<std::exception_ptr> decoy_errors(old_size); // an exception must not leave the parallel region
      // Only sequences of plain residue letters in parallel: AASequence::fromString() of a sequence with modification
      // syntax may register modifications (shared state) and warn; then serially, stopping at the first error, as before.
      bool plain_letters = true;
#pragma omp parallel for schedule(dynamic, 256) reduction(&&: plain_letters) if (old_size >= 256)
      for (SignedSize si = 0; si < static_cast<SignedSize>(old_size); ++si)
      {
        const std::string& sequence = db[static_cast<size_t>(si)].sequence;
        plain_letters = plain_letters && std::all_of(sequence.begin(), sequence.end(), [](char c) { return c >= 'A' && c <= 'Z'; });
      }
      const bool decoys_in_parallel = plain_letters && old_size >= 256;
      bool serial_stop = false; // (written only without the parallel region)
#pragma omp parallel if (decoys_in_parallel)
      {
        std::vector<std::pair<Size, Size>> target_spans, decoy_spans;
        // monotonic: without the parallel region (serial fallback) the proteins are processed in order
#pragma omp for schedule(monotonic: dynamic, 64)
        for (SignedSize si = 0; si < static_cast<SignedSize>(old_size); ++si)
        {
          const size_t i = static_cast<size_t>(si);
          if (serial_stop) { continue; }
          try
          {
            const std::string& sequence = db[i].sequence;
            if (keep_termini)
            {
              decoy_sequences[i] = keepTerminiDecoy(sequence, segments, segment_rule_ptr, target_spans);
              // the guarantee holds only if the decoy has the target's cleavage sites (enzymes whose rule reads K, R, P)
              digestSegments(segments, segment_rule_ptr, decoy_sequences[i], decoy_spans);
              site_mismatch[i] = decoy_spans != target_spans;
              continue;
            }
            // Keep an initial Met on generated decoys as well: otherwise clipping
            // expands the target N-terminal search space without its decoy counterpart.
            const bool preserve_met = clip_nterm_methionine && sequence.size() > 1 && sequence[0] == 'M';
            const std::string target = preserve_met ? sequence.substr(1) : sequence;
            std::string decoy = peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE
              ? decoy_generator.reverseProtein(AASequence::fromString(target)).toString()
              : decoy_generator.reversePeptides(AASequence::fromString(target), enzyme_).toString();
            if (preserve_met) { decoy.insert(decoy.begin(), 'M'); }
            decoy_sequences[i] = std::move(decoy);
          }
          catch (...)
          {
            decoy_errors[i] = std::current_exception();
            if (!decoys_in_parallel) { serial_stop = true; } // (serially: the first error stops, as before)
          }
        }
      }
      for (const std::exception_ptr& error : decoy_errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = db[i];
        e.sequence = std::move(decoy_sequences[i]);
        site_mismatches += site_mismatch[i];
        if (keep_termini || decoy_parity_report_) fnv1a(checksum, e.sequence);
        e.identifier = strategy.decoy_string + e.identifier;  // decoy_string is the prefix to add
        db.push_back(std::move(e));
      }
      std::vector<std::string>().swap(decoy_sequences);
      if (keep_termini)
      {
        if (site_mismatches > 0)
        {
          OPENMS_LOG_WARN << "[ProSE] decoy_method keep_termini: " << site_mismatches << " of " << old_size << " decoy "
                          << "proteins do not have the cleavage sites of their target (enzyme '" << enzyme_ << "'); their "
                          << "peptides are not paired with the target's." << std::endl;
        }
      }
      if (decoy_parity_report_)
      {
        OPENMS_LOG_INFO << "[ProSE] decoy parity: " << old_size << " decoys generated (" << (keep_termini ? "keep_termini" : "reverse")
                        << "), cleavage-site mismatches " << (keep_termini ? StringUtils::toStr(site_mismatches) : std::string("not checked"))
                        << ", FNV-1a 64 of the decoy sequences in target order " << std::hex << checksum << std::dec << std::endl;
      }
      Math::RandomShuffler shuffler(42);  // fixed seed for reproducible decoy ordering across runs/files
      shuffler.portable_random_shuffle(db.begin(), db.end());

      // Under the default decoys=auto, decoys are generated whenever the database lacks them —
      // including no-ProSE-FDR runs, because the decoys are typically consumed by a downstream or
      // external FDR step (FalseDiscoveryRate / Percolator) or a merged run. This roughly doubles
      // index build + search; for a genuinely target-only search, use decoys=ignore.
      if (decoy_mode_ == DecoyMode_::AUTO && fdr_psm_ == 0.0 && fdr_protein_ == 0.0)
      {
        OPENMS_LOG_INFO << "[ProSE] decoys=auto generated " << old_size << " decoys (the database had "
                        << "none) for downstream target-decoy FDR. Use '-Search:decoys ignore' for a "
                        << "target-only search if you do not need FDR." << std::endl;
      }
    }
    return db;
  }

  // =====================================================================
  // Strided calibration sample. Used by the chunked search paths so that
  // calibration always sees a representative, suitably-sized subset of the
  // full DB — independent of database_chunk_size_. Prevents the first-chunk
  // starvation mode where small chunk_size values (e.g. 500) produce a
  // calibration pool below calibration_min_psms_ and silently disable
  // calibration. See issue #9182.
  // =====================================================================
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildCalibrationSample_(
      const std::vector<FASTAFile::FASTAEntry>& full_db) const
  {
    // Calibration FI size is bounded by the user's memory budget — which is
    // exactly what they declared via database_chunk_size_. This keeps peak RSS
    // during calibration identical to peak RSS during the main search
    // (one chunk FI at a time), regardless of digestion mode.
    //
    // Critical for immunopeptidomics (non-specific digestion, MHC-I/II variant
    // databases): a naive fixed-size 5000-protein cal sample can generate
    // 50+ GB of fragment index because each protein produces thousands of
    // candidate peptides — defeating the memory savings chunking is designed
    // to deliver. Tying cal_target to chunk_size keeps the promise that
    // "user's memory budget = one chunk's FI."
    //
    // Strided (not contiguous / random) to give the calibration pool even
    // coverage of the full DB.
    //
    // Residual gap vs non-chunked calibration:
    //   - Trypsin / moderate chunk_size (e.g. 5000): cal sample is adequate;
    //     yield is within ~2% of non-chunked.
    //   - Small chunk_size (e.g. 100 proteins for MHC-II variants): cal pool
    //     can fall below calibration_min_psms_ → runCalibrationPass_ returns
    //     success=false, calibration silently no-ops. The main search uses
    //     user-configured tolerances.
    // Follow-up issue: pool calibration PSMs across the first K chunks of the
    // main search loop, so calibration quality becomes independent of
    // chunk_size.
    const Size cal_target = std::min(
        std::max<Size>(database_chunk_size_, Size(1)),
        full_db.size());
    if (cal_target >= full_db.size()) return full_db;

    const Size stride = std::max<Size>(1, full_db.size() / cal_target);
    std::vector<FASTAFile::FASTAEntry> cal_db;
    cal_db.reserve(cal_target);
    for (Size i = 0; i < full_db.size() && cal_db.size() < cal_target; i += stride)
    {
      cal_db.push_back(full_db[i]);
    }
    return cal_db;
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db) const
  {
    return prepareContext(fasta_db, false);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db, bool electron_ions) const
  {
    return prepareContext_(std::vector<FASTAFile::FASTAEntry>(fasta_db), electron_ions);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db, bool electron_ions,
      const std::function<const PeakMap*(Size)>& searched_spectra) const
  {
    SearchContext ctx;

    startProgress(0, 1, "Generate decoys...");
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    ctx.db = buildDecoyAugmentedDB_(std::move(fasta_db), strategy);
    ctx.decoy_string = strategy.decoy_string;
    ctx.decoy_is_prefix = strategy.is_prefix;
    ctx.have_decoys = strategy.have_decoys;
    endProgress();

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    Param this_params = fragmentIndexParameters_(electron_ions);
    ctx.fragment_index.setParameters(this_params);
    ctx.fragment_index.build(ctx.db, searched_spectra);
    ctx.electron_ions = electron_ions;
    endProgress();
    reportDecoyParity_(ctx);

    return ctx;
  }

  void ProSEAlgorithm::reportDecoyParity_(const SearchContext& ctx) const
  {
    if (!decoy_parity_report_) return;
    const std::vector<FASTAFile::FASTAEntry>& db = ctx.db;
    const bool clip = param_.getValue("peptide:clip_nterm_methionine").toBool();
    std::vector<char> decoy(db.size(), 0);
    if (!ctx.decoy_string.empty())
    {
      for (Size i = 0; i < db.size(); ++i)
      {
        decoy[i] = ctx.decoy_is_prefix ? StringUtils::hasPrefix(db[i].identifier, ctx.decoy_string)
                                       : StringUtils::hasSuffix(db[i].identifier, ctx.decoy_string);
      }
    }
    // subclass of an occurrence: 0 internal, 1 protein N-term starting with M, 2 protein N-term otherwise, 3 Met-clipped
    const auto subclass = [&](UInt32 p, Size start) -> int
    {
      const std::string& q = db[p].sequence;
      if (start == 0) return (!q.empty() && q[0] == 'M') ? 1 : 2;
      if (clip && start == 1 && q[0] == 'M') return 3;
      return 0;
    };
    Size entries[2][4] = {}, occurrences[2][4] = {}, acetyl[2][4] = {};
    // protein N-terminal acetylation (a variable modification such as Acetyl (Protein N-term)) of a kept entry
    std::vector<signed char> acetylated;
    const auto is_acetylated = [&](Size idx) -> bool
    {
      const FragmentIndex::Peptide& pe = ctx.fragment_index.getPeptides()[idx];
      if (pe.mod_bitmask_ == 0) return false;
      if (acetylated.empty()) acetylated.assign(ctx.fragment_index.getPeptides().size(), -1);
      if (acetylated[idx] < 0)
      {
        const AASequence seq = ctx.fragment_index.reconstructModifiedSequence(pe, db);
        acetylated[idx] = (seq.hasNTerminalModification()
                           && seq.getNTerminalModificationName().find("Acetyl") != std::string::npos) ? 1 : 0;
      }
      return acetylated[idx] == 1;
    };
    std::set<std::string_view> distinct[2][4];
    std::set<std::string_view> any_label[2];
    const std::vector<FragmentIndex::Peptide>& peptides = ctx.fragment_index.getPeptides();
    const auto count = [&](UInt32 p, Size start, Size length, bool unmodified)
    {
      const int l = decoy[p] ? 1 : 0, c = subclass(p, start);
      ++occurrences[l][c];
      if (unmodified)
      {
        const std::string_view seq(db[p].sequence.data() + start, length);
        distinct[l][c].insert(seq);
        any_label[l].insert(seq);
      }
    };
    for (Size idx = 0; idx < peptides.size(); ++idx)
    {
      const FragmentIndex::Peptide& pe = peptides[idx];
      const int l = decoy[pe.protein_idx] ? 1 : 0, c = subclass(pe.protein_idx, pe.sequence_.first);
      ++entries[l][c];
      count(pe.protein_idx, pe.sequence_.first, pe.sequence_.second, pe.mod_bitmask_ == 0);
      if (c != 0 && is_acetylated(idx)) ++acetyl[l][c];
    }
    for (const FragmentIndex::RemovedOccurrence& ro : ctx.fragment_index.getRemovedOccurrences())
    {
      const FragmentIndex::Peptide& kept = peptides[ro.peptide_idx];
      count(ro.protein_idx, ro.start, kept.sequence_.second, kept.mod_bitmask_ == 0);
      const int c = subclass(ro.protein_idx, ro.start);
      if (c != 0 && is_acetylated(ro.peptide_idx)) ++acetyl[decoy[ro.protein_idx] ? 1 : 0][c];
    }
    Size shared = 0;
    for (const std::string_view& seq : any_label[0]) shared += any_label[1].count(seq);
    static const char* names[4] = {"internal", "nterm_M", "nterm_other", "metclip"};
    OPENMS_LOG_INFO << "[ProSE] decoy parity (fragment index; target / decoy): subclass, peptidoform entries kept, "
                    << "protein occurrences of peptidoforms, distinct unmodified sequences" << std::endl;
    for (int c = 0; c < 4; ++c)
    {
      OPENMS_LOG_INFO << "[ProSE] decoy parity " << names[c] << ": entries " << entries[0][c] << " / " << entries[1][c]
                      << ", occurrences " << occurrences[0][c] << " / " << occurrences[1][c]
                      << ", sequences " << distinct[0][c].size() << " / " << distinct[1][c].size()
                      << ", N-terminally acetylated occurrences " << acetyl[0][c] << " / " << acetyl[1][c]
                      << (occurrences[0][c] == occurrences[1][c] && acetyl[0][c] == acetyl[1][c] ? " (occurrences equal)" : " (occurrences DIFFER)") << std::endl;
    }
    OPENMS_LOG_INFO << "[ProSE] decoy parity: distinct unmodified sequences target " << any_label[0].size() << ", decoy "
                    << any_label[1].size() << ", in both " << shared << std::endl;
  }

  // =====================================================================
  // Score all spectra against a single FragmentIndex, appending results to
  // annotated_hits. Used by both chunked and non-chunked search paths.
  // =====================================================================
  void ProSEAlgorithm::annotateFragmentPeaks_(const PeakMap& exp, PeptideIdentificationList& peptide_ids,
                                              double fragment_mass_tolerance,
                                              const std::string& fragment_mass_tolerance_unit_ppm) const
  {
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end()
      || std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end();
    if (!annotation_fragment_annotations) return;
    // the generators and the alignment of postProcessHits_()
    const SpectrumGenerators_ generators = spectrumGenerators_();
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }
    const bool sa_absolute = !sa.getParameters().getValue("is_relative_tolerance").toBool();
    const double sa_tolerance = (double)sa.getParameters().getValue("tolerance");
    const UInt scan_index_key = MetaInfoInterface::metaRegistry().getIndex("scan_index");

    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      try
      {
        PeptideIdentification& pi = peptide_ids[i];
        if (pi.getHits().empty()) continue;
        const Size scan_index = static_cast<Size>(static_cast<long long>(pi.getMetaValue(scan_index_key))); // (set by postProcessHits_())
        const MSSpectrum& spec = exp[scan_index];
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);
        const Size charge = spec.getPrecursors()[0].getCharge();
        AlignmentScratch_ alignment_scratch;
        for (PeptideHit& ph : pi.getHits())
        {
          // as postProcessHits_() for this hit (its charge is the charge used there)
          const Size used_charge = static_cast<Size>(ph.getCharge());
          const int max_frag_z = scoring_multiple_charges_ ? scoringMaxCharge_(static_cast<int>(used_charge))
                                                           : ((charge >= 2) ? std::min<int>(charge - 1, 2) : 1);
          MSSpectrum theoretical_spec;
          std::vector<std::pair<Size, Size>> alignment;
          tsg.getSpectrum(theoretical_spec, ph.getSequence(), 1, max_frag_z);
          if (sa_absolute)
          {
            alignAbsoluteTolerance_(alignment, theoretical_spec, spec, sa_tolerance, alignment_scratch);
          }
          else
          {
            sa.getSpectrumAlignment(alignment, theoretical_spec, spec);
          }
          const auto& ion_names = theoretical_spec.getStringDataArrays()[0];
          const auto& ion_charges = theoretical_spec.getIntegerDataArrays()[0];
          std::vector<PeptideHit::PeakAnnotation> peak_annotations;
          peak_annotations.reserve(alignment.size());
          for (const auto& [theo_idx, exp_idx] : alignment)
          {
            PeptideHit::PeakAnnotation pa;
            pa.mz = spec[exp_idx].getMZ();
            pa.intensity = spec[exp_idx].getIntensity();
            pa.annotation = ion_names[theo_idx];
            pa.charge = ion_charges[theo_idx];
            peak_annotations.push_back(pa);
          }
          ph.setPeakAnnotations(std::move(peak_annotations));
        }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
  }

  namespace
  {
    // The per-spectrum containers of scoreSpectraAgainstIndex_(), one set per thread for the whole run: the candidates
    // of the query and the theoretical spectrum of the scorer. Constructed and destroyed once per spectrum, they
    // returned their buffers to the allocator 8,000 times per thread and took fresh ones; kept per thread and cleared
    // between spectra, they keep the capacity the largest spectrum of the thread needed and allocate nothing else. The
    // query appends to the candidates, so they are cleared before every query. (A function-local thread_local, not one
    // in the loop body, so that the OpenMP data-sharing clauses of the loop do not see it.)
    struct ScoringScratch_
    {
      FragmentIndex::SpectrumMatchesTopN top_sms;
      PeakSpectrum theo_spectrum;
      std::vector<uint8_t> theo_ion_series;
    };
    ScoringScratch_& scoringScratch_()
    {
      thread_local ScoringScratch_ scratch;
      return scratch;
    }
  } // namespace

  void ProSEAlgorithm::scoreSpectraAgainstIndex_(
      const PeakMap& spectra,
      FragmentIndex& fi,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const SpectrumGenerators_& generators,
      double effective_fragment_tol,
      bool fragment_mass_tolerance_unit_ppm,
      bool open_search_mode,
      std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
      std::vector<CandidatePoolStats_>& pool_stats,
      const std::string& progress_label,
      const PeakMap* query_spectra) const
  {
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }
    if (query_spectra != nullptr && query_spectra->size() != spectra.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, query_spectra->size(),
        "Candidate-retrieval spectra must hold one entry per scoring spectrum.");
    }
    startProgress(0, spectra.size(), progress_label);
    // incremented by every thread and read by the master thread for the progress display: atomic for both
    std::atomic<Size> count_spectra{0};
    const double proton_mass_u = Constants::PROTON_MASS_U;
    const Size keep = std::max(searchTopHits_(), Size(2)); // keep ≥2 for delta score
    const bool deduplicate_chunks = param_.getValue("peptide:deduplicate").toBool() && database_chunk_size_ > 0;

    // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
    // leave the parallel region, which would terminate the process. The one of the SMALLEST failing scan index is
    // rethrown after the loop: what the scan-order loop on one thread reports. After a failure only the spectra before
    // it are still scored (to find an earlier failure); every other spectrum is skipped. (Spectra after the failing one
    // that came earlier in the scoring order may have been scored: the results are discarded with the exception.)
    std::exception_ptr scoring_error;
    const SignedSize n_spectra = static_cast<SignedSize>(spectra.size());
    std::atomic<SignedSize> first_failed{n_spectra};
    // The spectra are scored in the order of their precursor mass (scoringOrder_()): neighbouring iterations then probe
    // the same precursor-window slices of the fragment index, which stay in the caches. Every spectrum is scored on its
    // own and its results go to its scan index: the same results in any order.
    const std::vector<SignedSize> order = scoringOrder_(spectra, query_spectra);
    // While a background read lends threads of this thread, in slices, between which this thread takes them back once
    // the read has ended (see BackgroundSpectra). The slices run over positions in the scoring order.
    const SignedSize slice = readLendsThreads() ? std::max<SignedSize>(1, (n_spectra + 15) / 16) : std::max<SignedSize>(1, n_spectra);
    for (SignedSize slice_begin = 0; slice_begin < n_spectra; slice_begin += slice)
    {
    takeBackThreadsOfEndedRead();
    const SignedSize slice_end = std::min(n_spectra, slice_begin + slice);
    // Runs of neighbouring masses per core (dynamic chunks of up to 64 spectra), but at least 8 chunks per thread of the
    // slice so that small inputs and slices stay balanced. In peptide-mass mode (low_memory) a spectrum costs milliseconds
    // (every candidate in its precursor windows is realised and its fragments generated), not microseconds, and the cost
    // grows with the mass, so the last chunks of the ascending-mass order left the other threads idle (at 64 threads the
    // last 2 of 10 s ran on 37 -> 12 cores): one spectrum per chunk then, the dispatch costs nothing against such spectra.
#ifdef _OPENMP
    const int chunk = fi.isPeptideMassMode() ? 1 : static_cast<int>(std::clamp<SignedSize>(
      (slice_end - slice_begin) / (SignedSize(8) * std::max(1, omp_get_max_threads())), 1, 64));
#else
    [[maybe_unused]] const int chunk = 1; // (only the OpenMP directive uses it)
#endif
#pragma omp parallel for schedule(dynamic, chunk) default(none) shared(annotated_hits, pool_stats, query_spectra, count_spectra, fi, generators, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, effective_fragment_tol, keep, deduplicate_chunks, scoring_error, first_failed, order, slice_begin, slice_end, chunk)
    for (SignedSize position = slice_begin; position < slice_end; ++position)
    {
      const SignedSize scan_index = order[position];
      if (scan_index >= first_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      const Size spectra_done = count_spectra.fetch_add(1, std::memory_order_relaxed) + 1;
      IF_MASTERTHREAD { setProgress(spectra_done); }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      const TheoreticalSpectrumGenerator& spectrum_generator = generators.forSpectrum(exp_spectrum);
      ScoringScratch_& scratch = scoringScratch_();
      // The candidates of this spectrum: the query appends to the hits, so the hits of the thread's previous spectrum
      // are cleared first (clear() keeps the capacity)
      FragmentIndex::SpectrumMatchesTopN& top_sms = scratch.top_sms;
      top_sms.clear();
      // ions:by_activation: only electron-activated spectra are matched against the c and z+1 ions,
      // so that these ions do not change which candidates the other spectra keep
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_index] : exp_spectrum;
      fi.querySpectrum(query, db, top_sms, generators.electronIons(exp_spectrum));

      const bool peptide_mass_mode = fi.isPeptideMassMode();

      // Reused across the candidates of this spectrum and, per thread, across spectra. Avoids per-candidate heap
      // churn of a fresh PeakSpectrum + its DataArrays (TSG's add_metainfo fills
      // StringDataArrays with ion names — these are a notable allocation hot spot
      // when the candidate count per spectrum is in the tens/hundreds). Both generator calls below start by
      // clearing the spectrum, so what the previous candidate or spectrum left in it does not matter.
      PeakSpectrum& theo_spectrum = scratch.theo_spectrum;
      std::vector<uint8_t>& theo_ion_series = scratch.theo_ion_series;

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& sms_pep = fi.getPeptides()[sms.peptide_idx_];

        // Peptide-mass mode: the candidate is the prefix of the mother, with the variable modifications, that the index
        // matched to the precursor
        AASequence mod_candidate = peptide_mass_mode
          ? fi.reconstructRealizedSubSequence(sms_pep, db, sms.realized_length_, sms.subset_bitmask_)
          : fi.reconstructModifiedSequence(sms_pep, db);

        // Index construction already removes repeated occurrences within a
        // chunk. Across chunks, skip the same hypothesis before updating either
        // scores or pool statistics. Charge and isotope hypotheses stay distinct.
        if (deduplicate_chunks)
        {
          std::string key = mod_candidate.toString();
          key += "\t" + std::to_string(sms.precursor_charge_) + "\t" + std::to_string(sms.isotope_error_);
          if (! pool_stats[scan_index].seen_candidates.insert(std::move(key)).second) { continue; }
        }

        HyperScore::PSMDetail detail;
        double score;
        // The peaks of getSpectrum() with the ion series per peak instead of the ion names (no strings, no data arrays):
        // the same scores
        if (spectrum_generator.getIonSeriesSpectrum(theo_spectrum, theo_ion_series, mod_candidate, 1, scoringMaxCharge_(sms.precursor_charge_)))
        {
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                              theo_ion_series, mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                            theo_ion_series, detail);
        }
        else
        {
          // Clear peaks + data arrays (ion names / charges) before refilling for the
          // next candidate; getSpectrum appends to whatever is there. clear(false) keeps the capacity of the
          // peak vector (clear(true) would shrink it); the data arrays are emptied either way.
          theo_spectrum.clear(false);
          spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, scoringMaxCharge_(sms.precursor_charge_));
          // Note: TSG emits sorted output when add_metainfo=true (see the
          // sortByPositionPresorted() call at the tail of getSpectrum_); the extra
          // sortByPosition() pass here was a redundant O(N) scan per candidate.
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                              mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);
        }

        // Summarise the candidate before it can be dropped below or pruned at the
        // end of the loop: the pool-derived PSM features describe the whole search
        // space, so a candidate that scored 0 (no fragment matched) still counts
        // and still belongs in the null distribution the z-score is measured against.
        // Each scan_index is owned by exactly one thread, so this needs no guard.
        pool_stats[scan_index].add(score);

        if (score == 0) continue;

        AnnotatedHit_ ah;
        ah.sequence = std::move(mod_candidate);
        ah.score = score;
        // Account for the additional charge hypotheses in the ion-count fractions.
        double seq_length = static_cast<double>(ah.sequence.size()) * scoringMaxCharge_(sms.precursor_charge_);
        ah.prefix_fraction = static_cast<float>(detail.matched_prefix_ions / seq_length);
        ah.suffix_fraction = static_cast<float>(detail.matched_suffix_ions / seq_length);
        ah.mean_error = static_cast<float>(detail.mean_error);
        ah.matched_prefix_ions = static_cast<uint16_t>(detail.matched_prefix_ions);
        ah.matched_suffix_ions = static_cast<uint16_t>(detail.matched_suffix_ions);
        ah.isotope_error = sms.isotope_error_;
        ah.applied_charge = sms.precursor_charge_;
        ah.index_peptide = static_cast<UInt32>(sms.peptide_idx_);
        ah.delta_mass = 0.0;
        if (open_search_mode)
        {
          double theo_mh_plus = ah.sequence.getMZ(1);
          double exp_mz = exp_spectrum.getPrecursors()[0].getMZ();
          double exp_mh_plus = exp_mz * sms.precursor_charge_ - ((sms.precursor_charge_ - 1) * proton_mass_u);
          ah.delta_mass = exp_mh_plus - theo_mh_plus;
        }

        annotated_hits[scan_index].push_back(std::move(ah));
      }

      // Prune to top-N per spectrum to bound memory: up to
      // scoring:max_candidates_per_spectrum hits (each owning a heap AASequence)
      // would otherwise stay resident until postProcessHits_ truncates them.
      // Correct because all of this spectrum's candidates are appended above and
      // scores are independent of the chunk a candidate came from: a hit outside
      // the top-N here can never re-enter it.
      auto& hits = annotated_hits[scan_index];
#ifdef _OPENMP
      pool_stats[scan_index].scoring_thread = omp_get_thread_num();
#endif
      if (hits.size() > keep)
      {
        std::partial_sort(hits.begin(), hits.begin() + keep, hits.end(), AnnotatedHit_::hasBetterScore);
        hits.resize(keep);
        hits.shrink_to_fit();
      }
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_scoring_error)
        {
          if (scan_index < first_failed.load(std::memory_order_relaxed))
          {
            scoring_error = std::current_exception();
            first_failed.store(scan_index, std::memory_order_relaxed);
          }
        }
      }
    }
    }
    endProgress();
    if (scoring_error) std::rethrow_exception(scoring_error);
  }

  std::vector<SignedSize> ProSEAlgorithm::scoringOrder_(const PeakMap& spectra, const PeakMap* query_spectra)
  {
    // Key: the neutral precursor mass of the first precursor of the candidate-retrieval spectrum (an unknown or
    // non-positive charge counts as 1; only the order depends on it). Spectra without a finite key (no precursor, a
    // non-finite m/z) go to the end in scan order. (key, scan index) is a total order: any sort gives the same order.
    const SignedSize n = static_cast<SignedSize>(spectra.size());
    std::vector<SignedSize> order(n);
    if (n < 2)
    {
      for (SignedSize i = 0; i < n; ++i) order[i] = i;
      return order;
    }
    std::vector<std::pair<double, SignedSize>> keyed(n);
    std::vector<unsigned char> valid(n);
#pragma omp parallel for schedule(static) if (n >= 16384)
    for (SignedSize i = 0; i < n; ++i)
    {
      const MSSpectrum& s = query_spectra != nullptr ? (*query_spectra)[i] : spectra[i];
      double key = std::numeric_limits<double>::quiet_NaN();
      if (!s.getPrecursors().empty())
      {
        const Precursor& precursor = s.getPrecursors()[0];
        const int charge = precursor.getCharge() > 0 ? precursor.getCharge() : 1;
        key = (precursor.getMZ() - Constants::PROTON_MASS_U) * charge;
      }
      valid[i] = std::isfinite(key) ? 1 : 0;
      keyed[i] = {valid[i] ? key : 0.0, i};
    }
    // the spectra with a key first (in scan order), then the others (in scan order)
    SignedSize n_valid = 0;
    for (SignedSize i = 0; i < n; ++i) { if (valid[i]) keyed[n_valid++] = keyed[i]; }
    SignedSize next_invalid = n_valid;
    for (SignedSize i = 0; i < n; ++i) { if (!valid[i]) order[next_invalid++] = i; }
    const auto by_key = [](const std::pair<double, SignedSize>& a, const std::pair<double, SignedSize>& b)
    {
      return a.first < b.first || (!(b.first < a.first) && a.second < b.second);
    };
#ifdef PROSE_GNU_PARALLEL
    if (n_valid >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(keyed.begin(), keyed.begin() + n_valid, by_key); }
    else
#endif
    { std::sort(keyed.begin(), keyed.begin() + n_valid, by_key); }
    for (SignedSize p = 0; p < n_valid; ++p) order[p] = keyed[p].second;
    return order;
  }

  namespace
  {
    // chimeric:enabled: the peaks of @p spectrum that the theoretical ions of @p sequence (fragment charges up to
    // @p max_charge) claim: every peak within the fragment tolerance of an ion. Merged into the ascending peak indices
    // @p claimed.
    void claimPeaks(const MSSpectrum& spectrum, const AASequence& sequence, const TheoreticalSpectrumGenerator& generator,
                    int max_charge, double fragment_tol, bool fragment_tol_ppm, PeakSpectrum& theo, std::vector<UInt32>& claimed)
    {
      theo.clear(true);
      generator.getSpectrum(theo, sequence, 1, max_charge);
      for (const Peak1D& ion : theo)
      {
        const double tol = fragment_tol_ppm ? ion.getMZ() * fragment_tol * 1e-6 : fragment_tol;
        for (auto it = spectrum.MZBegin(ion.getMZ() - tol); it != spectrum.end() && it->getMZ() <= ion.getMZ() + tol; ++it)
        {
          claimed.push_back(static_cast<UInt32>(it - spectrum.begin()));
        }
      }
      std::sort(claimed.begin(), claimed.end());
      claimed.erase(std::unique(claimed.begin(), claimed.end()), claimed.end());
    }
  }

  MSSpectrum ProSEAlgorithm::residualSpectrum_(const MSSpectrum& spectrum, const std::vector<UInt32>& removed, const Precursor& precursor)
  {
    std::vector<Size> kept;
    kept.reserve(spectrum.size());
    auto next_removed = removed.begin();
    for (Size p = 0; p < spectrum.size(); ++p)
    {
      if (next_removed != removed.end() && *next_removed == p)
      {
        ++next_removed;
        continue;
      }
      kept.push_back(p);
    }
    MSSpectrum residual = spectrum;
    residual.select(kept);
    residual.removeMetaValue(CoIsolationDetector::META_MZ);
    residual.removeMetaValue(CoIsolationDetector::META_CHARGE);
    residual.setPrecursors({precursor});
    return residual;
  }

  void ProSEAlgorithm::searchCoIsolated_(const PeakMap& spectra,
                                         const std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
                                         FragmentIndex& fi,
                                         const std::vector<FASTAFile::FASTAEntry>& db,
                                         const SpectrumGenerators_& generators,
                                         double fragment_tol,
                                         bool fragment_tol_ppm,
                                         bool open_search_mode,
                                         std::vector<CoIsolatedPSM_>& psms,
                                         std::vector<std::vector<AnnotatedHit_>>& hits,
                                         std::vector<CandidatePoolStats_>& pool_stats) const
  {
    psms.clear();
    hits.clear();
    pool_stats.clear();
    const Size max_n = static_cast<Size>(static_cast<int>(param_.getValue("chimeric:max_coisolated")));
    const Size max_hits = std::min(report_top_hits_, static_cast<Size>(static_cast<int>(param_.getValue("chimeric:top_hits"))));

    // The spectra with co-isolated precursors and at least one hit: their best hit, the peaks claimed so far and
    // the winners found so far (indices into the round results below)
    struct Scan
    {
      Size spectrum = 0;
      DoubleList mzs;
      IntList charges;
      const AnnotatedHit_* primary = nullptr;
      std::vector<UInt32> claimed;
      std::vector<std::pair<Size, Size>> winners; // (round, position in that round's results)
    };
    std::vector<Scan> scans;
    for (Size i = 0; i < spectra.size() && i < annotated_hits.size(); ++i)
    {
      const MSSpectrum& spectrum = spectra[i];
      if (annotated_hits[i].empty() || spectrum.empty() || spectrum.getPrecursors().empty()
          || !spectrum.metaValueExists(CoIsolationDetector::META_MZ) || !spectrum.metaValueExists(CoIsolationDetector::META_CHARGE))
      {
        continue;
      }
      Scan scan;
      scan.spectrum = i;
      scan.mzs = spectrum.getMetaValue(CoIsolationDetector::META_MZ);
      scan.charges = spectrum.getMetaValue(CoIsolationDetector::META_CHARGE);
      const Size n = std::min({scan.mzs.size(), scan.charges.size(), max_n});
      scan.mzs.resize(n);
      scan.charges.resize(n);
      if (n == 0) continue;
      scan.primary = &*std::min_element(annotated_hits[i].begin(), annotated_hits[i].end(), AnnotatedHit_::hasBetterScore);
      scans.push_back(std::move(scan));
    }
    if (scans.empty()) return;

#pragma omp parallel for schedule(dynamic, 64) default(none) shared(scans, spectra, generators, fragment_tol, fragment_tol_ppm)
    for (SignedSize s = 0; s < static_cast<SignedSize>(scans.size()); ++s)
    {
      Scan& scan = scans[s];
      const MSSpectrum& spectrum = spectra[scan.spectrum];
      PeakSpectrum theo;
      claimPeaks(spectrum, scan.primary->sequence, generators.forSpectrum(spectrum), scoringMaxCharge_(scan.primary->applied_charge),
                 fragment_tol, fragment_tol_ppm, theo, scan.claimed);
    }

    // One round per co-isolated precursor rank; a round's residual spectra are scored in batches
    constexpr Size batch_size = 8192;
    struct Result
    {
      Size scan = 0;
      CoIsolatedPSM_ psm;
      std::vector<AnnotatedHit_> hits; // the winner first
      CandidatePoolStats_ stats;
      bool emitted = false;
    };
    std::vector<std::vector<Result>> rounds;
    Size n_searched = 0, n_dropped = 0;
    PeakSpectrum theo;
    for (Size k = 0; k < max_n; ++k)
    {
      std::vector<Size> in_round;
      for (Size s = 0; s < scans.size(); ++s)
      {
        if (k < scans[s].mzs.size() && scans[s].charges[k] > 0) in_round.push_back(s);
      }
      if (in_round.empty()) break;
      std::vector<Result>& results = rounds.emplace_back();
      for (Size begin = 0; begin < in_round.size(); begin += batch_size)
      {
        const Size end = std::min(begin + batch_size, in_round.size());
        PeakMap residuals;
        std::vector<CoIsolatedPSM_> batch_psms;
        for (Size b = begin; b < end; ++b)
        {
          const Scan& scan = scans[in_round[b]];
          const MSSpectrum& spectrum = spectra[scan.spectrum];
          CoIsolatedPSM_ psm;
          psm.parent = scan.spectrum;
          psm.precursor = spectrum.getPrecursors()[0];
          psm.precursor.setMZ(scan.mzs[k]);
          psm.precursor.setCharge(scan.charges[k]);
          psm.removed = scan.claimed;
          residuals.addSpectrum(residualSpectrum_(spectrum, psm.removed, psm.precursor));
          batch_psms.push_back(std::move(psm));
        }
        std::vector<std::vector<AnnotatedHit_>> batch_hits(residuals.size());
        std::vector<CandidatePoolStats_> batch_stats(residuals.size());
        scoreSpectraAgainstIndex_(residuals, fi, db, generators, fragment_tol, fragment_tol_ppm, open_search_mode,
                                  batch_hits, batch_stats, "Scoring co-isolated precursors...");
        n_searched += residuals.size();
        residuals.clear(true);
        for (Size b = begin; b < end; ++b)
        {
          Scan& scan = scans[in_round[b]];
          // the hits that are not the spectrum's best peptide, the best (the winner) first
          std::vector<AnnotatedHit_>& candidates = batch_hits[b - begin];
          candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                          [&scan](const AnnotatedHit_& h) { return h.sequence == scan.primary->sequence; }),
                           candidates.end());
          if (candidates.empty()) continue;
          std::sort(candidates.begin(), candidates.end(), AnnotatedHit_::hasBetterScore);
          if (candidates.size() > max_hits) candidates.resize(max_hits);
          const auto winner = candidates.begin();
          const MSSpectrum& spectrum = spectra[scan.spectrum];
          claimPeaks(spectrum, winner->sequence, generators.forSpectrum(spectrum), scoringMaxCharge_(winner->applied_charge),
                     fragment_tol, fragment_tol_ppm, theo, scan.claimed);
          // dropped if the spectrum already reports this peptide at this charge: its distinct precursor m/z would
          // make it a precursor of its own in Percolator's target-decoy competition and count it twice
          auto same = [&winner](const AnnotatedHit_& h) { return h.applied_charge == winner->applied_charge && h.sequence == winner->sequence; };
          bool repeated = std::any_of(annotated_hits[scan.spectrum].begin(), annotated_hits[scan.spectrum].end(), same);
          for (const auto& [round, position] : scan.winners)
          {
            repeated = repeated || same(rounds[round][position].hits.front());
          }
          Result result;
          result.scan = in_round[b];
          result.psm = std::move(batch_psms[b - begin]);
          result.hits = std::move(candidates);
          result.stats = std::move(batch_stats[b - begin]);
          result.emitted = !repeated;
          if (repeated) ++n_dropped;
          scan.winners.emplace_back(k, results.size());
          results.push_back(std::move(result));
        }
      }
    }

    // In spectrum order, within a spectrum in detection order
    for (const Scan& scan : scans)
    {
      for (const auto& [round, position] : scan.winners)
      {
        Result& result = rounds[round][position];
        if (!result.emitted) continue;
        psms.push_back(std::move(result.psm));
        hits.push_back(std::move(result.hits));
        pool_stats.push_back(std::move(result.stats));
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << n_searched << " co-isolated precursors searched, " << psms.size()
                    << " with a PSM (" << n_dropped << " more repeat a PSM of their spectrum and are dropped)." << std::endl;
  }

  // =====================================================================
  // In-memory search: thin wrapper that builds a fresh SearchContext per call.
  // For repeated searches against the same database, prefer the
  // context-taking overload below.
  //
  // When database:chunk_size > 0 and the DB exceeds that size, the database
  // is split into chunks: each chunk builds its own FragmentIndex, all spectra
  // are scored against each chunk, and hits are accumulated across chunks
  // before a single postProcessHits_ + PeptideIndexing + FDR pass. This
  // trades speed (N × spectrum-scoring passes) for memory (only one chunk's
  // FragmentIndex in memory at a time).
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      PeakMap& spectra,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // Reset per-run stats here so the chunked single-file path (searchChunked_,
    // which does not route through search(spectra, ctx, ...)) starts clean too.
    last_run_stats_ = RunStatistics{};
    last_percolator_ = PercolatorRescoring{};
    inference_only_peptide_ids_ = PeptideIdentificationList();

    // Chunking disabled → take the existing single-context path (decoys built
    // lazily by prepareContext). The ctx is locally owned and not reused,
    // so opt in to eager FI release (M1) before PeptideIndexing.
    // ions:by_activation: index c and z+1 ions up front if electron-activated spectra are searched
    const bool electron_ions = countElectronActivated_(spectra) > 0;
    if (database_chunk_size_ == 0)
    {
      // single use: index only the peptides these spectra can reach
      std::function<const PeakMap*(Size)> searched_spectra;
      if (restrictIndexToSpectra_()) { searched_spectra = [&spectra](Size) { return &spectra; }; }
      SearchContext ctx = prepareContext_(std::vector<FASTAFile::FASTAEntry>(fasta_db), electron_ions, searched_spectra);
      ctx.release_fragment_index_after_scoring = true;
      return search(spectra, ctx, protein_ids, peptide_ids);
    }

    // Chunking configured: build the decoy-augmented DB once up-front, then
    // decide based on the AUGMENTED size. If the target DB is 3000 proteins
    // and chunk_size is 5000 with decoys enabled, the augmented DB (6000)
    // still exceeds chunk_size — decide-before-augment would have skipped
    // chunking and built a full FI 2× the user's declared memory budget.
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    auto full_db = buildDecoyAugmentedDB_(fasta_db, strategy);
    if (full_db.size() <= database_chunk_size_)
    {
      // Augmented DB fits in one chunk — use the single-context path but skip
      // prepareContext's internal decoy re-generation by building ctx inline.
      SearchContext ctx;
      ctx.db = std::move(full_db);
      ctx.decoy_string = strategy.decoy_string;
      ctx.decoy_is_prefix = strategy.is_prefix;
      ctx.have_decoys = strategy.have_decoys;
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      startProgress(0, 1, "Building fragment index...");
      ctx.fragment_index.setParameters(fragmentIndexParameters_(electron_ions));
      ctx.fragment_index.build(ctx.db);
      ctx.electron_ions = electron_ions;
      endProgress();
      reportDecoyParity_(ctx);
      return search(spectra, ctx, protein_ids, peptide_ids);
    }
    return searchChunked_(spectra, full_db, strategy, protein_ids, peptide_ids);
  }

  // =====================================================================
  // Chunked search implementation. Takes a pre-built decoy-augmented DB
  // (from buildDecoyAugmentedDB_) and splits it into chunks for scoring.
  // Called from the single-file search() wrapper and the multi-file wrapper.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::searchChunked_(
      PeakMap& spectra,
      std::vector<FASTAFile::FASTAEntry>& full_db,
      const DecoyStrategy_& strategy,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    const Size n_chunks = (full_db.size() + database_chunk_size_ - 1) / database_chunk_size_;
    if (chimeric_enabled_)
    {
      OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
    }
    OPENMS_LOG_INFO << "[ProSE] Database chunking enabled: " << full_db.size()
                    << " proteins (incl. decoys), chunk_size=" << database_chunk_size_
                    << " → " << n_chunks << " chunks." << std::endl;

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    bool open_search_mode = isOpenSearchMode_();
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSearchedSpectra_(spectra, "", ion_evidence, evidence_ptr, query_ptr);

    // ions:by_activation: the chunk indices hold c and z+1 ions if electron-activated spectra are searched
    const Size n_electron_activated = countElectronActivated_(spectra);
    const bool electron_ions = n_electron_activated > 0;
    if (electron_ions)
    {
      OPENMS_LOG_INFO << "[ProSE] " << n_electron_activated << " of " << spectra.size()
                      << " spectra are electron-activated (ETD, ECD, EThcD or ETciD) and are also scored"
                      << " with c and z+1 ions." << std::endl;
    }

    // Effective tolerances — may be narrowed by calibration below. Kept asymmetric
    // (lower / upper separately) so the FragmentIndex query can exploit the full
    // calibrated window rather than the looser max() collapse.
    double effective_precursor_tol_lower = precursor_mass_tolerance_lower_;
    double effective_precursor_tol_upper = precursor_mass_tolerance_upper_;
    double effective_fragment_tol = fragment_mass_tolerance_;

    // Save originals so we can restore algo-level members after search (they may be
    // mutated below for the OSMA mod-match tolerance computation).
    const double orig_prec_tol_lower = precursor_mass_tolerance_lower_;
    const double orig_prec_tol_upper = precursor_mass_tolerance_upper_;
    bool calibration_applied = false;
    // Restore algo-level tolerance members on every return path, and when an exception leaves this function (e.g. a
    // worker exception that the scoring loop rethrows). Calibration may mutate them below; leaving them mutated would
    // poison subsequent searches that reuse this ProSEAlgorithm instance (e.g. the multi-file wrapper, or a caller
    // that catches the exception).
    auto restore_tolerances = [&]()
    {
      if (calibration_applied)
      {
        calibration_applied = false;
        precursor_mass_tolerance_lower_ = orig_prec_tol_lower;
        precursor_mass_tolerance_upper_ = orig_prec_tol_upper;
      }
    };
    const OnScopeExit restore_on_exit(restore_tolerances);

    // Optional calibration on a strided sample of the full DB. The sample size is
    // bounded (see buildCalibrationSample_) so calibration memory stays O(chunk)
    // regardless of database_chunk_size_. Strided rather than first-N so small
    // chunk_size values don't starve the calibration pool (#9182).
    if (calibration_enabled_ && !open_search_mode)
    {
      std::vector<FASTAFile::FASTAEntry> cal_db = buildCalibrationSample_(full_db);
      FragmentIndex cal_fi;
      cal_fi.setParameters(fragmentIndexParameters_(electron_ions));
      cal_fi.build(cal_db);

      CalibrationResult_ cal = runCalibrationPass_(spectra, cal_fi, cal_db, query_ptr);
      if (cal.success)
      {
        effective_fragment_tol = cal.fragment_tolerance;
        if (!cal.extreme_bias)
        {
          effective_precursor_tol_lower = cal.cal_lower;
          effective_precursor_tol_upper = cal.cal_upper;
          precursor_mass_tolerance_lower_ = cal.cal_lower;
          precursor_mass_tolerance_upper_ = cal.cal_upper;
          calibration_applied = true;
          OPENMS_LOG_INFO << "[ProSE] Calibration (chunked, strided sample): shift=" << cal.precursor_shift
                          << " " << precursor_mass_tolerance_unit_
                          << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                          << " fragment=" << cal.fragment_tolerance << std::endl;
        }
        else
        {
          OPENMS_LOG_WARN << "[ProSE] Calibration: extreme bias, precursor calibration discarded. "
                          << "Fragment calibration applied." << std::endl;
        }
      }
    }

    // 3. Prepare spectrum generators (once).
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // 4. Allocate per-spectrum hit accumulator (persists across chunks).
    // not reserved: most spectra keep few hits or none, and report:top_hits per spectrum ahead of time is a lot
    vector<vector<AnnotatedHit_>> annotated_hits(spectra.size());
    // Accumulates across chunks: every chunk contributes its own candidates for the
    // same spectrum, and the pool features must see all of them.
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    // 5. Chunk loop.
    const Size chunk_size = database_chunk_size_;
    Size chunk_idx = 0;
    for (Size start = 0; start < full_db.size(); start += chunk_size)
    {
      ++chunk_idx;
      const Size end = std::min(start + chunk_size, full_db.size());

      OPENMS_LOG_INFO << "[ProSE] Chunk " << chunk_idx << ": proteins "
                      << start << "–" << (end - 1) << " (" << (end - start) << " proteins)" << std::endl;

      // Build FragmentIndex for this chunk only.
      std::vector<FASTAFile::FASTAEntry> chunk_db(full_db.begin() + start, full_db.begin() + end);
      FragmentIndex chunk_fi;
      {
        Param fi_params = fragmentIndexParameters_(electron_ions);
        // Apply calibrated tolerances (if calibration succeeded above). Asymmetric
        // lower/upper preserved — collapsing to max() would re-open the tight side
        // of the calibrated window and admit spurious decoy candidates.
        fi_params.setValue("fragment:mass_tolerance", effective_fragment_tol);
        fi_params.setValue("precursor:mass_tolerance_lower", effective_precursor_tol_lower);
        fi_params.setValue("precursor:mass_tolerance_upper", effective_precursor_tol_upper);
        chunk_fi.setParameters(fi_params);
      }
      chunk_fi.build(chunk_db);

      // Score all spectra against this chunk's index.
      scoreSpectraAgainstIndex_(spectra, chunk_fi, chunk_db, generators,
                                effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                                open_search_mode, annotated_hits, pool_stats,
                                "Scoring chunk " + StringUtils::toStr(chunk_idx) + "...", query_ptr);

      // Prune to top-N per spectrum after each chunk to bound memory growth.
      // Without this, K chunks × T candidates per spectrum could accumulate
      // K*T hits before postProcessHits_ truncates — problematic for 100M+
      // peptide databases with many chunks. The pruning is correct because
      // scores are independent across chunks: a hit that fails to place in the
      // current top-N cannot improve when more chunks are added.
      const Size keep = std::max(searchTopHits_(), Size(2)); // keep ≥2 for delta score
#pragma omp parallel for default(none) shared(annotated_hits, keep)
      for (SignedSize si = 0; si < (SignedSize)annotated_hits.size(); ++si)
      {
        if (annotated_hits[si].size() > keep)
        {
          std::partial_sort(annotated_hits[si].begin(),
                            annotated_hits[si].begin() + keep,
                            annotated_hits[si].end(),
                            AnnotatedHit_::hasBetterScore);
          annotated_hits[si].resize(keep);
        }
      }
    } // end chunk loop
    query_spectra.clear(true);
    for (auto& stats : pool_stats)
    {
      std::unordered_set<std::string>().swap(stats.seen_candidates);
    }

    // 6. Post-process merged hits (sort, annotate, PeptideIndexing, FDR).
    //    This runs once on ALL hits accumulated across all chunks.
    //    Set mod-match tolerance for open-search modification analysis downstream.
    last_mod_match_tolerance_used_ = computeModMatchTolerance_();

    startProgress(0, 1, "Post-processing PSMs...");
    postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      searchTopHits_(),
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      std::max(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_),
      effective_fragment_tol,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      "", // no database filename for in-memory search
      evidence_ptr
      );
    // Release the hits and the additional peak list before PeptideIndexing.
    freeOnOwnerThreads(annotated_hits, [&pool_stats](SignedSize i) { return pool_stats[i].scoring_thread; }); // (where scored)
    freeInParallel(pool_stats);
    freeInParallel(evidence_spectra);
    endProgress();

    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);

    // 7. PeptideIndexing against the FULL database (not per-chunk).
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    // Use the effective decoy marker/position that was actually searched. A
    // non-empty string (decoy_prefix_ when target-only) avoids PeptideIndexing's
    // own auto-detection; missing_decoy_action=silent keeps it quiet when none match.
    param_pi.setValue("decoy_string", strategy.decoy_string.empty() ? decoy_prefix_ : strategy.decoy_string);
    param_pi.setValue("decoy_string_position", strategy.is_prefix ? "prefix" : "suffix");
    param_pi.setValue("enzyme:name", enzyme_);
    param_pi.setValue("enzyme:specificity",
                      EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(full_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      restore_tolerances();
      if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
        return ExitCodes::INPUT_FILE_EMPTY;
      else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
        return ExitCodes::UNEXPECTED_RESULT;
      else
        return ExitCodes::UNKNOWN_ERROR;
    }

    // 8. PSM-level FDR only. Protein-level FDR is the caller's responsibility
    //    (matching the non-chunked search(spectra, ctx, ...) semantics — that
    //    method also does PSM FDR only; the file-based search() and the single-file
    //    searchWithModificationAnalysis overloads apply protein FDR post-call).
    const bool has_decoys = strategy.have_decoys;

    // Optional per-run ion priors need the target/decoy labels and the native scores: after
    // PeptideIndexing, before FDR overwrites the scores.
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Pre-FDR stats (target/decoy counts + HyperScore distribution).
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
    last_percolator_ = applyPsmFDR_(protein_ids, peptide_ids, has_decoys, last_run_stats_, nullptr, &inference_only_peptide_ids_);

    restore_tolerances();

    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // Protein mapping from the fragment index (peptide:protein_mapping = index):
  // what PeptideIndexing::run() writes, without PeptideIndexing (the spans of the hits from the index, or with a
  // protein-terminal modification from a pass over the database for the hits' sequences)
  // =====================================================================
  namespace
  {
    bool isAmbiguousResidue(const char c)
    {
      return c == 'B' || c == 'J' || c == 'Z' || c == 'X';
    }

    // The residues AhoCorasickAmbiguous matches an ambiguous protein residue with: B = D or N, J = I or L,
    // Z = E or Q, X = any of its unambiguous amino acids (every other letter, O and U included)
    const std::string& residuesMatchedBy(const char ambiguous)
    {
      static const std::string b = "DN", j = "IL", z = "EQ";
      static const std::string x = []
      {
        std::string all;
        for (char c = 'A'; c <= 'Z'; ++c)
        {
          if (!AA(c).isAmbiguous()) all += c;
        }
        return all;
      }();
      switch (ambiguous)
      {
        case 'B': return b;
        case 'J': return j;
        case 'Z': return z;
        default: return x;
      }
    }

    // The first @p length (<= 8) residues of @p residues as a number that sorts like them
    uint64_t residueKey(const std::string_view residues, const Size length)
    {
      uint64_t key = 0;
      for (Size i = 0; i < length; ++i) key = (key << 8) | static_cast<unsigned char>(residues[i]);
      return key;
    }

    using Occurrence = std::pair<UInt32, UInt32>; // {protein index, position}

    /**
      @brief The spans that step 3 of buildProteinMapping_() looks up for one ambiguous residue of a protein: those
      with no ambiguous residue before it and 1 to aaa_max ambiguous residues from it on (so every span over ambiguous
      residues belongs to its first one).

      With r_0 = position < r_1 < ... the ambiguous residues from the position on, a span [start, end) with start <=
      r_0 contains exactly k of them if r_{k-1} < end <= r_k (r_k = the protein end if there is none). Each span is
      looked up once per combination of the residues its ambiguous ones stand for. Used for both the count (before
      any lookup, to decide on the fallback) and the lookups, so the two agree.
    */
    class AmbiguousSpans
    {
    public:
      /// PeptideIndexing's maximum of aaa_max
      static constexpr Size MAX_AMBIGUOUS = 10;

      AmbiguousSpans(const std::string& protein, const Size position, const Size max_length, const Size aaa_max) :
        protein_(protein), first_start_(position), aaa_max_(std::min(aaa_max, MAX_AMBIGUOUS))
      {
        // the first start: after the previous ambiguous residue, and close enough to reach the position
        while (first_start_ > 0 && position - first_start_ + 1 < max_length && !isAmbiguousResidue(protein[first_start_ - 1]))
        {
          --first_start_;
        }
        // the ambiguous residues from the position on that a span can reach: up to aaa_max + 1 (the last one only
        // bounds the spans with aaa_max)
        const Size reach = std::min(protein.size(), position + max_length);
        ambiguous_[0] = position;
        for (Size i = position + 1; i < reach && count_ <= aaa_max_; ++i)
        {
          if (isAmbiguousResidue(protein[i])) ambiguous_[count_++] = i;
        }
      }

      /// The largest number of ambiguous residues in a span
      Size maxAmbiguous() const { return std::min(aaa_max_, count_); }

      /// The k-th ambiguous residue from the position on (k = 0: the position)
      Size ambiguous(const Size k) const { return ambiguous_[k]; }

      /// The starts [first, last] of the spans of @p length with exactly @p k (1..maxAmbiguous()) ambiguous residues
      /// (first > last if there are none)
      std::pair<SignedSize, SignedSize> starts(const Size length, const Size k) const
      {
        const SignedSize end_after = static_cast<SignedSize>(ambiguous_[k - 1]);
        const SignedSize end_at_most = static_cast<SignedSize>(k < count_ ? ambiguous_[k] : protein_.size());
        const SignedSize len = static_cast<SignedSize>(length);
        return {std::max(static_cast<SignedSize>(first_start_), end_after + 1 - len),
                std::min(static_cast<SignedSize>(ambiguous_[0]), end_at_most - len)};
      }

      /// The sequences looked up for the spans of @p length (combinations saturate at 2^32, far above any budget)
      Size lookups(const Size length) const
      {
        Size total = 0, combinations = 1;
        for (Size k = 1; k <= maxAmbiguous(); ++k)
        {
          combinations = std::min(combinations * residuesMatchedBy(protein_[ambiguous_[k - 1]]).size(), Size(1) << 32);
          const auto [first, last] = starts(length, k);
          if (first <= last) total += combinations * static_cast<Size>(last - first + 1);
        }
        return total;
      }

    private:
      const std::string& protein_;
      Size first_start_;
      Size aaa_max_;
      std::array<Size, MAX_AMBIGUOUS + 1> ambiguous_{};
      Size count_ = 1; ///< entries of ambiguous_
    };
  }

  ProSEAlgorithm::ProteinMapping_ ProSEAlgorithm::buildProteinMapping_(const FragmentIndex& index,
                                                                      const std::vector<FASTAFile::FASTAEntry>& db,
                                                                      std::vector<UInt32> candidates,
                                                                      const Param& indexer_parameters)
  {
    ProteinMapping_ mapping;

    // PeptideIndexing settings the mapping reproduces (see ProteinMapping_)
    if (static_cast<Int>(indexer_parameters.getValue("mismatches_max")) != 0
        || indexer_parameters.getValue("IL_equivalent").toBool()
        || !indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool()
        || indexer_parameters.getValue("write_protein_sequence").toBool()
        || indexer_parameters.getValue("write_protein_description").toBool()
        || indexer_parameters.getValue("keep_unreferenced_proteins").toBool()
        || indexer_parameters.getValue("decoy_string").toString().empty())
    {
      mapping.fallback_reason = "PeptideIndexing settings other than no mismatches, no I/L equivalence, N-terminal methionine "
                                "cleavage, no protein sequences or descriptions, no unreferenced proteins and a given decoy string";
      return mapping;
    }
    const std::string enzyme_name = indexer_parameters.getValue("enzyme:name").toString();
    const std::string specificity = indexer_parameters.getValue("enzyme:specificity").toString();
    const Param index_parameters = index.getParameters();
    if (specificity != EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL]
        || enzyme_name == EnzymaticDigestion::UnspecificCleavage || enzyme_name == EnzymaticDigestion::NoCleavage
        || index_parameters.getValue("enzyme").toString() != enzyme_name
        || index_parameters.getValue("peptide:enzyme_specificity").toString() != specificity)
    {
      mapping.fallback_reason = "not a fully specific digestion with the enzyme of the fragment index";
      return mapping;
    }
    // The spans with the residues of a hit come from the index entries if they list them all; with a protein-terminal
    // modification they do not (the entries of a span depend on its position in its protein), and the spans are found
    // by a pass over the database instead (1. below).
    const bool index_occurrences = index.hasProteinOccurrences(db);
    if (!index_occurrences && !index.isBuiltFrom(db))
    {
      mapping.fallback_reason = "the fragment index does not list every protein occurrence of its peptides "
                                "(low_memory, or not built from this database)";
      return mapping;
    }
    const Size aaa_max = static_cast<Size>(static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(enzyme_name);
    enzyme.setSpecificity(EnzymaticDigestion::SPEC_FULL);

    // The candidates, once each, and the lengths of their sequences. A hit with an ambiguous residue leaves the
    // mapping to PeptideIndexing; checked first, since it needs no pass over the database.
#ifdef PROSE_GNU_PARALLEL
    if (candidates.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(candidates.begin(), candidates.end()); }
    else
#endif
    { std::sort(candidates.begin(), candidates.end()); }
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) return mapping; // nothing to map (applyProteinMapping_() decides)
    const std::vector<FragmentIndex::Peptide>& peptides = index.getPeptides();
    if (candidates.back() >= peptides.size())
    {
      mapping.fallback_reason = "a candidate the fragment index does not hold";
      return mapping;
    }
    const auto residues_of = [&](const UInt32 candidate)
    {
      const FragmentIndex::Peptide& peptide = peptides[candidate];
      return std::string_view(db[peptide.protein_idx].sequence).substr(peptide.sequence_.first, peptide.sequence_.second);
    };
    bool unambiguous_candidates = true;
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : unambiguous_candidates)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      const std::string_view residues = residues_of(candidates[i]);
      unambiguous_candidates = unambiguous_candidates && std::none_of(residues.begin(), residues.end(), isAmbiguousResidue);
    }
    if (!unambiguous_candidates)
    {
      mapping.fallback_reason = "a hit with an ambiguous residue";
      return mapping;
    }
    Size min_length = std::numeric_limits<Size>::max(), max_length = 0;
    for (const UInt32 candidate : candidates)
    {
      min_length = std::min<Size>(min_length, peptides[candidate].sequence_.second);
      max_length = std::max<Size>(max_length, peptides[candidate].sequence_.second);
    }
    std::vector<uint8_t> length_used(max_length + 1, 0);
    for (const UInt32 candidate : candidates) length_used[peptides[candidate].sequence_.second] = 1;

    // The database: letters A-Z only (PeptideIndexing removes '*' and skips other symbols, which shifts positions),
    // no stretch of more than aaa_max X (PeptideIndexing splits proteins at such stretches), and at most a budget of
    // sequences that 3. looks up for the spans over ambiguous residues (beyond it, PeptideIndexing does this better).
    // The lookups are counted here, without looking anything up, so that a fallback costs little more than this pass.
    constexpr Size lookup_budget = Size(1) << 22;
    std::array<uint8_t, 256> residue_class{}; // 0: unambiguous letter, 1: ambiguous letter, 2: anything else
    residue_class.fill(2);
    for (char c = 'A'; c <= 'Z'; ++c) residue_class[static_cast<unsigned char>(c)] = isAmbiguousResidue(c) ? 1 : 0;
    std::vector<uint8_t> with_spans(db.size(), 0); // has spans for 3.
    bool letters_only = true, short_stretches = true;
    std::atomic<Size> lookups{0}; // exact while at most lookup_budget (then counting stops)
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : letters_only, short_stretches)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      const std::string& protein = db[p].sequence;
      uint8_t classes = 0; // OR of its residue classes
      for (const char c : protein) classes |= residue_class[static_cast<unsigned char>(c)];
      letters_only = letters_only && classes < 2;
      if (classes != 1) continue;
      Size stretch = 0;
      for (const char c : protein)
      {
        stretch = (c == 'X') ? stretch + 1 : 0;
        short_stretches = short_stretches && stretch <= aaa_max;
      }
      Size protein_lookups = 0;
      for (Size position = 0; position < protein.size(); ++position)
      {
        if (!isAmbiguousResidue(protein[position])) continue;
        if (protein_lookups > lookup_budget || lookups.load(std::memory_order_relaxed) > lookup_budget) break;
        const AmbiguousSpans spans(protein, position, max_length, aaa_max);
        for (Size length = min_length; length <= max_length; ++length)
        {
          if (length_used[length]) protein_lookups += spans.lookups(length);
        }
      }
      with_spans[p] = protein_lookups > 0;
      lookups.fetch_add(std::min(protein_lookups, lookup_budget + 1), std::memory_order_relaxed);
    }
    if (!letters_only || !short_stretches)
    {
      mapping.fallback_reason = letters_only ? "the database has stretches of more than aaa_max X"
                                             : "the database has symbols other than the letters A-Z";
      return mapping;
    }
    if (lookups.load() > lookup_budget)
    {
      mapping.fallback_reason = "too many spans over ambiguous residues in the database";
      return mapping;
    }

    // 1. the spans of the digest: per candidate, then per sequence (from the index; else after the sequences, below)
    std::vector<std::vector<Occurrence>> candidate_occurrences(candidates.size());
    if (index_occurrences)
    {
#pragma omp parallel for schedule(dynamic, 64)
      for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
      {
        index.getProteinOccurrences(candidates[i], db, candidate_occurrences[i]);
      }
    }
    // One entry per sequence (the candidates of a sequence differ in their modifications only: they have the same
    // spans). The sequences go to ProteinMapping_::SHARDS hash tables by their hash, filled in parallel; the entries
    // are numbered table by table, in the order of the candidates.
    constexpr Size shards = ProteinMapping_::SHARDS;
    std::vector<Size> shard_of(candidates.size());
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      shard_of[i] = std::hash<std::string_view>{}(residues_of(candidates[i])) % shards;
    }
    std::vector<std::vector<Size>> shard_candidates(shards);
    for (Size i = 0; i < candidates.size(); ++i) shard_candidates[shard_of[i]].push_back(i);
    mapping.entries.resize(shards);
    std::vector<std::vector<Size>> shard_entries(shards); // per table: the candidate of each of its entries
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      auto& table = mapping.entries[shard];
      table.reserve(shard_candidates[shard].size());
      for (const Size i : shard_candidates[shard])
      {
        const auto [entry, added] = table.try_emplace(residues_of(candidates[i]), shard_entries[shard].size());
        if (added)
        {
          shard_entries[shard].push_back(i);
        }
        else
        { // another modified form of the same residues: the same spans (merged and made unique below)
          std::vector<Occurrence>& first = candidate_occurrences[shard_entries[shard][entry->second]];
          first.insert(first.end(), candidate_occurrences[i].begin(), candidate_occurrences[i].end());
        }
      }
    }
    std::vector<Size> shard_start(shards + 1, 0);
    for (Size shard = 0; shard < shards; ++shard) shard_start[shard + 1] = shard_start[shard] + shard_entries[shard].size();
    mapping.sequences.resize(shard_start[shards]);
    mapping.occurrences.resize(shard_start[shards]);
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      for (auto& [sequence, entry] : mapping.entries[shard]) entry += shard_start[shard];
      for (Size k = 0; k < shard_entries[shard].size(); ++k)
      {
        const Size i = shard_entries[shard][k];
        mapping.sequences[shard_start[shard] + k] = residues_of(candidates[i]);
        mapping.occurrences[shard_start[shard] + k] = std::move(candidate_occurrences[i]);
      }
    }

    std::vector<std::tuple<Size, UInt32, UInt32>> extra; // spans of 1. (database pass), 2. and 3.: {sequence entry, protein, position}

    // 1. without the occurrences of the index: every span of the database with the residues of a sequence that
    //    PeptideIndexing accepts (isValidProduct() as in its search, which matches the residues literally here: the
    //    spans over ambiguous residues are those of 3.). Positions are prefiltered by their first residues (a bit set over
    //    the sequences' first prefix_length residues as a base-26 number), then looked up by the key of their first
    //    key_length residues. The database holds the letters A-Z only (checked above).
    if (!index_occurrences)
    {
      // isValidProduct() tokenizes with the regex per call (~20 us); a simple rule is evaluated directly
      SegmentRule rule;
      const bool simple_rule = parseSegmentRule(ProteaseDB::getInstance()->getEnzyme(enzyme_name)->getRegEx(), rule);
      const Size key_length = std::min<Size>(min_length, 8);
      std::vector<std::pair<uint64_t, Size>> keyed(mapping.sequences.size());
#pragma omp parallel for schedule(static)
      for (SignedSize entry = 0; entry < static_cast<SignedSize>(keyed.size()); ++entry)
      {
        keyed[entry] = {residueKey(mapping.sequences[entry], key_length), static_cast<Size>(entry)};
      }
#ifdef PROSE_GNU_PARALLEL
      if (keyed.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::sort(keyed.begin(), keyed.end()); } // (distinct pairs)
      else
#endif
      { std::sort(keyed.begin(), keyed.end()); }
      const Size prefix_length = std::min<Size>(key_length, 5);
      Size prefix_count = 1;
      for (Size i = 0; i < prefix_length; ++i) prefix_count *= 26;
      const Size prefix_shift = prefix_count / 26; // 26^(prefix_length - 1)
      const auto prefix_of = [prefix_length](const char* residues)
      {
        Size code = 0;
        for (Size i = 0; i < prefix_length; ++i) code = code * 26 + static_cast<Size>(residues[i] - 'A');
        return code;
      };
      std::vector<uint64_t> prefix_bits((prefix_count + 63) / 64, 0);
      for (const std::string_view sequence : mapping.sequences)
      {
        const Size code = prefix_of(sequence.data());
        prefix_bits[code >> 6] |= uint64_t(1) << (code & 63);
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme; // as PeptideIndexing: one per thread
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
#pragma omp for schedule(dynamic, 16) nowait
        for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
        {
          const std::string& protein = db[p].sequence;
          if (protein.size() < min_length) continue;
          Size code = prefix_of(protein.data()) / 26; // the first prefix_length - 1 residues
          for (Size position = 0; position + min_length <= protein.size(); ++position)
          {
            code = (code % prefix_shift) * 26 + static_cast<Size>(protein[position + prefix_length - 1] - 'A');
            if (((prefix_bits[code >> 6] >> (code & 63)) & 1) == 0) continue;
            const uint64_t key = residueKey(std::string_view(protein).substr(position), key_length);
            for (auto k = std::lower_bound(keyed.begin(), keyed.end(), std::make_pair(key, Size(0)));
                 k != keyed.end() && k->first == key; ++k)
            {
              const std::string_view sequence = mapping.sequences[k->second];
              if (position + sequence.size() <= protein.size() && protein.compare(position, sequence.size(), sequence) == 0
                  && (simple_rule ? isValidFullySpecific(rule, protein, position, sequence.size())
                                  : thread_enzyme.isValidProduct(protein, static_cast<int>(position), static_cast<int>(sequence.size()), true, true, false)))
              {
                found.emplace_back(k->second, static_cast<UInt32>(p), static_cast<UInt32>(position));
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // 2. positions 1 and 2 of proteins that start with M: found through their first residues
    {
      // {key of the residues from the position on, protein, position}, sorted. Filled in parallel by chunks of
      // proteins, grouped by the first residue (the highest byte of the key), then each group sorted on its own.
      const Size key_length = std::min<Size>(min_length, 8);
      using Start = std::tuple<uint64_t, UInt32, UInt32>;
      const auto is_start = [&](const std::string& protein, Size position)
      {
        return position + min_length <= protein.size() && protein[0] == 'M';
      };
      constexpr Size letters = 26;
      const Size num_chunks = std::max<Size>(1, std::min<Size>(256, db.size() / 1024));
      std::vector<Size> offset(num_chunks * letters + 1, 0); // [letter][chunk]: count, then where the chunk fills it
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (is_start(protein, position)) ++offset[(protein[position] - 'A') * num_chunks + chunk + 1];
          }
        }
      }
      for (Size i = 1; i < offset.size(); ++i) offset[i] += offset[i - 1];
      std::vector<Start> starts(offset.back());
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (!is_start(protein, position)) continue;
            starts[offset[(protein[position] - 'A') * num_chunks + chunk]++] =
              Start(residueKey(std::string_view(protein).substr(position), key_length), static_cast<UInt32>(p), static_cast<UInt32>(position));
          }
        }
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize letter = 0; letter < static_cast<SignedSize>(letters); ++letter)
      {
        // after the fill, offset[letter * num_chunks + num_chunks - 1] is where the letter's group ends
        const Size begin = letter == 0 ? 0 : offset[letter * num_chunks - 1];
        const Size end = offset[(letter + 1) * num_chunks - 1];
        std::sort(starts.begin() + begin, starts.begin() + end);
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme; // as PeptideIndexing: one per thread
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
#pragma omp for schedule(dynamic, 256) nowait
        for (SignedSize entry = 0; entry < static_cast<SignedSize>(mapping.sequences.size()); ++entry)
        {
          const std::string_view sequence = mapping.sequences[entry];
          const uint64_t key = residueKey(sequence, key_length);
          for (auto s = std::lower_bound(starts.begin(), starts.end(), std::make_tuple(key, UInt32(0), UInt32(0)));
               s != starts.end() && std::get<0>(*s) == key; ++s)
          {
            const UInt32 protein = std::get<1>(*s), position = std::get<2>(*s);
            const std::string& protein_sequence = db[protein].sequence;
            if (position + sequence.size() <= protein_sequence.size() && protein_sequence.compare(position, sequence.size(), sequence) == 0
                && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(position), static_cast<int>(sequence.size()), true, true, false))
            {
              found.emplace_back(static_cast<Size>(entry), protein, position);
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // 3. spans over ambiguous residues (AmbiguousSpans; within the budget checked above): each looked up with every
    //    combination of the residues its ambiguous ones stand for
    {
      // work items: an ambiguous residue and a length with spans (a few residues with several ambiguous ones around
      // them would otherwise be the critical path); at most one per two lookups, so bounded by the budget
      struct Item { UInt32 protein, position, length; };
      std::vector<Item> items;
#pragma omp parallel
      {
        std::vector<Item> thread_items;
#pragma omp for schedule(dynamic, 16) nowait
        for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
        {
          if (!with_spans[p]) continue;
          const std::string& protein = db[p].sequence;
          for (Size position = 0; position < protein.size(); ++position)
          {
            if (!isAmbiguousResidue(protein[position])) continue;
            const AmbiguousSpans spans(protein, position, max_length, aaa_max);
            for (Size length = min_length; length <= max_length; ++length)
            {
              if (length_used[length] && spans.lookups(length) > 0)
              {
                thread_items.push_back({static_cast<UInt32>(p), static_cast<UInt32>(position), static_cast<UInt32>(length)});
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        items.insert(items.end(), thread_items.begin(), thread_items.end());
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme;
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
        std::string window;
        std::array<Size, AmbiguousSpans::MAX_AMBIGUOUS> choice{};
#pragma omp for schedule(dynamic, 1) nowait
        for (SignedSize i = 0; i < static_cast<SignedSize>(items.size()); ++i)
        {
          const Item& item = items[i];
          const std::string& protein_sequence = db[item.protein].sequence;
          const AmbiguousSpans spans(protein_sequence, item.position, max_length, aaa_max);
          for (Size k = 1; k <= spans.maxAmbiguous(); ++k)
          {
            const auto [first, last] = spans.starts(item.length, k);
            for (SignedSize start = first; start <= last; ++start)
            {
              window.assign(protein_sequence, static_cast<Size>(start), item.length);
              std::fill_n(choice.begin(), k, 0);
              while (true)
              {
                for (Size w = 0; w < k; ++w)
                {
                  window[spans.ambiguous(w) - static_cast<Size>(start)] = residuesMatchedBy(protein_sequence[spans.ambiguous(w)])[choice[w]];
                }
                const Size entry = mapping.find(window);
                if (entry < mapping.sequences.size()
                    && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(start), static_cast<int>(item.length), true, true, false))
                {
                  found.emplace_back(entry, item.protein, static_cast<UInt32>(start));
                }
                Size w = 0; // the next combination
                for (; w < k; ++w)
                {
                  if (++choice[w] < residuesMatchedBy(protein_sequence[spans.ambiguous(w)]).size()) break;
                  choice[w] = 0;
                }
                if (w == k) break;
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // all spans of every sequence, ascending and once (2. and 3. may find spans of 1. again)
    // (grouped by entry through counts, then appended per entry in parallel: in the order of `extra` within an entry,
    // which the sort below makes irrelevant)
    if (!extra.empty())
    {
      std::vector<Size> first(mapping.occurrences.size() + 1, 0);
      for (const auto& found : extra) ++first[std::get<0>(found) + 1];
      for (Size k = 1; k < first.size(); ++k) first[k] += first[k - 1];
      std::vector<Occurrence> grouped(extra.size());
      {
        std::vector<Size> next(first.begin(), first.end() - 1);
        for (const auto& [entry, protein, position] : extra) grouped[next[entry]++] = Occurrence(protein, position);
      }
      std::vector<std::tuple<Size, UInt32, UInt32>>().swap(extra);
      FirstError error; // (an allocation failure)
#pragma omp parallel for schedule(dynamic, 256)
      for (SignedSize k = 0; k < static_cast<SignedSize>(mapping.occurrences.size()); ++k)
      {
        if (first[k] == first[k + 1]) continue;
        try
        {
          mapping.occurrences[k].insert(mapping.occurrences[k].end(), grouped.begin() + first[k], grouped.begin() + first[k + 1]);
        }
        catch (...)
        {
          error.set(k);
        }
      }
      error.rethrow();
    }
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize k = 0; k < static_cast<SignedSize>(mapping.occurrences.size()); ++k)
    {
      std::vector<Occurrence>& occurrences = mapping.occurrences[k];
      std::sort(occurrences.begin(), occurrences.end());
      occurrences.erase(std::unique(occurrences.begin(), occurrences.end()), occurrences.end());
    }
    return mapping;
  }

  bool ProSEAlgorithm::applyProteinMapping_(const ProteinMapping_& mapping,
                                            const std::vector<FASTAFile::FASTAEntry>& db,
                                            const Param& indexer_parameters,
                                            std::vector<ProteinIdentification>& protein_ids,
                                            PeptideIdentificationList& peptide_ids,
                                            std::string& fallback_reason)
  {
    if (!mapping.fallback_reason.empty())
    {
      fallback_reason = mapping.fallback_reason;
      return false;
    }
    if (protein_ids.size() != 1)
    {
      fallback_reason = "not exactly one identification run";
      return false;
    }
    {
      // PeptideIndexing changes its cleavage rules for results of X! Tandem and MS-GF+
      std::string engine = protein_ids[0].getOriginalSearchEngineName();
      StringUtils::toUpper(engine);
      const ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
      if (engine == "XTANDEM" || engine == "MS-GF+" || engine == "MSGFPLUS"
          || search_parameters.metaValueExists("SE:XTandem") || search_parameters.metaValueExists("SE:MS-GF+"))
      {
        fallback_reason = "results of X! Tandem or MS-GF+";
        return false;
      }
    }

    // the sequence entry of every hit (PeptideIndexing maps the unmodified sequence)
    std::vector<Size> first_hit(peptide_ids.size() + 1, 0);
    for (Size i = 0; i < peptide_ids.size(); ++i) first_hit[i + 1] = first_hit[i] + peptide_ids[i].getHits().size();
    if (first_hit.back() == 0)
    {
      fallback_reason = "no peptide hits"; // PeptideIndexing's own handling (PEPTIDE_IDS_EMPTY)
      return false;
    }
    std::vector<Size> hit_entry(first_hit.back());
    bool all_mapped = true;
#pragma omp parallel for schedule(dynamic, 64) reduction(&& : all_mapped)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const Size entry = mapping.find(hits[h].getSequence().toUnmodifiedString());
        const bool mapped = entry < mapping.sequences.size() && !mapping.occurrences[entry].empty();
        all_mapped = all_mapped && mapped;
        if (mapped) hit_entry[first_hit[i] + h] = entry;
      }
    }
    if (!all_mapped)
    {
      fallback_reason = "a hit whose sequence was not mapped from the fragment index";
      return false;
    }

    const std::string decoy_string = indexer_parameters.getValue("decoy_string").toString();
    const bool decoy_prefix = indexer_parameters.getValue("decoy_string_position").toString() == "prefix";
    std::vector<uint8_t> is_decoy(db.size());
#pragma omp parallel for schedule(static)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      is_decoy[p] = decoy_prefix ? StringUtils::hasPrefix(db[p].identifier, decoy_string) : StringUtils::hasSuffix(db[p].identifier, decoy_string);
    }
    // the proteins of the hits
    std::vector<uint8_t> entry_used(mapping.occurrences.size(), 0);
    for (const Size entry : hit_entry) entry_used[entry] = 1;
    std::vector<uint8_t> protein_used(db.size(), 0);
    for (Size entry = 0; entry < entry_used.size(); ++entry)
    {
      if (!entry_used[entry]) continue;
      for (const Occurrence& occurrence : mapping.occurrences[entry]) protein_used[occurrence.first] = 1;
    }
    std::vector<UInt32> proteins; // in database order (PeptideIndexing's std::set of protein indices)
    bool any_decoy = false;
    for (Size p = 0; p < db.size(); ++p)
    {
      if (!protein_used[p]) continue;
      proteins.push_back(static_cast<UInt32>(p));
      any_decoy = any_decoy || is_decoy[p];
    }
    if (!any_decoy && indexer_parameters.getValue("missing_decoy_action").toString() != "silent")
    {
      fallback_reason = "no hit maps to a decoy"; // PeptideIndexing's warning or error
      return false;
    }

    // The meta value names in the order in which PeptideIndexing registers them (its first hit is mapped): MetaInfo
    // keeps, and idXML writes, the values of a hit ordered by registry index
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt target_decoy = registry.registerName("target_decoy");
    const UInt protein_references = registry.registerName("protein_references");
    const DataValue target("target"), decoy("decoy"), target_and_decoy("target+decoy"), unique("unique"), non_unique("non-unique");

    // the evidences of every hit, in PeptideIndexing's order (protein, position)
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        PeptideHit& hit = hits[h];
        const std::vector<Occurrence>& occurrences = mapping.occurrences[hit_entry[first_hit[i] + h]];
        const Size length = hit.getSequence().size();
        std::vector<PeptideEvidence> evidences;
        evidences.reserve(occurrences.size());
        bool in_target = false, in_decoy = false;
        Size protein_count = 0;
        UInt32 last_protein = std::numeric_limits<UInt32>::max();
        for (const auto& [protein, position] : occurrences)
        {
          const std::string& sequence = db[protein].sequence;
          const char before = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : sequence[position - 1];
          const char after = (position + length >= sequence.size()) ? PeptideEvidence::C_TERMINAL_AA : sequence[position + length];
          evidences.emplace_back(db[protein].identifier, static_cast<Int>(position), static_cast<Int>(position + length) - 1, before, after);
          if (protein != last_protein)
          {
            last_protein = protein;
            ++protein_count;
          }
          (is_decoy[protein] ? in_decoy : in_target) = true;
        }
        hit.setPeptideEvidences(std::move(evidences));
        hit.setMetaValue(target_decoy, (in_target && in_decoy) ? target_and_decoy : (in_target ? target : decoy));
        hit.setMetaValue(protein_references, protein_count == 1 ? unique : non_unique);
      }
    }

    // the protein hits: the proteins of the hits (keep_unreferenced_proteins = false), in database order
    std::vector<ProteinHit> protein_hits(proteins.size());
#pragma omp parallel for schedule(static)
    for (SignedSize k = 0; k < static_cast<SignedSize>(proteins.size()); ++k)
    {
      protein_hits[k].setAccession(db[proteins[k]].identifier);
      protein_hits[k].setMetaValue(target_decoy, is_decoy[proteins[k]] ? decoy : target);
    }
    protein_ids[0].getHits() = std::move(protein_hits);

    // PeptideIndexing's settings, as it records them
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(indexer_parameters.getValue("enzyme:name").toString());
    enzyme.setSpecificity(ProteaseDigestion::getSpecificityByName(indexer_parameters.getValue("enzyme:specificity").toString()));
    ProteinIdentification::SearchParameters search_parameters = protein_ids[0].getSearchParameters();
    search_parameters.setMetaValue("PeptideIndexer:decoy_string", decoy_string);
    search_parameters.setMetaValue("PeptideIndexer:decoy_string_position", decoy_prefix ? "prefix" : "suffix");
    search_parameters.setMetaValue("PeptideIndexer:enzyme", enzyme.getEnzymeName());
    search_parameters.setMetaValue("PeptideIndexer:enzyme_specificity", EnzymaticDigestion::NamesOfSpecificity[enzyme.getSpecificity()]);
    search_parameters.setMetaValue("PeptideIndexer:aaa_max", static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    search_parameters.setMetaValue("PeptideIndexer:mismatches_max", static_cast<Int>(indexer_parameters.getValue("mismatches_max")));
    search_parameters.setMetaValue("PeptideIndexer:IL_equivalent", indexer_parameters.getValue("IL_equivalent").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:allow_nterm_protein_cleavage",
                                   indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:unmatched_action", indexer_parameters.getValue("unmatched_action").toString());
    search_parameters.setMetaValue("PeptideIndexer:missing_decoy_action", indexer_parameters.getValue("missing_decoy_action").toString());
    protein_ids[0].setSearchParameters(std::move(search_parameters));
    return true;
  }

  // =====================================================================
  // In-memory search using a pre-built SearchContext (no index rebuild).
  // Takes the context by non-const reference because the underlying
  // FragmentIndex::querySpectrum() and PeptideIndexing::run() APIs are both
  // non-const, even though neither conceptually mutates shared state.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      PeakMap& spectra,
      SearchContext& ctx,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // Reset the per-run statistics bridge for this file. Callers copy it into
    // their SearchResult::stats after search() returns OK. The same for the pre-FDR PSMs: a copy
    // left by an earlier search must not end up in this file's result.
    last_run_stats_ = RunStatistics{};
    last_percolator_ = PercolatorRescoring{};
    pre_fdr_peptide_ids_ = PeptideIdentificationList();
    inference_only_peptide_ids_ = PeptideIdentificationList();

    // chimeric:enabled: find the precursors co-isolated with the selected one in the survey scans, which are not
    // searched themselves
    if (chimeric_enabled_ && !chimericSearch_())
    {
      OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
    }
    if (chimericSearch_())
    {
      std::vector<MSSpectrum>& all = spectra.getSpectra();
      const bool survey_scans = std::any_of(all.begin(), all.end(), [](const MSSpectrum& s) { return s.getMSLevel() == 1; });
      if (survey_scans)
      {
        const Size found = CoIsolationDetector::annotate(spectra, precursor_min_charge_, precursor_max_charge_,
                                                         param_.getValue("chimeric:mz_tolerance"), true,
                                                         param_.getValue("chimeric:max_kl"),
                                                         static_cast<Size>(static_cast<int>(param_.getValue("chimeric:max_coisolated"))));
        all.erase(std::remove_if(all.begin(), all.end(), [](const MSSpectrum& s) { return s.getMSLevel() == 1; }), all.end());
        OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << found << " co-isolated precursors in the isolation windows of "
                        << spectra.size() << " spectra." << std::endl;
      }
      else if (std::none_of(all.begin(), all.end(), [](const MSSpectrum& s) { return s.metaValueExists(CoIsolationDetector::META_MZ); }))
      {
        OPENMS_LOG_WARN << "[ProSE] chimeric:enabled needs the survey (MS1) scans, which the spectra do not include;"
                        << " no co-isolated precursors are searched." << std::endl;
      }
    }

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    bool open_search = isOpenSearchMode_();
    OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search ? "true" : "false")
                    << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] "
                    << precursor_mass_tolerance_unit_ << ")" << std::endl;

    takeBackThreadsOfEndedRead(); // here and between the other phases: see BackgroundSpectra
    startProgress(0, 1, "Filtering spectra...");
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSearchedSpectra_(spectra, "", ion_evidence, evidence_ptr, query_ptr);
    endProgress();

    // ions:by_activation: electron-activated spectra are also scored with c and z+1 ions, so the
    // index must hold them. A context without them is left unchanged, as other searches may share
    // it; this call builds its own index instead.
    const Size n_electron_activated = countElectronActivated_(spectra);
    FragmentIndex electron_index;
    FragmentIndex* index = &ctx.fragment_index;
    if (n_electron_activated > 0)
    {
      OPENMS_LOG_INFO << "[ProSE] " << n_electron_activated << " of " << spectra.size()
                      << " spectra are electron-activated (ETD, ECD, EThcD or ETciD) and are also scored"
                      << " with c and z+1 ions." << std::endl;
      if (!ctx.electron_ions)
      {
        OPENMS_LOG_WARN << "[ProSE] The prepared fragment index holds no c and z+1 ions; building one for this"
                        << " search. prepareContext(fasta_db, true) prepares a context that has them." << std::endl;
        startProgress(0, 1, "Building fragment index with c and z+1 ions...");
        electron_index.setParameters(fragmentIndexParameters_(true));
        electron_index.build(ctx.db);
        endProgress();
        index = &electron_index;
      }
    }

    // Reference the prepared (decoy-augmented) database and the fragment index to search.
    std::vector<FASTAFile::FASTAEntry>& db = ctx.db;
    FragmentIndex& fragment_index_ = *index;

    // Effective tolerances: may be overridden by calibration pass below.
    // The precursor scalar passed to postProcessHits_ is the widest bound
    // (legacy API); calibration may narrow it. Fragment tolerance is a
    // single scalar throughout.
    double effective_precursor_tol = std::max(precursor_mass_tolerance_lower_,
                                              precursor_mass_tolerance_upper_);
    double effective_fragment_tol = fragment_mass_tolerance_;

    // Save original FragmentIndex parameters AND algo-level tolerance members so we
    // can restore them after the search if calibration modifies query-time tolerances.
    // This avoids persistent mutation of the shared SearchContext and the algorithm
    // instance across multi-file calls (the multi-file wrapper reuses a single
    // ProSEAlgorithm instance, so member leaks corrupt later runs).
    const Param fi_params_original = fragment_index_.getParameters();
    const double orig_precursor_mass_tolerance_lower = precursor_mass_tolerance_lower_;
    const double orig_precursor_mass_tolerance_upper = precursor_mass_tolerance_upper_;
    bool fi_params_modified = false;
    // Restores the FragmentIndex parameters AND the algo-level tolerance members if calibration modified them, so the
    // shared SearchContext and the algorithm instance are both clean for subsequent per-file searches in the
    // multi-file wrapper: on the return paths below, and when an exception leaves this function (e.g. a worker
    // exception that the scoring or calibration loop rethrows) for a caller that catches it.
    auto restore_fi_params = [&]()
    {
      if (fi_params_modified)
      {
        fi_params_modified = false;
        precursor_mass_tolerance_lower_ = orig_precursor_mass_tolerance_lower;
        precursor_mass_tolerance_upper_ = orig_precursor_mass_tolerance_upper;
        fragment_index_.setParameters(fi_params_original);
      }
    };
    const OnScopeExit restore_on_exit([&restore_fi_params]() noexcept
    {
      try { restore_fi_params(); }
      catch (...) {} // only after another exception (e.g. std::bad_alloc copying the parameters): the members are restored
    });

    // --- Optional calibration pass ---
    if (calibration_enabled_ && !open_search)
    {
      startProgress(0, 1, "Running calibration pass...");
      StopWatch sw_cal; sw_cal.start();
      {
        OnScopeExit end_progress([this]() { endProgress(); }); // also when the calibration loop rethrows
        last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db, query_ptr);
      }
      sw_cal.stop();
      last_run_stats_.seconds_calibration = sw_cal.getClockTime();
      const CalibrationResult_& cal = last_calibration_result_;

      if (cal.success)
      {
        // Fragment-side calibration is independent of the precursor positive-magnitude
        // representability check — always apply it when calibration succeeds.
        Param fi_params = fi_params_original;
        fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
        effective_fragment_tol = cal.fragment_tolerance;

        if (!cal.extreme_bias)
        {
          // Precursor calibration is representable in the positive-magnitude schema —
          // apply the calibrated bounds.
          fi_params.setValue("precursor:mass_tolerance_lower", cal.cal_lower);
          fi_params.setValue("precursor:mass_tolerance_upper", cal.cal_upper);

          // Refresh algo-level member copies — OpenSearchModificationAnalysis reads them
          // via computeModMatchTolerance_() and would otherwise see stale pre-calibration
          // values. Restored to originals by restore_fi_params() below on return paths.
          precursor_mass_tolerance_lower_ = cal.cal_lower;
          precursor_mass_tolerance_upper_ = cal.cal_upper;
          effective_precursor_tol = std::max(cal.cal_lower, cal.cal_upper);

          OPENMS_LOG_INFO << "[ProSE] Calibration: shift=" << cal.precursor_shift
                          << " spread=" << cal.precursor_spread << " "
                          << precursor_mass_tolerance_unit_
                          << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                          << std::endl;
        }
        else
        {
          OPENMS_LOG_WARN << "[ProSE] Calibration: |shift|=" << std::abs(cal.precursor_shift)
                          << " > spread=" << cal.precursor_spread << " "
                          << precursor_mass_tolerance_unit_ << " - precursor calibration discarded. "
                          << "The true signed window ["
                          << (cal.precursor_shift - cal.precursor_spread) << ", "
                          << (cal.precursor_shift + cal.precursor_spread)
                          << "] lies entirely on one side of zero (not representable in the "
                          << "positive-magnitude schema without loosening). Fragment calibration "
                          << "still applied. Fix external calibration, or configure "
                          << "mass_tolerance_lower/_upper manually." << std::endl;
          // Precursor bounds unchanged; fragment_tolerance applied above.
        }

        fragment_index_.setParameters(fi_params);
        fi_params_modified = true;
      }
    }

    // Capture the mod-match tolerance reflecting the current (post-calibration)
    // member state. This fires unconditionally — even for closed searches that
    // won't invoke OpenSearchModificationAnalysis — so tests can observe whether
    // calibration wiring reaches the helper regardless of search mode. search()
    // itself does not run OSMA: the searchWithModificationAnalysis wrappers run it
    // on the PSMs returned here and read this field instead of re-computing (by
    // then restore_fi_params() has reset the members to user-configured values).
    last_mod_match_tolerance_used_ = computeModMatchTolerance_();

    // spectrum generators with the ion series of the FragmentIndex
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // storage for the PSMs of every spectrum (each grows on the thread that scores the spectrum)
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    bool open_search_mode = open_search;

    StopWatch sw_search; sw_search.start();
    scoreSpectraAgainstIndex_(spectra, fragment_index_, db, generators,
                              effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                              open_search_mode, annotated_hits, pool_stats,
                              "Scoring peptide models against spectra...", query_ptr);
    freeInParallel(query_spectra);

    // chimeric:enabled: search the co-isolated precursors in what the other PSMs of their spectrum leave of it. Their
    // hits and pool statistics are appended after those of the searched spectra; their residual spectra are built
    // only after the fragment index is released (below), and removed again after the post-processing.
    const Size n_searched_spectra = spectra.size();
    std::vector<CoIsolatedPSM_> coisolated;
    std::vector<Size> coisolated_parent;
    if (chimericSearch_())
    {
      std::vector<std::vector<AnnotatedHit_>> coisolated_hits;
      std::vector<CandidatePoolStats_> coisolated_stats;
      searchCoIsolated_(spectra, annotated_hits, fragment_index_, db, generators, effective_fragment_tol,
                        fragment_mass_tolerance_unit_ppm, open_search_mode, coisolated, coisolated_hits, coisolated_stats);
      annotated_hits.reserve(annotated_hits.size() + coisolated.size());
      pool_stats.reserve(pool_stats.size() + coisolated.size());
      coisolated_parent.reserve(coisolated.size());
      for (Size k = 0; k < coisolated.size(); ++k)
      {
        annotated_hits.push_back(std::move(coisolated_hits[k]));
        pool_stats.push_back(std::move(coisolated_stats[k]));
        coisolated_parent.push_back(coisolated[k].parent);
      }
    }

    // The PeptideIndexing that maps the hits to their proteins below.
    // The PeptideIndexer drops peptides whose termini do not match the configured
    // specificity, so it must agree with the search-time setting — otherwise
    // semi-specific / non-specific PSMs would be silently filtered out here.
    PeptideIndexing indexer;
    {
      Param param_pi = indexer.getParameters();
      // Use the effective decoy marker/position recorded in the context (the same
      // decoys that were searched), avoiding PeptideIndexing's own auto-detection.
      param_pi.setValue("decoy_string", ctx.decoy_string.empty() ? decoy_prefix_ : ctx.decoy_string);
      param_pi.setValue("decoy_string_position", ctx.decoy_is_prefix ? "prefix" : "suffix");
      param_pi.setValue("enzyme:name", enzyme_);
      param_pi.setValue("enzyme:specificity",
                        EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
      param_pi.setValue("missing_decoy_action", "silent");
      indexer.setParameters(param_pi);
    }
    // peptide:protein_mapping = index: the protein occurrences of the hits' candidates are read from the index
    // while it is still there (the index holds them all, PeptideIndexing searches the whole database for them)
    ProteinMapping_ protein_mapping;
    if (protein_mapping_from_index_)
    {
      std::vector<UInt32> candidates;
      for (const auto& hits : annotated_hits)
      {
        for (const auto& ah : hits) candidates.push_back(ah.index_peptide);
      }
      protein_mapping = buildProteinMapping_(fragment_index_, db, std::move(candidates), indexer.getParameters());
    }

    // M1: release the fragment index eagerly when the caller opted in (single-
    // use context). All downstream work (postProcessHits_, open-search mod
    // analysis, PeptideIndexing) is FI-independent, and the subsequent
    // Aho-Corasick pass inside PeptideIndexing::run() is the RSS high-water
    // mark of the whole search on large databases. Freeing here cuts
    // hundreds of MB of steady-state peak on human-proteome runs.
    // Not unconditional: external callers building ctx via prepareContext()
    // and calling search(spectra, ctx, ...) multiple times would break.
    if (ctx.release_fragment_index_after_scoring)
    {
      fragment_index_.clear();
    }

    takeBackThreadsOfEndedRead();
    if (!coisolated.empty())
    {
      spectra.reserve(spectra.size() + coisolated.size());
      for (CoIsolatedPSM_& psm : coisolated)
      {
        MSSpectrum residual = residualSpectrum_(spectra[psm.parent], psm.removed, psm.precursor);
        spectra.addSpectrum(std::move(residual));
        psm = CoIsolatedPSM_();
      }
      std::vector<CoIsolatedPSM_>().swap(coisolated);
    }

    startProgress(0, 1, "Post-processing PSMs...");
    // The fragment annotations of a hit are needed only for the output: when the rescoring (all but the best hit of a
    // spectrum) or the PSM FDR drops hits, they are made after that, for the hits that remain. Not if anything reads
    // them before (PeptDeep features) or keeps the hits of before (the pre-FDR copy). A search that fails in the protein
    // mapping makes them before it returns its error.
    // Only for hits of this search (identifications already in @p peptide_ids keep theirs), found again by their
    // scan_index (stored as unsigned int).
    // Not with chimeric:enabled: its residual spectra, which annotate the PSMs of co-isolated precursors, are gone by then.
    const bool defer_peak_annotations = (rescore_ || fdr_psm_ > 0.0) && !peptdeep_enable_ && !keep_pre_fdr_psms_
      && !chimericSearch_() && peptide_ids.empty() && spectra.size() <= std::numeric_limits<unsigned int>::max();
    ProSEAlgorithm::postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      searchTopHits_(),
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      effective_precursor_tol,
      effective_fragment_tol,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      "", // no database filename for in-memory search
      evidence_ptr,
      &coisolated_parent, // the residual spectra use the peak lists of their spectra
      defer_peak_annotations
      );
    // Release the hits and the additional peak list before PeptideIndexing.
    freeOnOwnerThreads(annotated_hits, [&pool_stats](SignedSize i) { return pool_stats[i].scoring_thread; }); // (where scored)
    freeInParallel(pool_stats);
    freeInParallel(evidence_spectra);
    endProgress();

    takeBackThreadsOfEndedRead();
    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);

    // chimeric:enabled: the PSMs of co-isolated precursors refer to the spectrum they were found in (scan_index and
    // spectrum reference), next to its own PSMs, and are told apart by chimeric_rank
    if (chimericSearch_())
    {
      Size n_coisolated_psms = 0;
      for (PeptideIdentification& pi : peptide_ids)
      {
        bool coisolated = false;
        if (pi.metaValueExists("scan_index"))
        {
          const Size scan = static_cast<unsigned int>(pi.getMetaValue("scan_index"));
          if (scan >= n_searched_spectra && scan - n_searched_spectra < coisolated_parent.size())
          {
            pi.setMetaValue("scan_index", static_cast<unsigned int>(coisolated_parent[scan - n_searched_spectra]));
            coisolated = true;
            ++n_coisolated_psms;
          }
        }
        for (PeptideHit& hit : pi.getHits()) { hit.setMetaValue(CHIMERIC_RANK, coisolated ? 2 : 1); }
      }
      std::stable_sort(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return static_cast<unsigned int>(a.getMetaValue("scan_index", 0u)) < static_cast<unsigned int>(b.getMetaValue("scan_index", 0u));
      });
      std::vector<MSSpectrum>& all = spectra.getSpectra();
      all.erase(all.begin() + static_cast<SignedSize>(n_searched_spectra), all.end());
      OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << n_coisolated_psms << " PSMs of co-isolated precursors." << std::endl;
    }
    sw_search.stop();
    last_run_stats_.seconds_search = sw_search.getClockTime();

    // map the hits to their proteins: from the index, or with PeptideIndexing (the same result)
    PeptideIndexing::ExitCodes indexer_exit = PeptideIndexing::ExitCodes::EXECUTION_OK;
    std::string fallback_reason;
    if (protein_mapping_from_index_
        && applyProteinMapping_(protein_mapping, db, indexer.getParameters(), protein_ids, peptide_ids, fallback_reason))
    {
      OPENMS_LOG_INFO << "[ProSE] Protein mapping from the fragment index: " << protein_mapping.sequences.size()
                      << " sequences, " << protein_ids[0].getHits().size() << " proteins." << std::endl;
    }
    else
    {
      if (protein_mapping_from_index_)
      {
        OPENMS_LOG_INFO << "[ProSE] Protein mapping with PeptideIndexing (" << fallback_reason << ")." << std::endl;
      }
      indexer_exit = indexer.run(db, protein_ids, peptide_ids);
    }
    protein_mapping = ProteinMapping_(); // released before the FDR

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      // (the hits returned with the error carry their fragment annotations, as before the deferral)
      if (defer_peak_annotations) annotateFragmentPeaks_(spectra, peptide_ids, effective_fragment_tol, fragment_mass_tolerance_unit_);
      restore_fi_params();
      if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
      {
        return ExitCodes::INPUT_FILE_EMPTY;
      }
      else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
      {
        return ExitCodes::UNEXPECTED_RESULT;
      }
      else
      {
        return ExitCodes::UNKNOWN_ERROR;
      }
    }

    // PSM-level FDR filtering. The context records whether decoys are present
    // (generated internally, or external decoys reused from the input FASTA).
    const bool has_decoys = ctx.have_decoys;

    // Optional per-run ion priors need the target/decoy labels and the native scores: after
    // PeptideIndexing, before FDR overwrites the scores.
    takeBackThreadsOfEndedRead();
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Capture pre-FDR stats now — BEFORE FDR, which drops pure-decoy hits
    // (FalseDiscoveryRate default add_decoy_peptides=false), may strip decoys
    // entirely, and overwrites each hit's HyperScore with its q-value.
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
    last_percolator_ = applyPsmFDR_(protein_ids, peptide_ids, has_decoys, last_run_stats_,
                                    keepPreRescoringPsms_() ? &pre_fdr_peptide_ids_ : nullptr, &inference_only_peptide_ids_);
    if (defer_peak_annotations)
    {
      takeBackThreadsOfEndedRead();
      // (counted as search time, as when postProcessHits_ made them)
      StopWatch sw_annotations;
      sw_annotations.start();
      annotateFragmentPeaks_(spectra, peptide_ids, effective_fragment_tol, fragment_mass_tolerance_unit_);
      // A file that Percolator did not rescore on its own may end with its PSMs of before the rescoring
      // (rescoreFilesTogether_()): they get theirs too.
      if (pool_rescoring_ && !last_percolator_.applied)
      {
        annotateFragmentPeaks_(spectra, pre_fdr_peptide_ids_, effective_fragment_tol, fragment_mass_tolerance_unit_);
      }
      sw_annotations.stop();
      last_run_stats_.seconds_search += sw_annotations.getClockTime();
    }

    restore_fi_params();

    takeBackThreadsOfEndedRead();
    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // Shared protein-FDR finalization for a COMPLETE protein set. Single source
  // of truth for both the file-based search() below and the ProSE TOPP tool's
  // single-file path (the two previously held byte-identical copies).
  // =====================================================================
  // static
  std::string ProSEAlgorithm::inferProteins(std::vector<ProteinIdentification>& protein_ids,
                                            const PeptideIdentificationList& peptide_ids,
                                            const PeptideIdentificationList& inference_only_psms)
  {
    if (protein_ids.empty()) { return ""; }
    ProteinIdentification& run = protein_ids[0];

    // Copies of the best hits with what the inference reads (sequence, charge, peptide evidences, target/decoy and
    // the scores), in the run: the inference switches their scores and filters them, while the PSMs stay as they are.
    // PSMs without hits are left out, as the score switch reads the first PSM.
    PeptideIdentificationList psms;
    psms.reserve(peptide_ids.size() + inference_only_psms.size());
    std::vector<std::string> keys;
    const auto add = [&](const PeptideIdentification& id)
    {
      if (id.getHits().empty()) { return; }
      const PeptideHit* best = &id.getHits()[0];
      for (const PeptideHit& hit : id.getHits())
      {
        if (id.isHigherScoreBetter() ? hit.getScore() > best->getScore() : hit.getScore() < best->getScore()) { best = &hit; }
      }
      PeptideHit hit(best->getScore(), best->getRank(), best->getCharge(), best->getSequence());
      hit.setPeptideEvidences(best->getPeptideEvidences());
      best->getKeys(keys);
      for (const std::string& key : keys)
      {
        if (key == Constants::UserParam::TARGET_DECOY || Scores::isKnownScoreType(key)) { hit.setMetaValue(key, best->getMetaValue(key)); }
      }
      PeptideIdentification& copy = psms.emplace_back();
      copy.setIdentifier(run.getIdentifier());
      copy.setScoreType(id.getScoreType());
      copy.setHigherScoreBetter(id.isHigherScoreBetter());
      copy.setHits({std::move(hit)});
    };
    for (const PeptideIdentification& id : peptide_ids)
    {
      if (id.getIdentifier() == run.getIdentifier()) { add(id); } // the PSMs of other runs are not the run's evidence
    }
    for (const PeptideIdentification& id : inference_only_psms) { add(id); }

    // The protein score: the best PEP of its PSMs, if they have one (Percolator's, also when its q-value is the main
    // score), as ProteomicsLFQ scores proteins; otherwise the best search engine score (e.g. 'ln(hyperscore)', also
    // after FDR:PSM made the q-value the main score). The q-values themselves rank the proteins coarsely: many PSMs
    // share the best q-value.
    // Every PSM must have the score: when Percolator rescored only some files of a pooled set, only their PSMs have a
    // PEP, and all of them are then scored by ProSE's own score, which all have ('ln(hyperscore)', the main score, or
    // a meta value after rescoring or FDR:PSM).
    std::string score_type;
    if (!psms.empty())
    {
      const IDScoreSwitcherAlgorithm switcher;
      const auto all_have = [&](IDScoreSwitcherAlgorithm::ScoreType type)
      {
        return std::all_of(psms.begin(), psms.end(), [&](const PeptideIdentification& id)
                           { return !switcher.findScoreType(id, type).score_name.empty(); });
      };
      if (all_have(IDScoreSwitcherAlgorithm::ScoreType::PEP)) { score_type = "PEP"; }
      // (after useSearchEngineScore_(), ProSE's score is the main score of the copies, which "RAW" finds)
      else if (useSearchEngineScore_(psms) || all_have(IDScoreSwitcherAlgorithm::ScoreType::RAW)) { score_type = "RAW"; }
    }
    BasicProteinInferenceAlgorithm bpia;
    Param bpia_param = bpia.getParameters();
    bpia_param.setValue("score_type", score_type);
    bpia.setParameters(bpia_param);
    bpia.run(psms, run);
    return score_type;
  }

  // static
  void ProSEAlgorithm::removeProteinsReferencedByNeither_(std::vector<ProteinIdentification>& protein_ids,
                                                         const PeptideIdentificationList& peptide_ids,
                                                         const PeptideIdentificationList& also_referencing)
  {
    // the accessions each run's PSMs of either list reference (as IDFilter::removeUnreferencedProteins(), by identifier)
    std::map<std::string, std::unordered_set<std::string>> referenced;
    for (const PeptideIdentificationList* list : {&peptide_ids, &also_referencing})
    {
      for (const PeptideIdentification& id : *list)
      {
        std::unordered_set<std::string>& accessions = referenced[id.getIdentifier()];
        for (const PeptideHit& hit : id.getHits())
        {
          for (const PeptideEvidence& evidence : hit.getPeptideEvidences()) { accessions.insert(evidence.getProteinAccession()); }
        }
      }
    }
    for (ProteinIdentification& run : protein_ids)
    {
      const std::unordered_set<std::string>& accessions = referenced[run.getIdentifier()];
      std::vector<ProteinHit>& hits = run.getHits();
      hits.erase(std::remove_if(hits.begin(), hits.end(), [&](const ProteinHit& hit) { return !accessions.contains(hit.getAccession()); }),
                 hits.end());
    }
  }

  // static
  void ProSEAlgorithm::removeUnreportedProteins(std::vector<ProteinIdentification>& protein_ids,
                                                const PeptideIdentificationList& peptide_ids,
                                                const PeptideIdentificationList& also_referencing)
  {
    if (also_referencing.empty()) { IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids); }
    else { removeProteinsReferencedByNeither_(protein_ids, peptide_ids, also_referencing); }
    for (ProteinIdentification& run : protein_ids)
    {
      IDFilter::updateProteinGroups(run.getIndistinguishableProteins(), run.getHits());
      IDFilter::updateProteinGroups(run.getProteinGroups(), run.getHits());
    }
  }

  // static
  bool ProSEAlgorithm::applyCompleteSetProteinFDR(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      const std::string& decoy_string,
      bool decoy_is_prefix,
      double protein_fdr,
      bool remove_decoys,
      const PeptideIdentificationList& inference_only_psms)
  {
    // Defensive: an empty decoy marker would make every protein match the prefix/suffix guard
    // below (hasPrefix(x, "") is always true) and corrupt picked-protein FDR. Callers gate on
    // have_decoys (so decoy_string is non-empty in practice), but this is a public static helper.
    if (decoy_string.empty())
    {
      OPENMS_LOG_WARN << "[ProSE] applyCompleteSetProteinFDR called with an empty decoy marker; "
                      << "skipping protein FDR (cannot identify decoys)." << std::endl;
      return false;
    }

    // Aggregate the best PSM score (PEP, or the search engine score) per peptide per protein over the reported PSMs and
    // the inference-only ones, then apply picked-protein FDR (Savitski et al. 2015) over the full target+decoy protein
    // set. The inference uses the best hit of each PSM, and so do the reported PSMs from here on.
    IDFilter::keepNBestHits(peptide_ids, 1);
    const std::string score_type = inferProteins(protein_ids, peptide_ids, inference_only_psms);

    const auto is_decoy = [&](const ProteinHit& ph)
    {
      return decoy_is_prefix ? StringUtils::hasPrefix(ph.getAccession(), decoy_string)
                             : StringUtils::hasSuffix(ph.getAccession(), decoy_string);
    };
    const auto& inferred = protein_ids[0].getHits();
    const Size decoy_proteins = static_cast<Size>(std::count_if(inferred.begin(), inferred.end(), is_decoy));
    const auto with_hits = [](const PeptideIdentification& id) { return !id.getHits().empty(); };
    const Size reported = static_cast<Size>(std::count_if(peptide_ids.begin(), peptide_ids.end(), with_hits));
    const Size inference_only = static_cast<Size>(std::count_if(inference_only_psms.begin(), inference_only_psms.end(), with_hits));
    OPENMS_LOG_INFO << "[ProSE] Protein inference on " << reported + inference_only << " PSMs (" << inference_only
                    << " of them up to FDR:protein_inference_PSM, not reported), protein score "
                    << (score_type.empty() ? "the main score" : score_type) << ": " << inferred.size() - decoy_proteins
                    << " target and " << decoy_proteins << " decoy proteins." << std::endl;

    // Picked-protein FDR needs identified decoy PROTEINS, not merely a decoy database: if PSM-level
    // filtering removed all decoy evidence, q-values would be estimated from targets only. The
    // callers gate on a DB-level / pre-inference "decoys present" check, which does not catch this
    // post-inference target-only edge -- guard it here.
    if (decoy_proteins == 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Protein FDR requested but no decoy proteins remain after inference "
                      << "(marker '" << decoy_string << "'); skipping picked-protein FDR to avoid "
                      << "target-only q-values." << std::endl;
      removeUnreportedProteins(protein_ids, peptide_ids); // the proteins only inference-only PSMs support
      return false;
    }

    FalseDiscoveryRate fdr;
    if (!remove_decoys)
    {
      // By default the picked-protein FDR drops the decoy proteins from the run; keep them, with q-values, so that the
      // filter below applies to them as to the targets (as mergeFileResults() does). The targets' q-values are the same.
      Param fdr_param = fdr.getParameters();
      fdr_param.setValue("add_decoy_proteins", "true");
      fdr.setParameters(fdr_param);
    }
    fdr.applyPickedProteinFDR(protein_ids[0], decoy_string, decoy_is_prefix);
    IDFilter::filterHitsByScore(protein_ids, protein_fdr);
    // The indistinguishable groups have q-values of their own (the picked FDR of the groups, which counts a group
    // once, however many targets it holds): a group whose proteins pass can still be above the cutoff. Filter the
    // groups at the cutoff too, and keep only the proteins of the groups that pass, so that every reported protein
    // and group is within it.
    std::vector<ProteinIdentification::ProteinGroup>& groups = protein_ids[0].getIndistinguishableProteins();
    if (!groups.empty())
    {
      IDFilter::filterGroupsByScore(groups, protein_fdr, protein_ids[0].isHigherScoreBetter());
      IDFilter::removeUngroupedProteins(groups, protein_ids[0].getHits());
    }

    // Decoy removal is the finalization step. removeDecoyHits strips decoy PSMs; decoy
    // proteins then fall out as unreferenced. Repair indistinguishable-protein and
    // protein-group references and drop peptide evidence pointing at removed proteins,
    // else idXML storage fails on dangling references.
    if (remove_decoys)
    {
      IDFilter::removeDecoyHits(peptide_ids);
      IDFilter::removeEmptyIdentifications(peptide_ids);
    }
    // Proteins (and groups) without a reported PSM go: the decoys once their PSMs are removed, and the proteins that
    // only inference-only PSMs support.
    removeUnreportedProteins(protein_ids, peptide_ids);
    // A PSM all of whose proteins failed the protein FDR goes too (as in ProteomicsLFQ): it would be reported without
    // any protein (no protein accession, start, end or flanking residues), so every reported PSM maps to a reported
    // protein. The proteins of the PSMs that stay are not affected.
    IDFilter::removeDanglingProteinReferences(peptide_ids, protein_ids, /*remove_peptides_without_reference=*/ true);
    IDFilter::removeEmptyIdentifications(peptide_ids);

    OPENMS_LOG_INFO << "[ProSE] Protein inference + picked-protein FDR: "
                    << protein_ids[0].getHits().size() << " proteins at "
                    << protein_fdr * 100 << "% FDR." << std::endl;
    return true;
  }

  // static
  bool ProSEAlgorithm::hasPSMs_(const PeptideIdentificationList& ids)
  {
    return std::any_of(ids.begin(), ids.end(), [](const PeptideIdentification& id) { return !id.getHits().empty(); });
  }

  // static
  void ProSEAlgorithm::tdcQValues_(PeptideIdentificationList& peptide_ids, bool plus_one)
  {
    // With plus_one: the best hit of each PSM (its score, oriented higher-better, and whether it is a decoy), which
    // FalseDiscoveryRate keeps, before it replaces the scores.
    std::vector<std::tuple<double, bool, Size>> best; // score, decoy, PSM
    if (plus_one)
    {
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty()) { continue; }
        id.sort();
        const PeptideHit& hit = id.getHits()[0];
        best.emplace_back(id.isHigherScoreBetter() ? hit.getScore() : -hit.getScore(), hit.isDecoy(), i);
      }
    }
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getParameters();
    fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
    fdr.setParameters(fdr_params);
    // Without a decoy PSM, FalseDiscoveryRate's D/T is 0 for every PSM, which would let any threshold pass PSMs
    // whose FDR nothing estimates. Use the +1-corrected estimate (D + 1) / T instead: 1/T for every PSM (the q-value
    // of the best of T targets, and of all below it), so that a threshold below 1/T reports no PSM, as an
    // unattainable protein FDR reports no protein.
    fdr.apply(peptide_ids); // keeps the best hit of each PSM, decoys included (add_decoy_peptides)
    Size targets = 0, decoys = 0;
    for (const PeptideIdentification& id : peptide_ids)
    {
      if (!id.getHits().empty()) { (id.getHits()[0].isDecoy() ? decoys : targets) += 1; }
    }
    if (decoys == 0 && targets > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] FDR:PSM: none of the " << targets << " PSMs matched a decoy; their q-value is the "
                      << "conservative estimate 1/" << targets << " = " << 1.0 / static_cast<double>(targets)
                      << " (no PSM passes a lower threshold)." << std::endl;
      if (!plus_one)
      {
        for (PeptideIdentification& id : peptide_ids)
        {
          for (PeptideHit& hit : id.getHits()) { hit.setScore(1.0 / static_cast<double>(targets)); }
        }
      }
    }
    if (plus_one)
    {
      // (D + 1) / T of the PSMs scoring at least as well (ties together), the minimum of it at and below each score
      // (Percolator's PosteriorEstimator::getQValues() without pi0), at most 1. Without a decoy, 1/T for every PSM.
      std::stable_sort(best.begin(), best.end(), [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
      std::vector<double> q(best.size());
      Size d = 0, t = 0;
      for (Size k = 0; k < best.size();)
      {
        Size end = k; // past the PSMs tied with PSM k
        for (; end < best.size() && std::get<0>(best[end]) == std::get<0>(best[k]); ++end) { (std::get<1>(best[end]) ? d : t) += 1; }
        std::fill(q.begin() + k, q.begin() + end, std::min(1.0, static_cast<double>(d + 1) / static_cast<double>(std::max<Size>(1, t))));
        k = end;
      }
      for (Size k = best.size(); k-- > 1;) { q[k - 1] = std::min(q[k - 1], q[k]); }
      for (Size k = 0; k < best.size(); ++k)
      {
        for (PeptideHit& hit : peptide_ids[std::get<2>(best[k])].getHits()) { hit.setScore(q[k]); }
      }
    }
  }

  void ProSEAlgorithm::annotatePsmQValues(PeptideIdentificationList& peptide_ids) const
  {
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getParameters();
    fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
    fdr.setParameters(fdr_params);

    // Group the spectra by the number of fragment charges their best hit was scored with. FalseDiscoveryRate
    // competes the best hit of each spectrum (stable score order, as sort() below).
    std::map<int, std::vector<Size>> groups;
    bool separate = fdr_psm_by_scored_charges_ && scoring_multiple_charges_;
    if (separate)
    {
      std::map<int, std::pair<bool, bool>> has_target_decoy;
      std::vector<Size> without_hits;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty())
        {
          without_hits.push_back(i);
          continue;
        }
        id.sort();
        const PeptideHit& best = id.getHits()[0];
        const int group = scoringMaxCharge_(best.getCharge());
        groups[group].push_back(i);
        auto& [has_target, has_decoy] = has_target_decoy[group];
        (best.isDecoy() ? has_decoy : has_target) = true;
      }
      separate = groups.size() > 1
                 && std::all_of(has_target_decoy.begin(), has_target_decoy.end(),
                                [](const auto& entry) { return entry.second.first && entry.second.second; });
      if (groups.size() > 1 && ! separate)
      {
        OPENMS_LOG_WARN << "[ProSE] FDR:PSM_groups: a group of PSMs scored with the same number of fragment charges has no "
                        << "target or no decoy PSM; computing q-values over all PSMs together." << std::endl;
      }
      if (separate && ! without_hits.empty()) // they carry no score, but get the q-value score type with the others
      {
        auto& first = groups.begin()->second;
        first.insert(first.end(), without_hits.begin(), without_hits.end());
      }
    }
    if (! separate)
    {
      tdcQValues_(peptide_ids);
      return;
    }
    for (const auto& [group, indices] : groups)
    {
      PeptideIdentificationList part;
      part.reserve(indices.size());
      for (Size i : indices) { part.push_back(std::move(peptide_ids[i])); }
      fdr.apply(part);
      for (Size k = 0; k < indices.size(); ++k) { peptide_ids[indices[k]] = std::move(part[k]); }
    }
  }

  namespace
  {
    // A helper thread that destroys hits, joined when this goes (see applyPsmFDR_()).
    struct HitRelease
    {
      std::thread thread;

      ~HitRelease()
      {
        if (thread.joinable()) { thread.join(); }
      }
    };

    // A helper thread that releases spectra, joined when this goes, after the release held before (see
    // releaseSpectraLater_())
    struct SpectraRelease
    {
      std::shared_ptr<void> before;
      std::thread thread;

      ~SpectraRelease()
      {
        if (thread.joinable()) { thread.join(); }
      }
    };
  }

  // static
  std::string ProSEAlgorithm::pinFeatureEnzyme(const std::string& enzyme, bool fixed_enzyme_termini)
  {
    // a ProteaseDB enzyme without a Percolator name keeps its name: PercolatorInfile then uses its cleavage rule
    const std::string name = pinEnzymeName(enzyme);
    return (name == "no_enzyme" ? enzyme : name) + (fixed_enzyme_termini ? "+termini" : "+keep_termini");
  }

  // static
  std::string ProSEAlgorithm::pinEnzymeName(const std::string& enzyme)
  {
    std::string name;
    for (const char c : enzyme)
    {
      if (c != ' ') name += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (name == "trypsin/p") return "trypsinp";
    if (name == "pepsina") return "pepsin";
    static const std::set<std::string> pin_enzymes = {"elastase", "pepsin", "proteinasek", "thermolysin", "chymotrypsin",
                                                      "lys-n", "lys-c", "arg-c", "asp-n", "glu-c", "trypsin", "trypsinp"};
    return pin_enzymes.contains(name) ? name : "no_enzyme";
  }

  // =====================================================================
  // rescore_reuse: cross-fitted Percolator models, trained once on PSMs of one or more files of a search and applied
  // to the files (static models: Fondrie & Noble, J Proteome Res 2020; joint models: mokapot, Fondrie & Noble, J
  // Proteome Res 2021; training subsampled over runs and scores relative to the decoys: PyProphet, Rosenberger et al.,
  // Nat Methods 2017).
  //
  // How the FDR stays valid: (1) No PSM is scored by a model that was trained on it, and no PSM's label enters the
  // transform of its own score. The spectra are split into two halves by a seeded hash of (file index, spectrum
  // index); one model is trained per half, and a PSM is scored by the model of the other half. All PSMs of a spectrum,
  // target and decoy, are in the same half, so Percolator's target-decoy competition (spectrum and experimental mass,
  // as without rescore_reuse) compares scores of the same model. A file that does not train ('static') is seen by
  // neither model. As in Percolator's cross-validation, the model that scores a PSM was trained on other spectra of
  // the same file (with 'joint' and 'own'): a PSM's label enters only models that score other PSMs. (2) The training
  // spectra are drawn by the hash, independently of labels and scores. (3) The two models' scores are put on one
  // scale, per model: relative to the mean and standard deviation of its scores of calibration decoys, decoys of its
  // own half that it was not trained on (5 % of the spectra of each half are kept out of the training for this);
  // these PSMs are scored by the other model. (4) q-values and PEPs are estimated per file by Percolator::score() on
  // the file's own PSMs; they are not pooled across files. (5) The domain-shift check of 'static' compares the feature
  // distributions of all PSMs of a file with those of the training PSMs, without their target/decoy labels, so the
  // choice between the static models and Percolator trained on the file does not depend on the labels it estimates
  // the FDR from.
  // =====================================================================
  namespace
  {
    UInt64 reuseMix(UInt64 x)
    {
      x += 0x9e3779b97f4a7c15ULL;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
      x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
      return x ^ (x >> 31);
    }

    // stream 1: the half of a spectrum (top bit); stream 2: its uniform draw for the training sample
    UInt64 reuseHash(UInt64 seed, UInt64 stream, Size file, Size spectrum)
    {
      return reuseMix(reuseMix(reuseMix(seed * 4 + stream) ^ static_cast<UInt64>(file)) ^ static_cast<UInt64>(spectrum));
    }

    int reuseHalf(UInt64 seed, Size file, Size spectrum)
    {
      return static_cast<int>(reuseHash(seed, 1, file, spectrum) >> 63);
    }

    double reuseUniform(UInt64 seed, Size file, Size spectrum)
    {
      return static_cast<double>(reuseHash(seed, 2, file, spectrum) >> 11) * 0x1.0p-53;
    }

    // The raw SVM score of a feature row: weights x features + bias (the last weight)
    double reuseRawScore(const std::vector<double>& weights, const double* features, Size n_features)
    {
      double score = weights.back();
      for (Size j = 0; j < n_features; ++j) score += weights[j] * features[j];
      return score;
    }

    // AUC of @p a against @p b: the probability that a value of a exceeds one of b (ties count half)
    double reuseAUC(std::vector<double> a, std::vector<double> b)
    {
      if (a.empty() || b.empty()) return 0.5;
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      double u = 0.0;
      Size below = 0, not_above = 0; // values of b below x, not above x
      for (const double x : a)
      {
        while (below < b.size() && b[below] < x) ++below;
        not_above = std::max(not_above, below);
        while (not_above < b.size() && b[not_above] <= x) ++not_above;
        u += static_cast<double>(below) + 0.5 * static_cast<double>(not_above - below);
      }
      return u / (static_cast<double>(a.size()) * static_cast<double>(b.size()));
    }

    // The Percolator parameters of ProSE's rescoring: @p training (parameters 'rescoring:...', 'rescore_solver'), then
    // @p extra (rescore_method 'reset': the RESET parameters)
    Param rescoringPercolatorParameters(Param pp, int threads, const ProSEAlgorithm::PercolatorTraining& training,
                                        const Param& extra)
    {
      pp.setValue("train_best_positive", "true");
      pp.setValue("post_processing_tdc", "true");
      pp.setValue("subset_max_train", training.subset_max_train);  // 0 (the default): train on all PSMs
      pp.setValue("subset_mode", training.subset_mode);
      pp.setValue("subset_curve_start", training.curve_start);
      pp.setValue("subset_curve_min_gain", training.curve_min_gain);
      pp.setValue("num_iterations", training.iterations);
      pp.setValue("tron_epsilon", training.tron_epsilon);
      pp.setValue("warm_start", training.warm_start ? "true" : "false");
      pp.setValue("single_cneg", training.single_cneg ? "true" : "false");
      pp.setValue("pep_method", "nonparametric");  // as the percolator executable
      // cross-validation trains 9 models at a time (3 folds x 3 Cneg values) and sorts on more threads (up to 12); the
      // results do not depend on the threads
      pp.setValue("num_threads", threads <= 0 ? 3 : std::min(threads, 128));
      pp.setValue("solver", training.solver);
      pp.setValue("use_pi0", training.use_pi0);
      pp.setValue("calc_mass_feature", training.calc_mass_feature ? "true" : "false");
      if (training.peptide_level) pp.setValue("peptide_level", "true");
      for (auto it = extra.begin(); it != extra.end(); ++it)
      {
        pp.setValue(it.getName(), it->value);
      }
      return pp;
    }

    // The enzyme name of the PIN features of the rescoring: Percolator's name of the search's enzyme, with
    // PercolatorInfile's "+termini" modifier for 'rescoring:fixed_enzyme_termini'.
    std::string rescoringEnzyme(const std::string& enzyme, const ProSEAlgorithm::PercolatorTraining& training)
    {
      return ProSEAlgorithm::pinFeatureEnzyme(enzyme, training.fixed_enzyme_termini);
    }

    // The standard PIN features, the search engine's own features, and the bookkeeping columns.
    StringList rescoringFeatureSet(const ProteinIdentification::SearchParameters& sp, int min_charge, int max_charge)
    {
      StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
      if (sp.metaValueExists("extra_features"))
      {
        StringList extra = ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString());
        feature_set.insert(feature_set.end(), extra.begin(), extra.end());
      }
      else
      {
        feature_set.push_back("score");
      }
      feature_set.push_back("Peptide");
      feature_set.push_back("Proteins");
      return feature_set;
    }

    // Training PSMs kept for the domain-shift check, and PSMs of a file compared with them
    constexpr Size REUSE_SHIFT_PSMS = 20000;
    // The spectra of a half whose uniform draw is at least this are calibration spectra: never trained on
    constexpr double REUSE_CALIBRATION_FROM = 0.95;
    // Calibration decoys used per model (every k-th)
    constexpr Size REUSE_CALIBRATION_DECOYS = 200000;
  }

  bool ProSEAlgorithm::reuseTrainingDraw(UInt64 seed, Size file, Size spectrum, double fraction)
  {
    return reuseUniform(seed, file, spectrum) < REUSE_CALIBRATION_FROM * fraction;
  }

  bool ProSEAlgorithm::reuseCalibrationDraw(UInt64 seed, Size file, Size spectrum, double fraction)
  {
    const double u = reuseUniform(seed, file, spectrum);
    return u >= REUSE_CALIBRATION_FROM && (u - REUSE_CALIBRATION_FROM) / (1.0 - REUSE_CALIBRATION_FROM) < fraction;
  }

  struct ProSEAlgorithm::ReuseModel_
  {
    std::string mode;                    ///< rescore_reuse:mode of the training
    UInt64 seed = 0;                     ///< seed of the random halves and training spectra (Percolator's)
    StringList features;                 ///< numeric features of the models
    int min_charge = 0;                  ///< charge range of the charge features
    int max_charge = 0;
    std::array<PercolatorModel, 2> models; ///< models[h]: trained on half h
    /// mean and standard deviation of models[h]'s scores of the calibration decoys of half h
    std::array<double, 2> center{0.0, 0.0};
    std::array<double, 2> scale{1.0, 1.0};
    std::array<Size, 2> training_psms{0, 0};
    std::vector<Size> training_files;    ///< indices of the files that trained
    std::vector<std::vector<double>> sample; ///< training PSMs (a sample, targets and decoys), for the domain-shift check
    std::vector<Size> used_features;     ///< features with a weight in a model, for the domain-shift check
    bool used_default_direction = false; ///< Percolator's sanity check replaced a learned model

    std::string describe() const
    {
      std::string files;
      for (Size f : training_files) files += (files.empty() ? "" : ",") + StringUtils::toStr(f + 1);
      return mode + " (cross-fitted models trained on " + StringUtils::toStr(training_psms[0]) + " + "
             + StringUtils::toStr(training_psms[1]) + " PSMs of file(s) " + files + ")";
    }
  };

  struct ProSEAlgorithm::ReuseScoring_
  {
    /// what the scorer did (written by it)
    struct Note
    {
      std::string model;
      bool used_default_direction = false;
      double shift = -1.0;
      std::string fallback; ///< why Percolator is to be trained on the run instead (then the scorer throws)
    };
    Percolator::PSMScorer scorer;
    int min_charge = 0;
    int max_charge = 0;
    std::shared_ptr<Note> note;
    std::string mode;
    UInt64 seed = 0;
    Size max_train = 0;
  };

  std::shared_ptr<const ProSEAlgorithm::ReuseModel_> ProSEAlgorithm::trainReuseModel_(const std::vector<ReuseRun_>& runs,
                                                                                       int threads) const
  {
    auto model = std::make_shared<ReuseModel_>();
    model->mode = reuse_mode_;
    {
      // The random halves and training spectra follow Percolator's seed (as its cross-validation folds do).
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, percolator_training_, Param()));
      model->seed = static_cast<UInt64>(static_cast<int>(perc.getParameters().getValue("seed")));
    }
    const UInt64 seed = model->seed;
    // The runs that rescoreWithPercolator_() would rescore: a run, at least 100 PSMs, decoys.
    std::vector<const ReuseRun_*> used;
    int min_charge = 10;
    int max_charge = 0;
    for (const ReuseRun_& run : runs)
    {
      if (run.protein_ids->empty() || run.peptide_ids->size() < 100) continue;
      bool has_decoys = false;
      int run_min_charge = 10;
      int run_max_charge = 0;
      for (const PeptideIdentification& pid : *run.peptide_ids)
      {
        for (const PeptideHit& hit : pid.getHits())
        {
          has_decoys = has_decoys || hit.isDecoy();
          run_min_charge = std::min(run_min_charge, hit.getCharge());
          run_max_charge = std::max(run_max_charge, hit.getCharge());
        }
      }
      if (!has_decoys) continue;
      min_charge = std::min(min_charge, run_min_charge);
      max_charge = std::max(max_charge, run_max_charge);
      used.push_back(&run);
    }
    if (used.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "No run with at least 100 PSMs and decoy PSMs to train the models on.");
    }
    model->min_charge = min_charge;
    model->max_charge = max_charge;
    const ProteinIdentification::SearchParameters& sp = used.front()->protein_ids->front().getSearchParameters();
    const StringList feature_set = rescoringFeatureSet(sp, min_charge, max_charge);
    const std::string enzyme = rescoringEnzyme(sp.digestion_enzyme.getName(), percolator_training_);

    // The training PSMs of a half: max_train shared equally by the files, a file with fewer PSMs gives all of them (the
    // hits of a spectrum count, as an estimate of its feature rows).
    const Size n_runs = used.size();
    std::vector<std::array<double, 2>> psms(n_runs, {0.0, 0.0});
    for (Size k = 0; k < n_runs; ++k)
    {
      const PeptideIdentificationList& ids = *used[k]->peptide_ids;
      for (Size s = 0; s < ids.size(); ++s)
      {
        psms[k][reuseHalf(seed, used[k]->file_index, s)] += static_cast<double>(ids[s].getHits().size());
      }
    }
    std::vector<std::array<double, 2>> fraction(n_runs, {1.0, 1.0});
    if (reuse_max_train_ > 0)
    {
      for (int h = 0; h < 2; ++h)
      {
        std::vector<Size> order(n_runs);
        std::iota(order.begin(), order.end(), Size(0));
        std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return psms[a][h] < psms[b][h]; });
        double budget = static_cast<double>(reuse_max_train_);
        Size left = n_runs;
        for (const Size k : order)
        {
          const double quota = budget / static_cast<double>(left--);
          const double trainable = psms[k][h] * REUSE_CALIBRATION_FROM;
          if (trainable > quota) { fraction[k][h] = quota / trainable; }
          budget -= std::min(trainable, quota);
        }
      }
    }
    // The calibration spectra of a half: a share of those with u >= REUSE_CALIBRATION_FROM, for about
    // 2 x REUSE_CALIBRATION_DECOYS PSMs per half.
    std::array<double, 2> calibration_fraction{1.0, 1.0};
    for (int h = 0; h < 2; ++h)
    {
      double half_psms = 0.0;
      for (Size k = 0; k < n_runs; ++k) { half_psms += psms[k][h]; }
      const double calibration_psms = half_psms * (1.0 - REUSE_CALIBRATION_FROM);
      if (calibration_psms > 2.0 * REUSE_CALIBRATION_DECOYS) { calibration_fraction[h] = 2.0 * REUSE_CALIBRATION_DECOYS / calibration_psms; }
    }
    std::array<RescoreInput, 2> train;
    StopWatch sw_rows;
    sw_rows.start();
    std::array<std::vector<std::vector<double>>, 2> calibration; // calibration decoys of each half
    for (Size k = 0; k < n_runs; ++k)
    {
      const Size file = used[k]->file_index;
      // Only the drawn spectra get feature rows: a copy of them (spectra[i]: the index of the copy's spectrum i).
      const PeptideIdentificationList& ids = *used[k]->peptide_ids;
      std::vector<Size> spectra;
      for (Size s = 0; s < ids.size(); ++s)
      {
        const int h = reuseHalf(seed, file, s);
        // (disjoint: a calibration spectrum is never a training spectrum)
        const bool drawn = reuseCalibrationDraw(seed, file, s, calibration_fraction[h])
                           || reuseTrainingDraw(seed, file, s, fraction[k][h]);
        if (drawn && !ids[s].getHits().empty()) { spectra.push_back(s); }
      }
      PeptideIdentificationList drawn_ids;
      drawn_ids.resize(spectra.size());
#pragma omp parallel for schedule(dynamic, 1024) num_threads(std::max(1, threads))
      for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i) { drawn_ids[i] = ids[spectra[i]]; }
      if (drawn_ids.empty()) { continue; }
      RescoreInput ri;
      Percolator::PSMRowLocations rows;
      Percolator::buildPSMRescoreInput(drawn_ids, feature_set, enzyme, min_charge, max_charge, ri, rows,
                                       percolator_training_.calc_mass_feature);
      if (model->features.empty()) { model->features = ri.feature_names; }
      else if (ri.feature_names != model->features)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The files of the search have different Percolator features.", "");
      }
      // (the rows of buildPSMRescoreInput are dense: the training rows are appended to dense matrices as well)
      const Size n_feat = ri.featureCount();
      for (Size r = 0; r < ri.rowCount(); ++r)
      {
        const Size spectrum = spectra[rows[r].first];
        const int h = reuseHalf(seed, file, spectrum);
        if (reuseUniform(seed, file, spectrum) >= REUSE_CALIBRATION_FROM)
        {
          if (ri.is_decoy[r]) { calibration[h].emplace_back(ri.row(r), ri.row(r) + n_feat); }
          continue;
        }
        RescoreInput& t = train[h];
        t.feature_matrix.insert(t.feature_matrix.end(), ri.row(r), ri.row(r) + n_feat);
        t.is_decoy.push_back(ri.is_decoy[r]);
        // The spectrum (not its scan number, which may be missing): Percolator's cross-validation folds and its best
        // PSM per spectrum go by file and scan number.
        t.scan_numbers.push_back(static_cast<int>(spectrum));
        t.spec_file_numbers.push_back(static_cast<int>(k));
        t.exp_masses.push_back(ri.exp_masses[r]);
        t.calc_masses.push_back(ri.calc_masses[r]);
      }
      model->training_files.push_back(file);
    }
    for (int h = 0; h < 2; ++h)
    {
      train[h].feature_names = model->features;
      model->training_psms[h] = train[h].rowCount();
      if (train[h].rowCount() == 0)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "No training PSMs in a half of the spectra.");
      }
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, percolator_training_, Param()));
      if (h == 0)
      {
        sw_rows.stop();
        OPENMS_LOG_INFO << "[ProSE] rescore_reuse: training PSMs of " << n_runs << " file(s) drawn in " << sw_rows.getClockTime()
                        << " s" << std::endl;
      }
      StopWatch sw_train;
      sw_train.start();
      model->models[h] = perc.train(train[h]);
      sw_train.stop();
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: model " << h + 1 << " trained on " << train[h].rowCount() << " PSMs in "
                      << sw_train.getClockTime() << " s" << std::endl;
      model->used_default_direction = model->used_default_direction || perc.usedDefaultDirection();
    }
    // The scale of each model: its scores of the calibration decoys of its half (not trained on; scored by the other model).
    for (int h = 0; h < 2; ++h)
    {
      const std::vector<std::vector<double>>& decoys = calibration[h];
      const Size step = std::max<Size>(1, (decoys.size() + REUSE_CALIBRATION_DECOYS - 1) / REUSE_CALIBRATION_DECOYS);
      double sum = 0.0, sum2 = 0.0;
      Size n = 0;
      for (Size r = 0; r < decoys.size(); r += step)
      {
        const double score = reuseRawScore(model->models[h].weights, decoys[r].data(), decoys[r].size());
        sum += score;
        sum2 += score * score;
        ++n;
      }
      if (n > 1)
      {
        const double mean = sum / static_cast<double>(n);
        const double sd = std::sqrt(std::max(0.0, (sum2 - sum * mean) / static_cast<double>(n - 1)));
        model->center[h] = mean;
        model->scale[h] = (sd > 0.0 && std::isfinite(sd)) ? sd : 1.0;
      }
    }
    // The domain-shift check: features with a weight, and a sample of the training PSMs (every k-th, without labels).
    for (Size j = 0; j < model->features.size(); ++j)
    {
      if (model->models[0].weights[j] != 0.0 || model->models[1].weights[j] != 0.0) { model->used_features.push_back(j); }
    }
    const Size n_train = train[0].rowCount() + train[1].rowCount();
    const Size step = std::max<Size>(1, (n_train + REUSE_SHIFT_PSMS - 1) / REUSE_SHIFT_PSMS);
    Size d = 0;
    for (const RescoreInput& t : train)
    {
      for (Size r = 0; r < t.rowCount(); ++r)
      {
        if (d++ % step == 0) { model->sample.emplace_back(t.row(r), t.row(r) + t.featureCount()); }
      }
    }
    return model;
  }

  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithReuseModel_(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      Size file_index,
      const ReuseModel_& model,
      bool check_shift,
      int threads,
      std::vector<std::vector<PeptideHit>>* dropped) const
  {
    ReuseScoring_ reuse;
    reuse.min_charge = model.min_charge;
    reuse.max_charge = model.max_charge;
    reuse.note = std::make_shared<ReuseScoring_::Note>();
    reuse.mode = model.mode;
    reuse.seed = model.seed;
    reuse.max_train = reuse_max_train_;
    const UInt64 seed = model.seed;
    const double max_shift = reuse_max_shift_;
    std::shared_ptr<ReuseScoring_::Note> note = reuse.note;
    // the reused models are Percolator models with the training options (RESET's parameters do not apply)
    const PercolatorTraining training = percolator_training_;
    reuse.scorer = [&model, note, seed, file_index, check_shift, max_shift, threads, training](
      const RescoreInput& ri, const Percolator::PSMRowLocations& rows) -> RescoreOutput
    {
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, training, Param()));
      // Percolator trained on the run, as without rescore_reuse: rescoreWithReuseModel_() rescores the run again without
      // the models (with the run's own feature set)
      const auto trained_on_run = [&](const std::string& why) -> RescoreOutput
      {
        note->fallback = why;
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "rescore_reuse: Percolator is to be trained on the run (" + why + ")", "");
      };
      if (ri.feature_names != model.features)
      {
        return trained_on_run("its Percolator features differ from those of the models");
      }
      const Size n = ri.rowCount();
      const Size n_feat = ri.featureCount();
      if (check_shift)
      {
        // all PSMs, without their labels (see above)
        const Size step = std::max<Size>(1, (n + REUSE_SHIFT_PSMS - 1) / REUSE_SHIFT_PSMS);
        double shift = 0.0;
        std::string worst;
        for (const Size j : model.used_features)
        {
          std::vector<double> a, b;
          for (Size r = 0; r < n; r += step) { a.push_back(ri.row(r)[j]); }
          for (const std::vector<double>& row : model.sample) { b.push_back(row[j]); }
          const double s = std::fabs(2.0 * reuseAUC(std::move(a), std::move(b)) - 1.0);
          if (s > shift) { shift = s; worst = model.features[j]; }
        }
        note->shift = shift;
        OPENMS_LOG_INFO << "[ProSE] rescore_reuse: domain shift of the PSMs " << shift << " (largest in feature '" << worst
                        << "'; limit " << max_shift << ")" << std::endl;
        if (shift > max_shift)
        {
          return trained_on_run("domain shift " + StringUtils::toStr(shift) + " in feature '" + worst + "'");
        }
      }
      // The cross-fitted score of every row: the model of the other half, on its scale.
      RescoreInput one;
      one.feature_names = {"rescore_reuse_score"};
      one.features.assign(n, std::vector<double>(1, 0.0));
#pragma omp parallel for schedule(static) num_threads(std::max(1, threads))
      for (SignedSize r = 0; r < static_cast<SignedSize>(n); ++r)
      {
        const int m = 1 - reuseHalf(seed, file_index, rows[r].first);
        one.features[r][0] = (reuseRawScore(model.models[m].weights, ri.row(r), n_feat) - model.center[m]) / model.scale[m];
      }
      one.is_decoy = ri.is_decoy;
      {
        double sum_t = 0.0, sum_d = 0.0;
        Size n_d = 0;
        for (Size r = 0; r < n; ++r) { (ri.is_decoy[r] ? sum_d : sum_t) += one.features[r][0]; n_d += ri.is_decoy[r] ? 1 : 0; }
        OPENMS_LOG_DEBUG << "[ProSE] rescore_reuse: " << n << " rows (" << n_d << " decoys) of file " << file_index + 1
                         << ", mean cross-fitted score of the targets " << sum_t / std::max<double>(1.0, double(n - n_d))
                         << ", of the decoys " << sum_d / std::max<double>(1.0, double(n_d)) << std::endl;
      }
      // The spectrum: best PSM per spectrum (target-decoy competition) as in training
      one.scan_numbers.resize(n);
      for (Size r = 0; r < n; ++r) { one.scan_numbers[r] = static_cast<int>(rows[r].first); }
      one.exp_masses = ri.exp_masses;
      one.calc_masses = ri.calc_masses;
      PercolatorModel unit;
      unit.feature_names = one.feature_names;
      unit.weights = {1.0, 0.0};
      try
      {
        RescoreOutput out = perc.score(one, unit);
        note->model = model.describe();
        note->used_default_direction = model.used_default_direction;
        return out;
      }
      catch (const Exception::BaseException& e)
      {
        return trained_on_run(std::string("scoring with the models failed: ") + e.what());
      }
    };
    PercolatorRescoring outcome = rescoreWithPercolator_(protein_ids, peptide_ids, threads, dropped, percolator_training_,
                                                         Param(), &reuse);
    if (note->fallback.empty())
    {
      outcome.model = note->model;
      return outcome;
    }
    // The run is rescored as without rescore_reuse (rescorePSMs() left the PSMs unchanged when the scorer threw).
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse: Percolator is trained on the file: " << note->fallback << std::endl;
    outcome = rescoreWithPercolator_(protein_ids, peptide_ids, threads, dropped, percolator_training_, Param());
    if (outcome.applied)
    {
      outcome.model = "trained on the file (" + note->fallback + ")";
      for (ProteinIdentification& run : protein_ids)
      {
        ProteinIdentification::SearchParameters run_sp = run.getSearchParameters();
        run_sp.setMetaValue("Percolator:rescore_reuse_mode", reuse.mode);
        run_sp.setMetaValue("Percolator:rescore_reuse_model", outcome.model);
        if (note->shift >= 0.0) { run_sp.setMetaValue("Percolator:rescore_reuse_shift", note->shift); }
        run.setSearchParameters(run_sp);
      }
    }
    return outcome;
  }

  // =====================================================================
  // In-process Percolator rescoring of one search run (parameter 'rescore'). Single
  // source of truth for the ProSE TOPP tool and other callers (e.g. the
  // WebAssembly bindings); it was a copy in the tool before.
  // =====================================================================
  // static
  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithPercolator(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      int threads,
      Size max_train_psms)
  {
    PercolatorTraining training;
    training.subset_max_train = static_cast<int>(std::min<Size>(max_train_psms, std::numeric_limits<int>::max()));
    return rescoreWithPercolator_(protein_ids, peptide_ids, threads, nullptr, training, Param());
  }

  // static
  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithPercolator_(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      int threads,
      std::vector<std::vector<PeptideHit>>* dropped,
      const PercolatorTraining& training,
      const Param& extra,
      const ReuseScoring_* reuse)
  {
    PercolatorRescoring outcome;
    if (protein_ids.empty())
    {
      outcome.skipped_reason = "No search run (protein identification) to rescore the PSMs of.";
      return outcome;
    }
    if (peptide_ids.size() < 100)
    {
      outcome.skipped_reason = "Only " + StringUtils::toStr(peptide_ids.size())
                               + " PSMs, Percolator rescoring needs at least 100.";
      return outcome;
    }

    bool has_decoys = false;
    int min_charge = 10;
    int max_charge = 0;
    for (const PeptideIdentification& pid : peptide_ids)
    {
      for (const PeptideHit& hit : pid.getHits())
      {
        has_decoys = has_decoys || hit.isDecoy();
        min_charge = std::min(min_charge, hit.getCharge());
        max_charge = std::max(max_charge, hit.getCharge());
      }
    }
    if (!has_decoys)
    {
      // The search's decoy mode is not known here: name both causes.
      outcome.skipped_reason = "No decoy PSMs for the target/decoy competition of Percolator. A search without "
                               "decoys needs decoys 'auto' or 'generate', or a FASTA with decoy proteins; with "
                               "decoys, no spectrum of this run matched one (e.g. a small run).";
      return outcome;
    }

    // The charge features of the models with rescore_reuse (a charge outside their range has none)
    if (reuse != nullptr)
    {
      min_charge = reuse->min_charge;
      max_charge = reuse->max_charge;
    }
    const ProteinIdentification::SearchParameters& sp = protein_ids.front().getSearchParameters();
    const StringList feature_set = rescoringFeatureSet(sp, min_charge, max_charge);

    Percolator perc;
    Param pp = rescoringPercolatorParameters(perc.getDefaults(), threads, training, extra);
    perc.setParameters(pp);
    if (training.subset_max_train > 0 && peptide_ids.size() > static_cast<Size>(training.subset_max_train))
    {
      OPENMS_LOG_INFO << "[ProSE] Percolator trains on " << training.subset_max_train << " of the " << peptide_ids.size()
                      << " PSMs (Search:percolator:max_train_psms) and scores all of them." << std::endl;
    }

    try
    {
      // (only the best hit of each spectrum is kept below: the others need no Percolator meta values)
      if (reuse == nullptr)
      {
        perc.rescorePSMs(peptide_ids, feature_set, rescoringEnzyme(sp.digestion_enzyme.getName(), training),
                         min_charge, max_charge, "q-value", /*best_hit_only=*/true);
      }
      else
      {
        perc.rescorePSMsWithScorer(peptide_ids, feature_set, rescoringEnzyme(sp.digestion_enzyme.getName(), training),
                                   min_charge, max_charge, "q-value", reuse->scorer, /*best_hit_only=*/true);
      }
    }
    catch (const Exception::BaseException& e)
    {
      outcome.skipped_reason = std::string("Percolator failed: ") + e.what();
      return outcome;
    }
    // The best hit per spectrum, as IDFilter::keepNBestHits(peptide_ids, 1), with the spectra in parallel: releasing
    // the other hits (and their meta values) takes most of the time. With the threads of the rescoring, as Percolator's
    // own loops (at most 9 threads, ProSE's former cap on Percolator's threads, released them 1.4 times slower at 16).
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
#else
    const int max_threads = 1;
#endif
    [[maybe_unused]] const int trim_threads = std::max(1, std::min(threads <= 0 ? 3 : threads, max_threads));
    // With @p dropped, the other hits are moved there (to be destroyed later, see applyPsmFDR_()) if memory allows,
    // otherwise destroyed here.
    if (dropped != nullptr)
    {
      try { dropped->resize(peptide_ids.size()); }
      catch (...) { dropped = nullptr; }
    }
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 256) num_threads(trim_threads)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      try
      {
        PeptideIdentification& pid = peptide_ids[i];
        pid.sort();
        std::vector<PeptideHit>& hits = pid.getHits();
        if (hits.size() > 1)
        {
          if (dropped != nullptr)
          {
            // moves no hit if the allocation fails (then they are destroyed below)
            try { (*dropped)[i].assign(std::make_move_iterator(hits.begin() + 1), std::make_move_iterator(hits.end())); }
            catch (...) {}
          }
          hits.resize(1);
        }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
    if (reuse != nullptr ? reuse->note->used_default_direction : perc.usedDefaultDirection())
    {
      outcome.warning = "Percolator's sanity check rejected the learned model (on the held-out folds it found no target "
                        "PSM at q < " + StringUtils::toStr(static_cast<double>(pp.getValue("test_fdr")))
                        + ", or fewer than the best single feature; typical for blank runs and runs with few confident "
                          "PSMs), so the PSMs are scored with that feature alone.";
    }

    markRescoredByPercolator_(protein_ids, pp, training, &perc, reuse);
    outcome.applied = true;
    return outcome;
  }

  // static
  void ProSEAlgorithm::markRescoredByPercolator_(std::vector<ProteinIdentification>& protein_ids, const Param& pp,
                                                 const PercolatorTraining& training, const Percolator* perc,
                                                 const ReuseScoring_* reuse)
  {
    // Mark the runs as post-processed by Percolator, as PercolatorAdapter does.
    for (ProteinIdentification& run : protein_ids)
    {
      run.setSearchEngine("Percolator");
      run.setSearchEngineVersion("3.08-vendored");
      run.setMetaValue("percolator", "ProSE");
      ProteinIdentification::SearchParameters run_sp = run.getSearchParameters();
      run_sp.setMetaValue("Percolator:peptide_level_fdrs", "false");
      run_sp.setMetaValue("Percolator:protein_level_fdrs", "false");
      run_sp.setMetaValue("Percolator:testFDR", static_cast<double>(pp.getValue("test_fdr")));
      run_sp.setMetaValue("Percolator:trainFDR", static_cast<double>(pp.getValue("train_fdr")));
      run_sp.setMetaValue("Percolator:maxiter", static_cast<int>(pp.getValue("num_iterations")));
      run_sp.setMetaValue("Percolator:subset_max_train", static_cast<int>(pp.getValue("subset_max_train")));
      // (only with a cap and an out-of-fold mode, or the learning curve: the output without a cap, or in_sample, stays
      // as it was)
      const std::string subset_mode = pp.getValue("subset_mode").toString();
      if (subset_mode == "learning_curve" ||
          (static_cast<int>(pp.getValue("subset_max_train")) > 0 && subset_mode != "in_sample"))
      {
        run_sp.setMetaValue("Percolator:subset_mode", subset_mode);
      }
      else run_sp.removeMetaValue("Percolator:subset_mode");
      // the learning curve: "validation=<rows>; <sample rows>:<validation targets per fold>; ...; chosen=<rows per fold>"
      if (subset_mode == "learning_curve" && perc != nullptr && !perc->getCappedTrainingInfo().curve_sizes.empty())
      {
        const Percolator::CappedTrainingInfo& capped = perc->getCappedTrainingInfo();
        std::string curve = "validation=" + std::to_string(capped.validation_rows);
        for (size_t k = 0; k < capped.curve_sizes.size(); ++k)
        {
          curve += "; " + std::to_string(capped.curve_sizes[k]) + ":";
          for (size_t f = 0; f < capped.curve_counts[k].size(); ++f)
          {
            curve += (f ? "/" : "") + std::to_string(capped.curve_counts[k][f]);
          }
        }
        curve += "; chosen=";
        for (size_t f = 0; f < capped.chosen_size.size(); ++f) curve += (f ? "/" : "") + std::to_string(capped.chosen_size[f]);
        run_sp.setMetaValue("Percolator:subset_curve", curve);
      }
      else run_sp.removeMetaValue("Percolator:subset_curve");
      run_sp.setMetaValue("Percolator:seed", static_cast<int>(pp.getValue("seed")));
      run_sp.setMetaValue("Percolator:post_processing_tdc", "true");
      run_sp.setMetaValue("Percolator:train_best_positive", "true");
      if (training.solver != "tron") run_sp.setMetaValue("Percolator:solver", training.solver);
      else run_sp.removeMetaValue("Percolator:solver");
      if (training.peptide_level && reuse == nullptr && perc != nullptr && !perc->getResetOutput())
      {
        run_sp.setMetaValue("Percolator:peptide_q_values", "percolator_peptide_q_value, percolator_peptide_pep");
      }
      else
      {
        run_sp.removeMetaValue("Percolator:peptide_q_values");
      }
      if (const auto rr = perc != nullptr ? perc->getResetOutput() : nullptr)
      {
        run_sp.setMetaValue("Percolator:algorithm", "reset");
        run_sp.setMetaValue("RESET:salt", StringUtils::toStr(rr->salt));
        run_sp.setMetaValue("RESET:protocol", "PSM-only");
        run_sp.setMetaValue("RESET:learner", rr->learner);
        if (!rr->learner_settings.empty()) run_sp.setMetaValue("RESET:learner_settings", rr->learner_settings);
        run_sp.setMetaValue("RESET:fdr_level", "peptide (stem): RESET:peptide_q_value; PSM q-values heuristic");
        run_sp.setMetaValue("RESET:target_stems_1pct", static_cast<int>(rr->discoveries_1pct));
        if (pp.getValue("reset:decoys").toString() == "two")
        {
          run_sp.setMetaValue("RESET:decoys", "two (second set: " + pp.getValue("reset:second_decoy_prefix").toString() + ")");
          run_sp.setMetaValue("RESET:decoy_stems_second_set", static_cast<int>(rr->decoy_stems_second_set));
          run_sp.setMetaValue("RESET:decoy_stems", static_cast<int>(rr->decoy_stems));
          run_sp.setMetaValue("RESET:decoy_stems_both_sets", static_cast<int>(rr->decoy_stems_both_sets));
        }
      }
      if (reuse != nullptr)
      {
        run_sp.setMetaValue("Percolator:rescore_reuse_mode", reuse->mode);
        run_sp.setMetaValue("Percolator:rescore_reuse_model", reuse->note->model);
        run_sp.setMetaValue("Percolator:rescore_reuse_max_train", static_cast<int>(reuse->max_train));
        run_sp.setMetaValue("Percolator:rescore_reuse_seed", static_cast<int>(reuse->seed));
        if (reuse->note->shift >= 0.0) { run_sp.setMetaValue("Percolator:rescore_reuse_shift", reuse->note->shift); }
      }
      run.setSearchParameters(run_sp);
    }
  }

  void ProSEAlgorithm::collectInferenceOnlyPsms_(const PeptideIdentificationList& peptide_ids, double psm_fdr,
                                                 PeptideIdentificationList& inference_only_psms) const
  {
    const double cutoff = fdr_protein_inference_psm_;
    if (cutoff <= psm_fdr) { return; }
    for (const PeptideIdentification& id : peptide_ids)
    {
      std::vector<PeptideHit> hits;
      for (const PeptideHit& hit : id.getHits())
      {
        if (hit.getScore() > psm_fdr && hit.getScore() <= cutoff) { hits.push_back(hit); }
      }
      if (hits.empty()) { continue; }
      PeptideIdentification& copy = inference_only_psms.emplace_back(id);
      copy.setHits(std::move(hits));
    }
  }

  bool ProSEAlgorithm::applyDeferredPsmFDR(const std::vector<ProteinIdentification>& protein_ids,
                                           PeptideIdentificationList& peptide_ids,
                                           bool rescored,
                                           double psm_fdr,
                                           PeptideIdentificationList* inference_only_psms) const
  {
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    if (psm_fdr <= 0.0 || protein_ids.empty() || peptide_ids.empty())
    {
      return false;
    }
    // Rescored PSMs carry Percolator q-values. The others need q-values from their search engine
    // scores, and so decoys to estimate them with.
    if (!rescored)
    {
      const auto& hits = protein_ids[0].getHits();
      const bool has_decoys = std::any_of(hits.begin(), hits.end(), [](const ProteinHit& ph)
      {
        return ph.metaValueExists("target_decoy") && ph.getMetaValue("target_decoy").toString() == "decoy";
      });
      if (!has_decoys)
      {
        return false;
      }
      annotatePsmQValues(peptide_ids); // FDR:PSM_groups as in the search's own FDR:PSM
    }
    // The PSMs the filter removes but the protein inference uses (FDR:protein_inference_PSM), on the same q-values.
    if (inference_only_psms != nullptr) { collectInferenceOnlyPsms_(peptide_ids, psm_fdr, *inference_only_psms); }
    // Filters target and decoy PSMs alike by q-value; decoy removal is the protein FDR's finalization.
    IDFilter::filterHitsByScore(peptide_ids, psm_fdr);
    return true;
  }

  // rescoring:stratum_q_values (ProSEAlgorithm::stratumQValues() in the header)
  namespace
  {
    // rescoring:stratum_q_values: grouped target-decoy competition (as mokapot's grouped confidence): the rescored best
    // hit of every spectrum belongs to the stratum 'protein_nterm' (an evidence at a protein N-terminus, or only
    // Met-clipped evidences) or 'other'; within each stratum the hits are ordered by percolator_score (ties: decoys
    // first), q = (decoys + 1) / max(targets, 1) from the top, made monotone from the bottom. Meta values
    // 'stratum' and 'stratum_q_value'; the main score stays the pooled q-value.
    void annotateStratumQValues(PeptideIdentificationList& peptide_ids)
    {
      // In parallel (serially ~1 s for a timsTOF run, by the meta value lookups by name): per spectrum, the stratum of
      // its best hit (-1: no hit, 2: not rescored), its score and whether it is a decoy.
      const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
      std::vector<int8_t> stratum_of(n, -1);
      std::vector<double> score_of(n, 0.0);
      std::vector<uint8_t> decoy_of(n, 0);
      MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      const UInt score_index = registry.getIndex("percolator_score");
      const bool score_registered = score_index != UInt(-1); // (else no hit has the meta value)
      FirstError collect_error;
#pragma omp parallel for schedule(dynamic, 1024)
      for (SignedSize i = 0; i < n; ++i)
      {
        if (peptide_ids[i].getHits().empty()) continue;
        const PeptideHit& hit = peptide_ids[i].getHits()[0];
        if (!score_registered || !hit.metaValueExists(score_index))
        {
          stratum_of[i] = 2;
          continue;
        }
        try
        {
          bool nterm = false, all_clipped = !hit.getPeptideEvidences().empty();
          for (const PeptideEvidence& ev : hit.getPeptideEvidences())
          {
            nterm = nterm || ev.getAABefore() == PeptideEvidence::N_TERMINAL_AA;
            all_clipped = all_clipped && ev.getAABefore() == 'M' && ev.getStart() == 1;
          }
          score_of[i] = static_cast<double>(hit.getMetaValue(score_index));
          decoy_of[i] = hit.isDecoy() ? 1 : 0;
          stratum_of[i] = (nterm || all_clipped) ? 1 : 0;
        }
        catch (...)
        {
          collect_error.set(i);
        }
      }
      if (collect_error.error)
      {
        // as the serial loop: the best hits before the failing one that are not rescored have their q-value
        for (SignedSize i = 0; i < collect_error.index; ++i)
        {
          if (stratum_of[i] == 2) peptide_ids[i].getHits()[0].setMetaValue("stratum_q_value", 1.0);
        }
        collect_error.rethrow();
      }
      // The rows of each stratum in spectrum order. The meta value names are registered in the order in which the
      // serial annotation set them first (MetaInfo orders the values of a hit by registry index): 'stratum_q_value'
      // first if a best hit is not rescored (its q-value was set in a first pass), else 'stratum'.
      std::vector<SignedSize> rows[2];
      bool not_rescored = false;
      for (SignedSize i = 0; i < n; ++i)
      {
        if (stratum_of[i] == 2) not_rescored = true;
        else if (stratum_of[i] >= 0) rows[stratum_of[i]].push_back(i);
      }
      UInt stratum_index = 0, q_index = 0;
      if (not_rescored) q_index = registry.registerName("stratum_q_value");
      if (!rows[0].empty() || !rows[1].empty())
      {
        stratum_index = registry.registerName("stratum");
        q_index = registry.registerName("stratum_q_value");
      }
      FirstError write_error;
      if (not_rescored)
      {
        const DataValue one(1.0); // not rescored (q = 1 in the pooled list as well)
#pragma omp parallel for schedule(static)
        for (SignedSize i = 0; i < n; ++i)
        {
          try
          {
            if (stratum_of[i] == 2) peptide_ids[i].getHits()[0].setMetaValue(q_index, one);
          }
          catch (...)
          {
            write_error.set(i);
          }
        }
        write_error.rethrow();
      }
      static const char* names[2] = {"other", "protein_nterm"};
      for (int g = 0; g < 2; ++g)
      {
        // the rows (in spectrum order) sorted as before: the same stable sort with the same comparison (also for NaN)
        std::vector<SignedSize>& r = rows[g];
        if (std::any_of(r.begin(), r.end(), [&score_of](const SignedSize i) { return std::isnan(score_of[i]); }))
        {
          // With NaN the comparison is no strict weak order, and the permutation depends on the algorithm's path (e.g.
          // on the temporary buffer it gets, counted in elements): the former sort of rows of the former size, exactly
          // (the former Row held a PeptideHit* where this one holds the spectrum index). (Only the buffer that the heap
          // grants can differ, as it could for the former code between runs.)
          struct Row { double score; bool decoy; SignedSize spectrum; };
          static_assert(sizeof(SignedSize) == sizeof(PeptideHit*) && alignof(SignedSize) == alignof(PeptideHit*), "rows of the former size");
          std::vector<Row> former(r.size());
          for (Size i = 0; i < r.size(); ++i) former[i] = {score_of[r[i]], decoy_of[r[i]] != 0, r[i]};
          std::stable_sort(former.begin(), former.end(), [](const Row& a, const Row& b)
                           { return a.score > b.score || (a.score == b.score && a.decoy && !b.decoy); });
          for (Size i = 0; i < r.size(); ++i) r[i] = former[i].spectrum;
        }
        else
        {
          // a strict weak order: every stable sort gives the same (unique) result
          const auto before = [&score_of, &decoy_of](const SignedSize a, const SignedSize b)
                              { return score_of[a] > score_of[b] || (score_of[a] == score_of[b] && decoy_of[a] && !decoy_of[b]); };
#ifdef PROSE_GNU_PARALLEL
          if (r.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::stable_sort(r.begin(), r.end(), before); }
          else
#endif
          { std::stable_sort(r.begin(), r.end(), before); }
        }
        std::vector<double> q(r.size());
        Size targets = 0, decoys = 0;
        for (Size i = 0; i < r.size(); ++i)
        {
          (decoy_of[r[i]] ? decoys : targets) += 1;
          q[i] = std::min(1.0, static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1)));
        }
        for (Size i = r.size(); i-- > 1;) q[i - 1] = std::min(q[i - 1], q[i]);
        const DataValue name(names[g]);
#pragma omp parallel for schedule(static)
        for (SignedSize i = 0; i < static_cast<SignedSize>(r.size()); ++i)
        {
          try
          {
            PeptideHit& hit = peptide_ids[r[i]].getHits()[0];
            hit.setMetaValue(stratum_index, name);
            hit.setMetaValue(q_index, q[i]);
          }
          catch (...)
          {
            write_error.set(i);
          }
        }
        write_error.rethrow();
      }
    }
  }

  void ProSEAlgorithm::stratumQValues(PeptideIdentificationList& peptide_ids)
  {
    annotateStratumQValues(peptide_ids);
  }

  void ProSEAlgorithm::applyRescoredPsmFDR_(std::vector<ProteinIdentification>& protein_ids,
                                            PeptideIdentificationList& peptide_ids,
                                            bool has_decoys,
                                            RunStatistics& stats,
                                            bool rescored,
                                            PeptideIdentificationList* inference_only_psms) const
  {
    if (rescored && stratum_q_values_)
    {
      annotateStratumQValues(peptide_ids);
      for (ProteinIdentification& run : protein_ids) run.getSearchParameters().setMetaValue("rescoring:stratum_q_values", "true");
    }
    // A run that Percolator did not rescore reports at most report:top_hits hits per spectrum, which can be fewer than
    // the rescoring:top_hits that the search kept for Percolator.
    if (!rescored && rescore_top_hits_ > 0) keepTopHits_(peptide_ids, report_top_hits_);
    // The hits per spectrum that Percolator rescored; not recorded with rescoring:top_hits 0 (= report:top_hits, as
    // earlier ProSE versions, which recorded nothing).
    if (rescored && rescore_top_hits_ > 0)
    {
      for (ProteinIdentification& run : protein_ids)
      {
        run.getSearchParameters().setMetaValue("Percolator:top_hits", static_cast<int>(rescoreTopHits_()));
      }
    }
    // PSM-level FDR: Percolator's q-values, or target-decoy q-values of the search engine scores (for a run that
    // Percolator did not rescore, also without a decoy PSM, see annotatePsmQValues()), as applyPsmFDR_() without
    // rescoring. The filter is never lifted: a threshold that the run cannot reach reports no PSM.
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    if (fdr_psm_ > 0.0 && (rescored || has_decoys))
    {
      StopWatch sw_fdr; sw_fdr.start();
      if (!rescored && !peptide_ids.empty()) { annotatePsmQValues(peptide_ids); }
      // The PSMs the filter removes but the protein inference uses (FDR:protein_inference_PSM), on the same q-values.
      if (inference_only_psms != nullptr) { collectInferenceOnlyPsms_(peptide_ids, fdr_psm_, *inference_only_psms); }
      // Filters target and decoy PSMs alike by q-value; decoy removal is the protein FDR's finalization.
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);
      stats.fdr_applied = true;
      stats.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      stats.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << std::endl;
    }
    // The proteins of the run are those its PSMs reference (see applyPsmFDR_()). A multi-file search that may rescore
    // its files together does this at the end (rescoreFilesTogether_()) for the runs that Percolator did not rescore on
    // their own (their PSMs may change); a rescored run is final here (and can be written while the next file is searched).
    if (!pool_rescoring_ || rescored)
    {
      removeUnreportedProteins(protein_ids, peptide_ids,
                               inference_only_psms != nullptr ? *inference_only_psms : PeptideIdentificationList());
    }
  }

  void ProSEAlgorithm::finishDeferredRescoring_(MultiFileSearchResult& mfres, const std::vector<Size>& deferred) const
  {
    reuse_defer_ = false;
    if (deferred.empty()) { return; }
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    std::vector<ReuseRun_> runs;
    for (const Size d : deferred)
    {
      SearchResult& result = mfres.per_file[d];
      if (result.exit_code == ExitCodes::EXECUTION_OK) { runs.push_back({d, &result.protein_ids, &result.peptide_ids}); }
    }
    std::shared_ptr<const ReuseModel_> model;
    StopWatch sw; sw.start();
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse " << reuse_mode_ << ": training the Percolator models on " << runs.size()
                    << " file(s)..." << std::endl;
    try
    {
      model = trainReuseModel_(runs, threads);
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << model->describe() << " in " << sw.getClockTime() << " s" << std::endl;
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "[ProSE] rescore_reuse: cannot train the models (" << e.what()
                      << "); every file is rescored by Percolator trained on it, as with 'per_file'." << std::endl;
    }
    // 'static': the files after the training files are rescored with these models after their search (applyPsmFDR_())
    reuse_model_ = reuse_mode_ == "static" ? model : nullptr;
    for (const Size d : deferred)
    {
      SearchResult& result = mfres.per_file[d];
      if (result.exit_code != ExitCodes::EXECUTION_OK) { continue; }
      OPENMS_LOG_INFO << "[ProSE] Rescoring the PSMs of " << result.stats.input_file << " with Percolator..." << std::endl;
      StopWatch sw_file;
      sw_file.start();
      result.percolator = model != nullptr
                          ? rescoreWithReuseModel_(result.protein_ids, result.peptide_ids, d, *model, false, threads, nullptr)
                          : rescoreWithPercolator_(result.protein_ids, result.peptide_ids, threads, nullptr, percolator_training_,
                                                   Param());
      // (the PSMs only the protein inference uses: with FDR:protein, as applyPsmFDR_() collects them)
      applyRescoredPsmFDR_(result.protein_ids, result.peptide_ids, mfres.have_decoys, result.stats, result.percolator.applied,
                           fdr_protein_ > 0.0 ? &result.inference_only_peptide_ids : nullptr);
      // The statistics of the final PSMs (the search collected them before its rescoring); the spectra are gone.
      const Size ms2_spectra = result.stats.ms2_spectra;
      result.stats.charge_histogram.clear();
      result.stats.missed_cleavage_histogram.clear();
      collectRunStatistics_(PeakMap(), result.protein_ids, result.peptide_ids, result.stats);
      result.stats.ms2_spectra = ms2_spectra;
      sw_file.stop();
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << result.stats.input_file << " rescored in " << sw_file.getClockTime() << " s"
                      << std::endl;
    }
    sw.stop();
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << deferred.size() << " file(s) rescored in " << sw.getClockTime() << " s" << std::endl;
  }

  // static
  void ProSEAlgorithm::keepTopHits_(PeptideIdentificationList& peptide_ids, Size n)
  {
    // The hits of a spectrum are in the order of postProcessHits_() (best search engine score first): the first n
    // are its top n. Nothing changes if no spectrum has more (e.g. always with rescoring:top_hits 0).
    if (std::none_of(peptide_ids.begin(), peptide_ids.end(),
                     [n](const PeptideIdentification& pid) { return pid.getHits().size() > n; }))
    {
      return;
    }
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      if (hits.size() > n) hits.erase(hits.begin() + n, hits.end()); // (destroys hits only: cannot throw)
    }
  }

  namespace
  {
    // rescoring:terminal_features: 'protein_nterm' = 1 if an evidence of the hit is at a protein N-terminus; 'met_clipped'
    // = 1 if every evidence is the peptide after an initial Met (aa_before M, start 1: in the search space only through
    // peptide:clip_nterm_methionine). Every hit (targets and decoys); the names join the rescoring features.
    void annotateTerminalFeatures(std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids)
    {
      static const std::string nterm_name = "protein_nterm", clip_name = "met_clipped";
      // In parallel (serially ~0.8 s for a timsTOF run, by the meta value lookups by name), with the names registered
      // as by the first hit of the serial loop (MetaInfo orders the values of a hit by registry index)
      if (std::any_of(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& pid) { return !pid.getHits().empty(); }))
      {
        MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
        const UInt nterm_index = registry.registerName(nterm_name);
        const UInt clip_index = registry.registerName(clip_name);
        const DataValue one(1.0), zero(0.0);
        FirstError error; // (an allocation failure)
#pragma omp parallel for schedule(dynamic, 1024)
        for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
        {
          try
          {
            for (PeptideHit& hit : peptide_ids[i].getHits())
            {
              bool nterm = false, all_clipped = !hit.getPeptideEvidences().empty();
              for (const PeptideEvidence& ev : hit.getPeptideEvidences())
              {
                nterm = nterm || ev.getAABefore() == PeptideEvidence::N_TERMINAL_AA;
                all_clipped = all_clipped && ev.getAABefore() == 'M' && ev.getStart() == 1;
              }
              hit.setMetaValue(nterm_index, nterm ? one : zero);
              hit.setMetaValue(clip_index, (!nterm && all_clipped) ? one : zero);
            }
          }
          catch (...)
          {
            error.set(i);
          }
        }
        error.rethrow();
      }
      for (ProteinIdentification& prot : protein_ids)
      {
        ProteinIdentification::SearchParameters sp = prot.getSearchParameters();
        StringList features;
        if (sp.metaValueExists("extra_features"))
        {
          const std::string existing = sp.getMetaValue("extra_features").toString();
          if (!existing.empty()) features = ListUtils::create<std::string>(existing);
        }
        for (const std::string& name : {nterm_name, clip_name})
        {
          if (std::find(features.begin(), features.end(), name) == features.end()) features.push_back(name);
        }
        sp.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
        sp.setMetaValue("rescoring:terminal_features", "true");
        prot.setSearchParameters(std::move(sp));
      }
    }
  }

  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::applyPsmFDR_(std::vector<ProteinIdentification>& protein_ids,
                                                                   PeptideIdentificationList& peptide_ids,
                                                                   bool has_decoys,
                                                                   RunStatistics& stats,
                                                                   PeptideIdentificationList* pre_fdr_psms,
                                                                   PeptideIdentificationList* inference_only_psms) const
  {
    PercolatorRescoring rescoring;
    // The PSMs between FDR:PSM and FDR:protein_inference_PSM are only needed for the protein FDR.
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    PeptideIdentificationList* inference_only = fdr_protein_ > 0.0 ? inference_only_psms : nullptr;
    if (terminal_features_) annotateTerminalFeatures(protein_ids, peptide_ids); // also for a .pin without rescoring
    // The copy for a Percolator .pin (setKeepPreFdrPsms()) holds the search engine scores, every reported hit and the
    // decoys: it is the input of a later rescoring, so it is taken before the rescoring as well as before FDR:PSM.
    if (pre_fdr_psms != nullptr)
    {
      *pre_fdr_psms = peptide_ids;
      // (with rescoring:top_hits above report:top_hits the search kept more hits, for Percolator). A multi-file search
      // that may rescore its files together keeps the hits Percolator rescores as well; rescoreFilesTogether_() uses
      // the first rescoring:top_hits of them and cuts the .pin copy to report:top_hits.
      if (rescore_top_hits_ > 0)
      {
        keepTopHits_(*pre_fdr_psms, pool_rescoring_ ? std::max(report_top_hits_, rescoreTopHits_()) : report_top_hits_);
      }
    }
    if (rescore_)
    {
      // Percolator rescores the first rescoring:top_hits hits of each spectrum (the search kept more for the .pin copy).
      // Not with decoys 'ignore': the search kept report:top_hits, and Percolator skips the run (no decoys). Nothing is
      // cut with rescoring:top_hits 0 (identifications that the caller's list held before the search included).
      if (rescore_top_hits_ > 0 && decoy_mode_ != DecoyMode_::IGNORE) keepTopHits_(peptide_ids, rescoreTopHits_());
      // Percolator gets all target and decoy PSMs; the PSM FDR follows on its q-values (or, if it
      // skipped the run, on the search engine scores).
      startDeferredRead(); // the next file of a multi-file search, with threads Percolator leaves idle
      if (reuse_defer_)
      {
        // rescore_reuse joint/static: rescored (and FDR:PSM applied) by finishDeferredRescoring_() once the models are trained
        rescoring.skipped_reason = "The rescoring waits for the models of the multi-file search (rescore_reuse).";
        return rescoring;
      }
#ifdef _OPENMP
      const int threads = omp_get_max_threads();
#else
      const int threads = 1;
#endif
      OPENMS_LOG_INFO << "[ProSE] Rescoring the PSMs with Percolator..." << std::endl;
      // The hits that the rescoring drops (all but the best per spectrum) are destroyed on a helper thread while the
      // search and its output continue: destroying the 65,000 dropped hits of an Astral run took the rescoring's 16
      // threads 42-60 ms (allocator contention; 150 ms on one thread) while the search waited. The phases after the
      // rescoring keep 11-14 of 16 threads busy. The helper is joined when the holder goes: at the next rescoring of
      // this object, or with this object (its last copy). Not on Windows, where a holder destroyed while a DLL is
      // unloaded would join under the loader lock, which the ending helper thread may need.
      released_hits_.reset(); // the release of the search before (multi-file search)
#ifdef OPENMS_WINDOWSPLATFORM
      constexpr bool release_later = false;
#else
      const bool release_later = threads > 1;
#endif
      std::vector<std::vector<PeptideHit>> dropped;
      // rescore_method 'reset' with a directory as reset:stem_table: one table per run, named after its spectra file
      Param extra = percolator_extra_;
      if (extra.exists("reset:stem_table"))
      {
        std::string path = extra.getValue("reset:stem_table").toString();
        if (!path.empty() && (path.back() == '/' || File::isDirectory(path)))
        {
          std::string base = reset_run_name_.empty() ? std::string("run") : File::basename(reset_run_name_);
          if (const auto dot = base.rfind('.'); dot != std::string::npos && dot > 0) base.erase(dot);
          extra.setValue("reset:stem_table", (path.back() == '/' ? path : path + "/") + base + ".reset_stems.tsv");
        }
      }
      if (reuse_model_ != nullptr && reuse_mode_ == "static")
      {
        // rescore_reuse static: a file after the training files
        rescoring = rescoreWithReuseModel_(protein_ids, peptide_ids, reuse_file_index_, *reuse_model_, true, threads,
                                           release_later ? &dropped : nullptr);
      }
      else if (reuse_mode_ == "own")
      {
        std::shared_ptr<const ReuseModel_> own;
        try
        {
          own = trainReuseModel_({ReuseRun_{reuse_file_index_, &protein_ids, &peptide_ids}}, threads);
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_WARN << "[ProSE] rescore_reuse: cannot train the models (" << e.what()
                          << "); Percolator is trained on the run as with 'per_file'." << std::endl;
        }
        rescoring = own != nullptr
                    ? rescoreWithReuseModel_(protein_ids, peptide_ids, reuse_file_index_, *own, false, threads,
                                             release_later ? &dropped : nullptr)
                    : rescoreWithPercolator_(protein_ids, peptide_ids, threads, release_later ? &dropped : nullptr,
                                             percolator_training_, Param());
      }
      else
      {
        rescoring = rescoreWithPercolator_(protein_ids, peptide_ids, threads, release_later ? &dropped : nullptr,
                                           percolator_training_, extra);
      }
      if (!dropped.empty())
      {
        try
        {
          auto release = std::make_shared<HitRelease>();
          // With a quarter of the threads (on one thread, the 2 million hits that the rescoring of a timsTOF run drops took
          // 2.5 s, longer than the rest of the search: the holder then waited for it)
          release->thread = std::thread([hits = std::move(dropped), release_threads = std::max(1, threads / 4)]() mutable
          {
#ifdef _OPENMP
            omp_set_num_threads(release_threads);
#endif
            (void)release_threads;
            freeInParallel(hits);
          });
          released_hits_ = std::move(release);
        }
        catch (...) {} // if the thread cannot be started, the hits are destroyed in this function
      }
      takeBackThreadsOfEndedRead();
      applyRescoredPsmFDR_(protein_ids, peptide_ids, has_decoys, stats, rescoring.applied, inference_only);
      return rescoring;
    }

    if (fdr_psm_ > 0.0 && (rescoring.applied || has_decoys))
    {
      // PSM-level FDR: Percolator's q-values, or target-decoy q-values of the search engine scores (also for a run
      // that Percolator did not rescore, and for one without a decoy PSM, see annotatePsmQValues()); target and
      // decoy PSMs are filtered alike by the threshold. No decoy-specific stripping happens here (decoupled from
      // decoy removal); the decoys that pass are kept because downstream protein-level FDR and cross-file merging
      // need them. Categorical decoy removal happens only at protein-FDR finalization (applyCompleteSetProteinFDR()).
      // The filter is never lifted: a threshold that the run cannot reach reports no PSM.
      StopWatch sw_fdr; sw_fdr.start();
      if (!rescoring.applied && !peptide_ids.empty()) { annotatePsmQValues(peptide_ids); }
      if (inference_only != nullptr) { collectInferenceOnlyPsms_(peptide_ids, fdr_psm_, *inference_only); }
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);
      stats.fdr_applied = true;
      stats.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      stats.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << std::endl;
    }

    // The proteins of the run are those its PSMs reference: PeptideIndexing added the proteins of all hits, and the
    // filters above removed PSMs and hits. Those of the PSMs for the protein inference stay until the protein FDR.
    // (A multi-file search that may rescore its files together does this at the end, rescoreFilesTogether_(): the
    // PSMs may change.)
    if (!pool_rescoring_)
    {
      removeUnreportedProteins(protein_ids, peptide_ids, inference_only != nullptr ? *inference_only : PeptideIdentificationList());
    }
    return rescoring;
  }

  namespace
  {
    // Starts the OpenMP threads of the calling thread, if they do not run yet, and returns their
    // number: 1 without OpenMP, with a single OpenMP thread, and inside a parallel region.
    Size startOpenMPThreads()
    {
      Size threads = 0;
#pragma omp parallel reduction(+ : threads)
      {
        ++threads;
      }
      return threads;
    }

    // The index build waits for the spectra that are read meanwhile, to index only the peptides they can reach, if
    // generating the fragments it skips takes longer than the rest of the read: with few threads (measured on
    // 8,000-spectrum files: faster for all instrument types up to 4 threads, slower for some from 6 threads on) and a
    // read that is short against the build (the read is serial, the build parallel): at most
    // MAX_SPECTRA_BYTES_PER_PEPTIDE / threads bytes of spectra file per peptide (measured with 9 M peptides: faster
    // with 11, 27 and 41 bytes per peptide at 2 threads and with 11 at 4 threads, even with 27, slower with 41).
    // With more threads, the single-file search always waits (and helps with the rest of the read, see
    // readerThreads()): the read then mostly ends while the peptides are generated and sorted, and the fragment passes over the third or so of the peptides that no spectrum reaches take longer
    // than the wait (measured at 16 threads on Orbitrap, Astral, timsTOF and TMT files of 8,000 spectra: 0.1-0.2 s
    // shorter per file; waiting for a read with a quarter of the threads was up to 0.1 s slower for some).
    constexpr Size MAX_THREADS_WAITING_FOR_SPECTRA = 4;
    constexpr UInt64 MAX_SPECTRA_BYTES_PER_PEPTIDE = 100;

    // FASTAFile::load() with the OpenMP threads: the file is cut into pieces at starts of entries,
    // the threads read the pieces with FASTAFile::readNext(), and the pieces are joined in file
    // order. Same entries as load(): readNext() reads an entry from its '>' up to a line break
    // followed by '>'. A piece starts at a '>' after a line break whose line holds no '>' (not a
    // header line, so usually the last line of a sequence), where the reader of load() starts an
    // entry as well, and the reader of the piece before stops when it arrives there. Whenever this
    // does not work out (small file, no such place, a reader that fails or does not arrive at the
    // next piece, e.g. as it skipped an entry without sequence there), load() reads the file and
    // reports errors as before.
    void readFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries)
    {
      std::vector<std::streamoff> pieces{0}; // where each piece starts; the first one where readStart() starts
#ifdef _OPENMP
      if (!omp_in_parallel() && omp_get_max_threads() > 1)
      {
        std::ifstream in(filename, std::ios::binary);
        in.seekg(0, std::ios::end);
        const std::streamoff size = in.tellg();
        // several pieces per thread: the threads do not read equally fast
        const std::streamoff count = std::min<std::streamoff>(4 * omp_get_max_threads(), size / (std::streamoff(1) << 18));
        std::string window(Size(1) << 16, '\0');
        for (std::streamoff i = 1; i < count && in.good(); ++i)
        {
          const std::streamoff from = size / count * i;
          in.seekg(from);
          in.read(window.data(), static_cast<std::streamsize>(window.size()));
          const std::string_view text(window.data(), static_cast<Size>(in.gcount()));
          in.clear(); // reading up to the end of the file is fine
          // lines of the window that are complete: between two line breaks
          for (Size line = text.find('\n'); line != std::string_view::npos && line + 1 < text.size();)
          {
            const Size next = text.find('\n', line + 1);
            if (next == std::string_view::npos || next + 1 >= text.size()) { break; }
            if (text[next + 1] == '>' && text.substr(line + 1, next - line - 1).find('>') == std::string_view::npos)
            {
              if (from + static_cast<std::streamoff>(next) + 1 > pieces.back()) { pieces.push_back(from + static_cast<std::streamoff>(next) + 1); }
              break;
            }
            line = next;
          }
        }
      }
#endif
      if (pieces.size() < 2)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      std::vector<std::vector<FASTAFile::FASTAEntry>> read(pieces.size());
      bool complete = true;
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize i = 0; i < static_cast<SignedSize>(pieces.size()); ++i)
      {
        bool ok = false;
        try // exceptions must not leave the parallel region
        {
          FASTAFile file;
          FASTAFile::FASTAEntry entry;
          file.readStart(filename);
          ok = i == 0 || file.setPosition(pieces[i]);
          if (ok && i + 1 == static_cast<SignedSize>(pieces.size()))
          { // the last piece: up to the end of the file, as load()
            while (file.readNext(entry)) { read[i].push_back(std::move(entry)); }
          }
          else if (ok)
          {
            const std::streampos end = pieces[i + 1];
            while (ok && file.position() < end)
            {
              ok = file.readNext(entry);
              if (ok) { read[i].push_back(std::move(entry)); }
            }
            ok = ok && file.position() == end;
          }
        }
        catch (...)
        {
          ok = false;
        }
        if (!ok)
        {
#pragma omp critical (ProSEAlgorithm_loadFASTA)
          complete = false;
        }
      }
      if (!complete)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      Size total = 0;
      for (const auto& piece : read) { total += piece.size(); }
      entries.clear();
      entries.reserve(total);
      for (auto& piece : read)
      {
        entries.insert(entries.end(), std::make_move_iterator(piece.begin()), std::make_move_iterator(piece.end()));
      }
    }

    // Lets the thread that waits for a spectra reader (BackgroundSpectra::get()) parse chunks of its loadMzMLChunked()
    // with its own OpenMP threads, which would idle otherwise, instead of only waiting for the end of the read.
    class ChunkHelp
    {
    public:
      /// By the reader: from now on, help() runs @p work (which parses chunks until none is left) ...
      void open(std::function<void()> work)
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!done_) { work_ = std::move(work); }
        changed_.notify_all();
      }

      /// ... until the reader closes: waits for the helpers that run it. Also by the reader once it reads no more chunks.
      void close() noexcept
      {
        std::unique_lock<std::mutex> lock(mutex_);
        work_ = nullptr;
        done_ = true;
        changed_.notify_all();
        changed_.wait(lock, [this]() { return helpers_ == 0; });
      }

      /// By the waiting thread: runs the work (with the OpenMP threads of this thread) if the reader opens it, or returns
      /// once the reader closes. The thread then waits for the spectra.
      void help()
      {
        std::function<void()> work;
        {
          std::unique_lock<std::mutex> lock(mutex_);
#ifdef _OPENMP
          waiting_threads_ = omp_get_max_threads();
#else
          waiting_threads_ = 1;
#endif
          changed_.wait(lock, [this]() { return work_ || done_; });
          if (!work_) { return; }
          work = work_;
          ++helpers_;
        }
        const OnScopeExit leave([this]()
        {
          std::lock_guard<std::mutex> lock(mutex_);
          --helpers_;
          changed_.notify_all();
        });
        work();
      }

      /// By the reader: whether the reader took its threads from the thread that waits (see readBorrowsThreads())
      void setLent(const bool lent) { lent_ = lent; }

      /// By the reader, after the chunks: the OpenMP threads for its remaining parallel work, its own @p reader_threads
      /// and, once the thread that waits for the spectra does (help() was called; it does nothing else until it has
      /// them), those of that thread, too (if they were lent: the threads of both; otherwise as many as the larger team)
      int joinThreads(const int reader_threads)
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (waiting_threads_ == 0) { return reader_threads; }
        return lent_ ? reader_threads + waiting_threads_ : std::max(reader_threads, waiting_threads_);
      }

    private:
      std::mutex mutex_;
      std::condition_variable changed_;
      std::function<void()> work_;
      int helpers_ = 0;
      bool done_ = false;
      int waiting_threads_ = 0; // OpenMP threads of the thread that waits for the spectra; 0: it does not wait yet
      bool lent_ = false;
    };

    // The database of a search: readFASTA(), without repeated entries (ProSEAlgorithm::removeDuplicateProteins(), which
    // names the file in its warning or error), before the database is used.
    void loadFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries)
    {
      readFASTA(filename, entries);
      ProSEAlgorithm::removeDuplicateProteins(entries, filename);
    }

    // Removes the spectra of @p spectra that @p keep does not accept (none without @p keep)
    void removeUnkept(std::vector<MSSpectrum>& spectra, const std::function<bool(const MSSpectrum&)>& keep)
    {
      if (!keep) { return; }
      spectra.erase(std::remove_if(spectra.begin(), spectra.end(), [&keep](const MSSpectrum& s) { return !keep(s); }), spectra.end());
    }

    // MzMLFile::load() with @p threads OpenMP threads (at least 2). The spectra are cut into chunks at
    // their <spectrum> start tags; every chunk is parsed together with the file's header and closing
    // tags by an MzMLFile of its own (MzMLFile::loadBuffer), and the spectra are joined in file order.
    // The last chunk runs up to </run>, so it also holds the chromatograms. Each thread reads its chunks
    // from the file itself, so the data held beyond the spectra is at most about 9 MB per thread,
    // whatever the size of the file. The header is parsed first on this thread: it gives the experimental
    // settings, registers its meta value names in file order and initialises the XML parser.
    // Returns false if the file does not have the layout the chunks rely on (compressed; a comment,
    // CDATA section or DOCTYPE before the end of the run; not exactly one spectrum list; a count
    // attribute other than the number of spectra found; few spectra) or a chunk fails to parse. The
    // caller then reads the file with MzMLFile::load(), which also reports its errors.
    // With @p keep only the spectra it accepts are kept (by each chunk, right after parsing it).
    bool loadMzMLChunked(const std::string& filename, const PeakFileOptions& options, int threads, PeakMap& exp,
                         ChunkHelp* help = nullptr, int chunk_threads = 0,
                         const std::function<bool(const MSSpectrum&)>& keep = {})
    {
      std::ifstream in(filename, std::ios::binary);
      if (!in) { return false; }
      in.seekg(0, std::ios::end);
      const Size size = static_cast<Size>(in.tellg());
      char magic[2] = {0, 0};
      in.seekg(0);
      in.read(magic, 2);
      if (!in || (magic[0] == 'B' && magic[1] == 'Z') || (magic[0] == '\x1f' && magic[1] == '\x8b') || (magic[0] == 'P' && magic[1] == 'K'))
      {
        return false; // compressed (bzip2, gzip, zip) or too short
      }

      // 1. Where the markup of interest starts, found by the threads in blocks of the file.
      enum class Mark : unsigned char { spectrum, list_start, list_end, run_end, declaration };
      const auto is_tag = [](std::string_view text, std::string_view name)
      {
        return text.size() > name.size() && text.substr(0, name.size()) == name
               && std::string_view(" \t\r\n>/").find(text[name.size()]) != std::string_view::npos;
      };
      constexpr Size block = Size(1) << 20;
      constexpr Size lookahead = 16; // longer than the longest name compared below, "</spectrumList" plus one
      const Size n_blocks = (size + block - 1) / block;
      std::vector<std::vector<std::pair<Size, Mark>>> marks(n_blocks);
      bool ok = true;
#pragma omp parallel num_threads(threads)
      {
        std::ifstream file; // opened by the thread's first block, in the try below
        std::string buffer;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(n_blocks); ++sb)
        {
          bool read = false;
          try // exceptions (e.g. std::bad_alloc of the buffers) must not leave the parallel region
          {
            if (!file.is_open()) { file.open(filename, std::ios::binary); }
            const Size b = static_cast<Size>(sb);
            const Size from = b * block;
            const Size end = std::min(block, size - from); // marks starting before end belong to this block
            buffer.resize(std::min(block + lookahead, size - from));
            file.seekg(static_cast<std::streamoff>(from));
            file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
            read = static_cast<bool>(file);
            for (Size i = 0; read && i < end; ++i)
            {
              const char* p = static_cast<const char*>(std::memchr(buffer.data() + i, '<', end - i));
              if (p == nullptr) { break; }
              i = static_cast<Size>(p - buffer.data());
              const std::string_view text(p, std::min(lookahead, buffer.size() - i));
              if (is_tag(text, "<spectrumList")) { marks[b].emplace_back(from + i, Mark::list_start); }
              else if (is_tag(text, "<spectrum")) { marks[b].emplace_back(from + i, Mark::spectrum); }
              else if (is_tag(text, "</spectrumList")) { marks[b].emplace_back(from + i, Mark::list_end); }
              else if (is_tag(text, "</run")) { marks[b].emplace_back(from + i, Mark::run_end); }
              else if (text.substr(0, 2) == "<!") { marks[b].emplace_back(from + i, Mark::declaration); }
            }
          }
          catch (...)
          {
            read = false;
          }
          if (!read)
          {
#pragma omp critical (ProSEAlgorithm_loadMzMLChunked)
            ok = false;
            file.clear();
          }
        }
      }
      constexpr Size none = std::numeric_limits<Size>::max();
      Size list_start = none, list_end = none, run_end = none;
      std::vector<Size> starts; // of the spectra
      for (const auto& block_marks : marks)
      {
        for (const auto& [pos, mark] : block_marks)
        {
          switch (mark)
          {
            case Mark::list_start: ok = ok && list_start == none; list_start = pos; break;
            case Mark::spectrum: ok = ok && list_start != none && list_end == none; starts.push_back(pos); break;
            case Mark::list_end: ok = ok && list_start != none && list_end == none; list_end = pos; break;
            case Mark::run_end: if (run_end == none) { ok = ok && list_end != none; run_end = pos; } break;
            case Mark::declaration: ok = ok && run_end != none; break;
          }
        }
      }
      marks = {};
      if (!ok || run_end == none || starts.size() < 64) { return false; }

      // 2. The header (everything before the first spectrum) and the count attribute of its
      //    <spectrumList>, which each chunk sets to its own number of spectra.
      std::string header(starts[0], '\0');
      in.seekg(0);
      in.read(header.data(), static_cast<std::streamsize>(header.size()));
      if (!in) { return false; }
      Size count_begin = none, count_end = none;
      for (Size p = list_start + std::string_view("<spectrumList").size();;)
      {
        p = header.find_first_not_of(" \t\r\n", p);
        if (p == std::string::npos) { return false; }
        if (header[p] == '>') { break; }
        const Size name_end = header.find_first_of("= \t\r\n>", p);
        Size q = header.find_first_not_of(" \t\r\n", name_end);
        if (name_end == std::string::npos || q == std::string::npos || header[q] != '=') { return false; }
        q = header.find_first_not_of(" \t\r\n", q + 1);
        if (q == std::string::npos || (header[q] != '"' && header[q] != '\'')) { return false; }
        const Size value_end = header.find(header[q], q + 1);
        if (value_end == std::string::npos) { return false; }
        if (std::string_view(header).substr(p, name_end - p) == "count") { count_begin = q + 1; count_end = value_end; }
        p = value_end + 1;
      }
      if (count_begin == none || header.compare(count_begin, count_end - count_begin, std::to_string(starts.size())) != 0)
      {
        return false;
      }
      const std::string end_of_run = std::string("</run></mzML>") + (header.find("<indexedmzML") != std::string::npos ? "</indexedmzML>" : "");
      const auto chunk_text = [&](Size n_spectra, std::string& text)
      {
        text.assign(header, 0, count_begin);
        text += std::to_string(n_spectra);
        text.append(header, count_end, std::string::npos);
      };

      // 3. Chunks of about equal size: two per thread (the threads do not parse equally fast; every chunk
      //    costs the set-up of a parser and handler), of at most 8 MB, and of at least 16 spectra. A read in the
      //    background counts chunk_threads instead (see chunkThreads()).
      constexpr Size max_chunk_bytes = Size(8) << 20;
      const Size n = starts.size();
      const Size bytes = list_end - starts[0];
      const Size n_chunks = std::min(n / 16, std::max(2 * static_cast<Size>(std::max(threads, chunk_threads)), (bytes + max_chunk_bytes - 1) / max_chunk_bytes));
      std::vector<Size> first{0}; // the first spectrum of each chunk
      for (Size c = 1; c < n_chunks; ++c)
      {
        const Size f = static_cast<Size>(std::lower_bound(starts.begin(), starts.end(), starts[0] + bytes / n_chunks * c) - starts.begin());
        if (f > first.back() && f < n) { first.push_back(f); }
      }
      first.push_back(n);

      // A chunk decodes its spectra on its thread; in batches of the default size, so that it does not
      // hold the encoded data of all its spectra at once.
      PeakFileOptions chunk_options = options;
      chunk_options.setMaxDataPoolSize(PeakFileOptions().getMaxDataPoolSize());
      try
      {
        std::string text;
        chunk_text(0, text);
        text += "</spectrumList>" + end_of_run;
        MzMLFile mzml;
        mzml.getOptions() = chunk_options;
        mzml.loadBuffer(text, exp);
      }
      catch (...)
      {
        return false;
      }
      // The spectra that each chunk keeps, until all are parsed and put in their places (below)
      const Size n_chunks_used = first.size() - 1;
      std::vector<std::vector<MSSpectrum>> parts(n_chunks_used);
      std::vector<Size> kept(n_chunks_used, 0);
      std::vector<MSChromatogram> chromatograms;
      // The chunks are handed out one by one, to the threads of this reader and to those of a helper (another
      // contention group: no OpenMP critical between the two).
      std::atomic<Size> next_chunk{0};
      std::atomic<bool> chunks_ok{true};
      const auto parse_chunks = [&]([[maybe_unused]] const int team)
      {
#pragma omp parallel num_threads(team)
      {
        std::ifstream file; // opened by the thread's first chunk, in the try below
        std::string chunk;
        for (Size c = next_chunk++; c < n_chunks_used; c = next_chunk++)
        {
          bool parsed = false;
          try // exceptions (a parse error, std::bad_alloc of the chunk text, ...) must not leave the parallel region
          {
            if (!file.is_open()) { file.open(filename, std::ios::binary); }
            const bool last = c + 1 == n_chunks_used;
            const Size from = starts[first[c]];
            const Size to = last ? run_end : starts[first[c + 1]];
            chunk_text(first[c + 1] - first[c], chunk);
            const Size at = chunk.size();
            chunk.resize(at + (to - from));
            file.seekg(static_cast<std::streamoff>(from));
            file.read(chunk.data() + at, static_cast<std::streamsize>(to - from));
            parsed = static_cast<bool>(file);
            if (parsed)
            {
              chunk += last ? end_of_run : "</spectrumList>" + end_of_run;
              MzMLFile mzml;
              mzml.getOptions() = chunk_options;
              PeakMap part;
              mzml.loadBuffer(chunk, part);
              parsed = part.size() <= first[c + 1] - first[c];
              if (parsed) { removeUnkept(part.getSpectra(), keep); }
              if (parsed)
              {
                for (MSSpectrum& spectrum : part) { spectrum.sortByPosition(); } // see loadMS2Spectra()
                kept[c] = part.size();
                if (last) { chromatograms = std::move(part.getChromatograms()); }
                parts[c] = std::move(part.getSpectra());
              }
            }
          }
          catch (...)
          {
            parsed = false;
          }
          if (!parsed)
          {
            chunks_ok = false;
            file.clear();
          }
        }
      }
      };
      {
        // closed (waiting for the chunks of the helper) before what the chunks write to goes, also on an exception
        const OnScopeExit end_help([help]() { if (help != nullptr) { help->close(); } });
#ifdef _OPENMP
        // At least 32 parsers together (each holds a chunk of up to 8 MB): the 32 chunks of a file of up to 256 MB in
        // one round. A reader alone has at most 16 threads (no faster with more), but a thread that waits for the
        // spectra has nothing else to do (at 64 threads, where the reader has 16 threads of its own and the index build
        // waits for it: the chunks of astral_B1 parsed in 91 ms instead of 136 ms, those of hfx_B1 in 62 ms as before).
        // Up to all the threads of the search (joinThreads(): the reader's and the waiting thread's, as many as it was
        // given), as a large file has many more chunks than parsers: the 600 chunks of a 4.9 GB timsTOF HeLa run were
        // parsed in 2.2 s by 32 parsers and in 1.4 s by 64 (once StringManager::parseInt() no longer serialised them).
        constexpr int min_parsers = 32;
        if (help != nullptr)
        {
          help->open([&]()
          {
            const int team = std::min(omp_get_max_threads(), std::max(min_parsers, help->joinThreads(threads)) - threads);
            if (team >= 1) { parse_chunks(team); }
          });
        }
#endif
        parse_chunks(threads);
      }
      if (!ok || !chunks_ok) { return false; }

      // The rest in parallel, also with the threads of a thread that waits for the spectra: the spectra are put into
      // their places in the list, in retention time order as MSExperiment::sortSpectra(false) puts them (equal ones in
      // file order), unless a retention time is NaN (then in file order: the order of sortSpectra() depends on its sort,
      // which loadMS2Spectra() then runs), and MSExperiment::updateRanges(). (A list in file order, sorted afterwards,
      // took 0.29 s on one thread for the 619,000 spectra of a timsTOF HeLa run, which are not in retention time order:
      // 0.2 s to move them along the cycles of the permutation.)
#ifdef _OPENMP
      const int own_threads = omp_get_max_threads();
      omp_set_num_threads(help != nullptr ? help->joinThreads(threads) : threads);
      const OnScopeExit restore([own_threads]() { omp_set_num_threads(own_threads); });
#endif
      std::vector<Size> offset(n_chunks_used + 1, 0); // the first spectrum of each chunk among the spectra kept
      for (Size c = 0; c < n_chunks_used; ++c) { offset[c + 1] = offset[c] + kept[c]; }
      const Size kept_total = offset.back();
      std::vector<double> rts(kept_total);
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize c = 0; c < static_cast<SignedSize>(n_chunks_used); ++c)
      {
        for (Size j = 0; j < kept[c]; ++j) { rts[offset[c] + j] = parts[c][j].getRT(); }
      }
      const bool any_nan = std::any_of(rts.begin(), rts.end(), [](double rt) { return std::isnan(rt); });
      const bool in_order = std::is_sorted(rts.begin(), rts.end());
      std::vector<Size> place; // the place of the spectrum with file position i; empty: i
      if (!any_nan && !in_order)
      {
        const std::vector<Size> order = orderByKeyInParallel(rts);
        place.resize(kept_total);
#pragma omp parallel for schedule(static)
        for (SignedSize k = 0; k < static_cast<SignedSize>(kept_total); ++k) { place[order[k]] = static_cast<Size>(k); }
      }
      std::vector<MSSpectrum> spectra;
      resizeWithParallelFaults(spectra, kept_total);
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize c = 0; c < static_cast<SignedSize>(n_chunks_used); ++c)
      {
        for (Size j = 0; j < kept[c]; ++j)
        {
          const Size i = offset[c] + j;
          spectra[place.empty() ? i : place[i]] = std::move(parts[c][j]);
        }
        std::vector<MSSpectrum>().swap(parts[c]);
      }
      exp.setSpectra(std::move(spectra));
      exp.setChromatograms(std::move(chromatograms));
      exp.setLoadedFileType(filename);
      exp.setLoadedFilePath(filename);
      exp.updateRanges();
      return true;
    }

    // The MS2 spectra (with @p survey_scans also the MS1 spectra) of a spectrum file (mzML, Bruker .d or Thermo .raw), sorted by RT, read with
    // @p threads OpenMP threads: an mzML file in chunks parsed in parallel (loadMzMLChunked) with two
    // or more (or @p chunk_threads of two or more, with a helper: @p help, see ChunkHelp), otherwise with FileHandler.
    // With @p keep only the spectra it accepts (see ProSEAlgorithm::loadSpectra()).
    PeakMap loadMS2Spectra(const std::string& filename, int threads, ChunkHelp* help = nullptr, int chunk_threads = 0,
                           bool survey_scans = false, const std::function<bool(const MSSpectrum&)>& keep = {})
    {
      PeakMap spectra;
      FileHandler f;
      f.getOptions().clearMSLevels();
      f.getOptions().addMSLevel(2);
      if (survey_scans) { f.getOptions().addMSLevel(1); } // chimeric:enabled (see ProSEAlgorithm::search(spectra, ctx, ...))
      // MzMLHandler decodes the spectra in a parallel region every maxDataPoolSize spectra (100 by
      // default): with several threads, fewer and larger regions
      if (threads > 1) { f.getOptions().setMaxDataPoolSize(1000); }
      const bool chunked = (threads > 1 || (help != nullptr && chunk_threads > 1)) && FileHandler::getTypeByFileName(filename) == FileTypes::MZML
                           && loadMzMLChunked(filename, f.getOptions(), threads, spectra, help, chunk_threads, keep);
      if (chunked)
      {
        ChromatogramTools().convertSpectraToChromatograms<PeakMap>(spectra, true); // as FileHandler::loadExperiment()
        // No malloc_trim(0) here: the pages that the reader threads freed in their arenas (109-165 MB with 16 readers)
        // stay with the process. Returning them cost 5-6 % wall time on 20 small files at 64 threads (int14: 0.744 vs
        // 0.700 of int12-c2, 3 rounds) for 0.20 GiB less max RSS.
      }
      else if (keep && FileHandler::getTypeByFileName(filename) == FileTypes::MZML)
      {
        // streamed, so that only the kept spectra are held
        spectra.clear(true);
        MSDataTransformingConsumer consumer;
        // (moved: the reader does not use a spectrum once its consumer has it)
        consumer.setSpectraProcessingFunc([&spectra, &keep](MSSpectrum& spectrum) { if (keep(spectrum)) { spectra.addSpectrum(std::move(spectrum)); } });
        f.transformExperiment(filename, &consumer, {FileTypes::MZML});
      }
      else
      {
        f.loadExperiment(filename, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW});
        removeUnkept(spectra.getSpectra(), keep); // (Bruker .d, Thermo .raw)
      }
      // by retention time, and the peaks of each spectrum by m/z: the chunks of loadMzMLChunked() sorted the peaks of
      // their spectra on their threads (11 ms on one thread for the 16 million peaks of an Astral run); spectra listed in
      // retention time order (as mzML files usually list them) are not sorted again (see inRetentionTimeOrder())
      if (!chunked || !inRetentionTimeOrder(spectra)) { spectra.sortSpectra(!chunked); }
      return spectra;
    }

    // The threads that a BackgroundConsumer (below) lends from this thread while it runs: the read policies below count
    // them as this thread's (the threads this thread was given), so that a read decides as without the consumer.
    thread_local int consumer_lent_threads = 0;

    // The OpenMP threads of this thread, with those that a BackgroundConsumer has lent
    [[maybe_unused]] int threadsGiven()
    {
#ifdef _OPENMP
      return omp_get_max_threads() + consumer_lent_threads;
#else
      return 1;
#endif
    }

    // The OpenMP threads a spectrum file is read with, of those of this thread: at most 16 (the chunked read gets no
    // faster with more). A read in the foreground gets all of them. A read in the background (BackgroundSpectra) takes
    // its threads from this thread, which helps with the rest of the read when it takes the spectra: next to the index
    // build an eighth of them, at least one (measured at 16 and 64 threads on 8,000-spectrum files of six instrument
    // types, against a reader with a quarter or all of the threads next to all threads of this thread: no more busy
    // threads than these, 12-17% less CPU time, 3% more wall time at 16 threads and 4-5% less at 64; a fixed share
    // without the help: 7-10% more wall time at 16 threads, as either side idles its share once it is done); next to the
    // search of another file a quarter, at least one, which this thread takes back between the phases of the search
    // once the read has ended.
    int readerThreads([[maybe_unused]] bool background, [[maybe_unused]] bool next_to_index_build = true)
    {
#ifdef _OPENMP
      const int threads = threadsGiven();
      if (!background) { return std::min(threads, 16); }
      return std::clamp(threads / (next_to_index_build ? 8 : 4), 1, 16);
#else
      return 1;
#endif
    }

    // The threads whose number sets the chunks of a read in the background (see loadMzMLChunked()): as many as such a
    // read had before it took its threads from those of this thread (a quarter, at least two; or all, at most 16, if
    // the index build waits for it), so that a file is cut into the same chunks as before. The order in which the
    // parsers of the chunks register the names of the meta values they meet first depends on the chunks.
    int chunkThreads([[maybe_unused]] bool awaited)
    {
#ifdef _OPENMP
      const int threads = threadsGiven();
      return awaited ? std::min(threads, 16) : std::clamp(threads / 4, std::min(threads, 2), 16);
#else
      return 1;
#endif
    }

    // The spectra of a file that a helper thread reads (loadMS2Spectra()) with some of the OpenMP threads of this thread
    // (see readerThreads()), while this thread runs its parallel regions with the others: together they keep to the
    // threads they were given (-threads). With all of its threads next to the reader, this thread kept up to twice as
    // many busy (29 cores at 16 threads while the FASTA file was read and the peptides generated and sorted), which
    // costs time wherever the process may not use more cores than it was given (a CPU quota, a shared machine). This
    // thread gets its threads back when it takes the spectra (get()) or no longer waits for them (destruction, start()
    // of the next read). The helper sets its OpenMP threads itself: a new thread starts with the default of the process
    // (OMP_NUM_THREADS or all cores), not with what was set for this one (e.g. by TOPPBase from -threads).
    // Whether a read in the background takes its threads from this thread (see BackgroundSpectra): only if this thread
    // has fewer threads than the process has physical cores (SysInfo::getAvailablePhysicalCores()). With as many threads
    // as cores, or more, the reader runs next to them with threads of its own (chunkThreads()), as before the budget
    // (measured on 64 cores, 20 files: at 64 and 128 threads 1-2% more wall time and 1-7% more CPU time than with the
    // budget; at 16 threads unchanged).
    bool readBorrowsThreads()
    {
#ifdef _OPENMP
      static const int cores = static_cast<int>(SysInfo::getAvailablePhysicalCores());
      return threadsGiven() < cores;
#else
      return false;
#endif
    }

    class BackgroundSpectra;
    thread_local BackgroundSpectra* lending_read = nullptr;  // the read that lends threads of this thread, if any
    thread_local BackgroundSpectra* deferred_read = nullptr; // the read that waits for startDeferredRead(), if any

    class BackgroundSpectra
    {
    public:
      BackgroundSpectra() = default;

      /// @p survey_scans: the MS1 spectra too (chimeric:enabled, see loadMS2Spectra())
      explicit BackgroundSpectra(const std::string& filename, const bool awaited = false, const bool survey_scans = false)
      {
        start(filename, awaited, 0, survey_scans);
      }

      BackgroundSpectra(const BackgroundSpectra&) = delete;
      BackgroundSpectra& operator=(const BackgroundSpectra&) = delete;

      ~BackgroundSpectra()
      {
        if (deferred_read == this) { deferred_read = nullptr; }
        restoreThreads_();
      } // then ~future waits for the helper

      /// Starts to read @p filename, with @p reader_threads (default: readerThreads(true)) of the threads of this thread,
      /// after the read before (if any) has ended and this thread has its threads back. @p awaited: the index build waits
      /// for the spectra (see chunkThreads()). Used by the thread that created it only.
      void start(const std::string& filename, const bool awaited = false, const int reader_threads = 0,
                 const bool survey_scans = false)
      {
        if (deferred_read == this) { deferred_read = nullptr; }
        deferred_ = nullptr;
        restoreThreads_();
        const int chunk_threads = chunkThreads(awaited);
        const bool borrow = readBorrowsThreads();
        const int threads = !borrow ? chunk_threads : reader_threads > 0 ? reader_threads : readerThreads(true);
        spectra_ = std::future<PeakMap>(); // waits for the helper of a read before
        help_ = std::make_shared<ChunkHelp>();
        help_->setLent(borrow);
        spectra_ = std::async(std::launch::async, [&filename, threads, chunk_threads, survey_scans, help = help_]()
        {
          const OnScopeExit done([&help]() { help->close(); }); // get() waits no longer for chunks to parse
#ifdef _OPENMP
          omp_set_num_threads(threads);
#endif
          return loadMS2Spectra(filename, threads, help.get(), chunk_threads, survey_scans);
        });
#ifdef _OPENMP
        if (!borrow) { return; }
        own_threads_ = omp_get_max_threads();
        omp_set_num_threads(std::max(1, own_threads_ - threads));
        // taken back early (takeBackThreadsOfEndedRead()) only outside of parallel regions, where this thread runs a
        // single task, whose ICV this is
        if (omp_get_level() == 0) { lending_read = this; }
#endif
      }

      /// Like start(@p filename, false, readerThreads(true, false)), but only when this thread calls startDeferredRead()
      /// (or get()): during the rescoring of the search before, whose Percolator training uses at most 9 threads.
      void startDeferred(const std::string& filename, const bool survey_scans = false)
      {
#ifdef _OPENMP
        if (omp_get_level() != 0 || !readBorrowsThreads()) { start(filename, false, readerThreads(true, false), survey_scans); return; } // see start()
#endif
        restoreThreads_();
        spectra_ = std::future<PeakMap>(); // waits for the helper of a read before
        help_ = nullptr;
        deferred_ = &filename;
        deferred_survey_scans_ = survey_scans;
        deferred_read = this;
      }

      /// Starts a deferred read (startDeferred())
      void startNow()
      {
        if (deferred_ != nullptr) { start(*deferred_, false, readerThreads(true, false), deferred_survey_scans_); }
      }

      bool valid() const { return spectra_.valid() || deferred_ != nullptr; }

      /// Whether the spectra are read (valid() only)
      bool ready() const { return deferred_ == nullptr && spectra_.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }

      /// Gives this thread its threads back if the spectra are read
      void takeBackThreadsIfRead()
      {
        if (own_threads_ > 0 && ready()) { restoreThreads_(); }
      }

      /// The spectra (valid() only; waits for them; rethrows what the reader threw). Meanwhile, the threads of this
      /// thread that the reader leaves parse chunks of the file with it (see ChunkHelp).
      PeakMap get()
      {
        startNow();
        {
          const OnScopeExit restore([this]() { restoreThreads_(); }); // also if help() throws
          if (help_) { help_->help(); }
        }
        return spectra_.get();
      }

    private:
      void restoreThreads_() noexcept
      {
#ifdef _OPENMP
        if (own_threads_ > 0) { omp_set_num_threads(own_threads_); }
#endif
        own_threads_ = 0;
        if (lending_read == this) { lending_read = nullptr; }
      }

      std::shared_ptr<ChunkHelp> help_;
      std::future<PeakMap> spectra_;
      int own_threads_ = 0; // of this thread, while the helper reads; 0: not reduced
      const std::string* deferred_ = nullptr; // the file of a deferred read
      bool deferred_survey_scans_ = false;    // ... read with its MS1 spectra
    };

    // The consumer of the multi-file search's results (ProSEAlgorithm::setFileResultConsumer()), one result at a time, on a
    // helper thread with two of the OpenMP threads of this thread (one if it has fewer than 8), which lends them (lowers
    // its own) until it takes them back where it starts a phase after the helper has ended
    // (takeBackThreadsOfEndedRead()), when it starts the next consumer, or at the end (finish(), destruction). Only if a read in the background never takes threads from
    // this thread (readBorrowsThreads() false): this thread's number of threads is then changed here only, so it can be
    // set back as it was. Otherwise, or with a single thread, or outside of the top level, the consumer runs here.
    class BackgroundConsumer;
    thread_local BackgroundConsumer* lending_consumer = nullptr; // the consumer that lends threads of this thread, if any

    class BackgroundConsumer
    {
    public:
      BackgroundConsumer() = default;
      BackgroundConsumer(const BackgroundConsumer&) = delete;
      BackgroundConsumer& operator=(const BackgroundConsumer&) = delete;
      ~BackgroundConsumer() { finish(); }

      /// Calls @p consumer(@p index, @p result) on a helper thread or here (see above), after the consumer before.
      void start(const ProSEAlgorithm::FileResultConsumer& consumer, const Size index, const ProSEAlgorithm::SearchResult& result)
      {
        finish();
        int helper_threads = 0;
#ifdef _OPENMP
        const int threads = omp_get_max_threads();
        // Two threads, not a quarter: the consumer (the idXML output of ProSE) has the search of the next file to run
        // next to, and with glibc's malloc (default arena settings) each further thread of its team can populate an
        // arena of its own, whose freed memory may remain resident. Search of 6 timsTOF HeLa runs in one process at
        // 64 threads (kim, 64 cores; 16 / 2 helper threads):
        // max RSS 25.9 / 24.5 GiB, wall time 38.2 / 37.2 s, identical output (the idXML of the largest run written in
        // about 2.5 s instead of 1.5 s, still during the search of the next one). The next consumer waits for the one
        // before, so searches shorter than that write wait for it: that run followed by two small ones (0.4 GB) took
        // 17.4 instead of 15.4 s (4 / 8 helper threads: 15.9 / 15.6 s, and less of the memory saved).
        if (omp_get_level() == 0 && threads > 1 && !readBorrowsThreads()) { helper_threads = std::min(2, std::max(1, threads / 4)); }
#endif
        if (helper_threads > 0)
        {
          try
          {
            done_ = false;
            helper_ = std::thread([this, &consumer, index, &result, helper_threads]()
            {
#ifdef _OPENMP
              omp_set_num_threads(helper_threads);
#endif
              try { consumer(index, result); } catch (...) {} // (the consumer handles its errors)
              done_.store(true, std::memory_order_release);
            });
#ifdef _OPENMP
            own_threads_ = threads;
            omp_set_num_threads(threads - helper_threads);
            consumer_lent_threads = helper_threads;
            lending_consumer = this;
#endif
            return;
          }
          catch (...) {} // no thread: here
        }
        try { consumer(index, result); } catch (...) {}
      }

      /// Gives this thread its threads back if the consumer has ended (at the top level only)
      void takeBackThreadsIfDone()
      {
        if (helper_.joinable() && done_.load(std::memory_order_acquire)) { finish(); }
      }

      /// Waits for the consumer and gives this thread its threads back
      void finish() noexcept
      {
        if (helper_.joinable()) { helper_.join(); }
#ifdef _OPENMP
        if (own_threads_ > 0) { omp_set_num_threads(own_threads_); consumer_lent_threads = 0; }
#endif
        own_threads_ = 0;
        if (lending_consumer == this) { lending_consumer = nullptr; }
      }

    private:
      std::thread helper_;
      std::atomic<bool> done_{false};
      int own_threads_ = 0; // of this thread, while the helper runs; 0: not lowered
    };

    void takeBackThreadsOfEndedRead()
    {
#ifdef _OPENMP
      if (lending_read != nullptr && omp_get_level() == 0) { lending_read->takeBackThreadsIfRead(); } // see start()
      if (lending_consumer != nullptr && omp_get_level() == 0) { lending_consumer->takeBackThreadsIfDone(); } // see BackgroundConsumer
#endif
    }

    void startDeferredRead()
    {
      if (deferred_read != nullptr) { deferred_read->startNow(); }
    }

    bool readLendsThreads() { return lending_read != nullptr; }

    // Whether the first MBs of an mzML file name an electron-based activation (the terms that
    // MzMLHandler reads as ECD, ETD, ETciD or EThcD). A cheap prediction of what the spectra will
    // hold, used only to decide when to build the index, not what the index holds.
    bool mzMLHeadNamesElectronActivation(const std::string& filename)
    {
      std::string head(Size(8) << 20, '\0');
      std::ifstream in(filename, std::ios::binary);
      in.read(head.data(), static_cast<std::streamsize>(head.size()));
      head.resize(static_cast<Size>(in.gcount()));
      // MS:1000250, MS:1000598, MS:1003182 and MS:1002631: one pass for what they share
      const std::string_view text(head);
      for (Size pos = text.find("MS:100"); pos != std::string_view::npos; pos = text.find("MS:100", pos + 1))
      {
        const std::string_view number = text.substr(pos + 6, 4);
        if (number == "0250" || number == "0598" || number == "3182" || number == "2631")
        {
          return true;
        }
      }
      return false;
    }
  }

  PeakMap ProSEAlgorithm::loadSpectra(const std::string& filename, const std::function<bool(const MSSpectrum&)>& keep)
  {
    return loadMS2Spectra(filename, readerThreads(false), nullptr, 0, false, keep);
  }

  // =====================================================================
  // File-based search: thin I/O wrapper that delegates to in-memory search
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      const std::string& in_spectra, const std::string& in_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    reset_run_name_ = in_spectra; // names the RESET stem table of the run (reset:stem_table a directory)
    // load MS2 map
    PeakMap spectra;

    vector<FASTAFile::FASTAEntry> fasta_db;
    DecoyStrategy_ strategy; // decoys of the searched database, for protein FDR below
    bool strategy_resolved = false;
    ExitCodes ec;
    const Size threads = startOpenMPThreads();
    if (database_chunk_size_ == 0 && FileHandler::getTypeByFileName(in_spectra) == FileTypes::MZML && threads > 1)
    {
      // Unchunked, multi-threaded search of an mzML file: a helper thread reads the spectra while
      // this thread reads the FASTA file and builds the fragment index; then the search continues
      // as search(spectra, fasta_db, ...) does. The two sides share no data: the spectra are used
      // here only after get(), and the registries of OpenMS are written only by the mzML reader
      // in the meantime (FASTA reading, decoys and the index build register no meta value names),
      // so they end up as after reading the spectra first. Other readers (.d, .raw) have not
      // been checked for this and keep the order below.
      // The OpenMP threads of this thread run before the helper starts (see the condition), as
      // they did when the spectra were read first (decoding them was the first parallel region):
      // started next to a busy helper, they more often end up on another NUMA node than this
      // thread, which slows down the index build and the scoring.
      BackgroundSpectra spectra_ready(in_spectra, false, chimericSearch_());

      // load FASTA
      loadFASTA(in_db, fasta_db);

      // ions:by_activation: the index needs c and z+1 ions if electron-activated spectra are
      // searched, which is known once the spectra are read. Wait for them if they are read already
      // (or failed to be read: get() rethrows) or if the file names an electron-based activation
      // early on. Otherwise build the index without these ions meanwhile, which is what most data
      // need, and rebuild it should such spectra turn up after all (as the multi-file search does
      // when a later file has them): the index searched is the one prepareContext(fasta_db, true)
      // builds.
      const bool spectra_read = spectra_ready.ready()
                                || (ions_by_activation_ && mzMLHeadNamesElectronActivation(in_spectra));
      if (spectra_read) { spectra = spectra_ready.get(); }
      bool spectra_waited = spectra_read;
      // The single-use index holds only the peptides the spectra can reach if they are read by the time the
      // peptides are generated, or if waiting for them pays off (see MAX_THREADS_WAITING_FOR_SPECTRA).
      std::function<const PeakMap*(Size)> searched_spectra;
      if (restrictIndexToSpectra_())
      {
        const UInt64 spectra_bytes = File::fileSize(in_spectra); // UInt64(-1) if unknown: no waiting
        searched_spectra = [&, spectra_bytes](Size peptides) -> const PeakMap*
        {
          if (!spectra_waited)
          {
            const bool wait = threads <= MAX_THREADS_WAITING_FOR_SPECTRA
                              && spectra_bytes <= MAX_SPECTRA_BYTES_PER_PEPTIDE * peptides / threads;
            if (!wait && !spectra_ready.ready())
            {
              return nullptr;
            }
            spectra = spectra_ready.get();
            spectra_waited = true;
          }
          return &spectra;
        };
      }
      // The context takes the entries over instead of copying them. fasta_db is not read
      // afterwards: protein FDR below takes the decoy marker from the context, which holds what
      // resolveDecoyStrategy_(fasta_db) returned.
      SearchContext ctx = prepareContext_(std::move(fasta_db), spectra_read && countElectronActivated_(spectra) > 0, searched_spectra);
      strategy.have_decoys = ctx.have_decoys;
      strategy.decoy_string = ctx.decoy_string;
      strategy.is_prefix = ctx.decoy_is_prefix;
      strategy_resolved = true;
      if (!spectra_read)
      {
        if (!spectra_waited)
        {
          spectra = spectra_ready.get();
          spectra_waited = true;
        }
        if (countElectronActivated_(spectra) > 0)
        {
          startProgress(0, 1, "Building fragment index with c and z+1 ions...");
          ctx.fragment_index.clear();
          ctx.fragment_index.setParameters(fragmentIndexParameters_(true));
          ctx.fragment_index.build(ctx.db, searched_spectra);
          ctx.electron_ions = true;
          endProgress();
        }
      }
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      ec = search(spectra, ctx, protein_ids, peptide_ids);
    }
    else
    {
      spectra = loadMS2Spectra(in_spectra, readerThreads(false), nullptr, 0, chimericSearch_());

      // load FASTA
      loadFASTA(in_db, fasta_db);

      // delegate to in-memory search
      ec = search(spectra, fasta_db, protein_ids, peptide_ids);
    }

    if (ec != ExitCodes::EXECUTION_OK)
    {
      return ec;
    }

    // Protein inference + picked-protein FDR for single-file search.
    // Must run before decoy removal so both target and decoy proteins
    // receive aggregated scores from BPIA. Resolve the decoy strategy from the
    // same input FASTA the search used so the marker/position match.
    // The strategy is only needed for protein FDR: do not scan the accessions again without it.
    if (fdr_protein_ > 0.0)
    {
      if (!strategy_resolved) { strategy = resolveDecoyStrategy_(fasta_db); } // else: known from the context
      if (strategy.have_decoys)
      {
        // Single input file = complete experiment, so picked-protein FDR is valid. Use the resolved
        // decoy marker (prefix or suffix, detected by DecoyHelper in resolveDecoyStrategy_) so the
        // shared finalization recognises the same decoys that were searched.
        applyCompleteSetProteinFDR(protein_ids, peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_,
                                   true, inference_only_peptide_ids_);
      }
    }
    inference_only_peptide_ids_ = PeptideIdentificationList();

    // patch file-specific metadata
    protein_ids[0].getSearchParameters().db = in_db;
    protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // In-memory searchWithModificationAnalysis
  // =====================================================================
  ProSEAlgorithm::SearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      PeakMap& spectra,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const std::string& output_base_name) const
  {
    SearchResult result;
    result.is_open_search = isOpenSearchMode_();

    result.exit_code = search(spectra, fasta_db, result.protein_ids, result.peptide_ids);

    // Carry whatever statistics the search managed to populate, even on failure.
    result.stats = last_run_stats_;
    result.percolator = last_percolator_;
    result.inference_only_peptide_ids = std::move(inference_only_peptide_ids_);
    inference_only_peptide_ids_ = PeptideIdentificationList();

    if (result.exit_code != ExitCodes::EXECUTION_OK)
    {
      return result;
    }

    if (result.is_open_search)
    {
      OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for open search results..." << std::endl;

      OpenSearchModificationAnalysis mod_analyzer;

      std::string output_file;
      if (!output_base_name.empty())
      {
        output_file = output_base_name + "_ModificationAnalysis.idXML";
      }

      // Read the post-calibration value captured during the internal search()
      // call. Do NOT re-compute here — by now search()'s restore_fi_params() has
      // reset the tolerance members to user-configured values.
      result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
        result.peptide_ids,
        last_mod_match_tolerance_used_,
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        output_file
      );
    }
    else
    {
      OPENMS_LOG_INFO << "[ProSE] Closed search mode - modification analysis skipped" << std::endl;
    }

    // A single input is a complete experiment: protein FDR as in search(file).
    if (fdr_protein_ > 0.0 && !result.protein_ids.empty() && !result.peptide_ids.empty())
    {
      const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
      if (strategy.have_decoys)
      {
        applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_,
                                   true, result.inference_only_peptide_ids);
        updateFinalStats(result.stats, result.peptide_ids, enzyme_, result.stats.fdr_applied);
      }
    }
    result.inference_only_peptide_ids = PeptideIdentificationList(); // used by the protein FDR only

    return result;
  }

  // =====================================================================
  // Multi-file searchWithModificationAnalysis (in-memory FASTA).
  //
  // Builds the FragmentIndex once via prepareContext() and reuses it across
  // all input spectrum files. Each input file produces its own SearchResult
  // (with per-file modification analysis); a final aggregate SearchResult is
  // computed by pooling all per-file PSMs and running modification analysis
  // once on the pooled set. build_pooled_aggregate == false suppresses keeping
  // that pooled PSM copy (see the header for the exact contract).
  // =====================================================================
  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      const std::vector<std::string>& in_spectra_files,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    return searchWithModificationAnalysis_(in_spectra_files, fasta_db, nullptr, output_base_names, aggregate_base_name,
                                           build_pooled_aggregate);
  }

  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis_(
      const std::vector<std::string>& in_spectra_files,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db_given,
      const std::string* in_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    // The database: fasta_db_given, or read from *in_db below, once the first spectrum file is read in the background.
    std::vector<FASTAFile::FASTAEntry> fasta_db_read;
    const std::vector<FASTAFile::FASTAEntry>& fasta_db = in_db == nullptr ? fasta_db_given : fasta_db_read;
    if (!output_base_names.empty() && output_base_names.size() != in_spectra_files.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "output_base_names must be empty or have exactly one entry per spectrum file (got "
        + StringUtils::toStr(output_base_names.size()) + " entries for " + StringUtils::toStr(in_spectra_files.size())
        + " spectrum files).");
    }

    MultiFileSearchResult mfres;
    mfres.aggregate.is_open_search = isOpenSearchMode_();
    released_hits_.reset(); // the releases of a search before (see the end of the loop below), so that they do not pile up
    // 'rescore' with several files: keep each file's PSMs before the rescoring until all files are searched, for
    // rescoreFilesTogether_() (in case Percolator cannot rescore some file on its own)
    pool_rescoring_ = rescore_ && in_spectra_files.size() > 1;
    const OnScopeExit end_pool_rescoring([this]() { pool_rescoring_ = false; });

    if (in_spectra_files.empty())
    {
      if (in_db != nullptr) { loadFASTA(*in_db, fasta_db_read); } // reports a database that cannot be read, as before
      mfres.aggregate.exit_code = ExitCodes::INPUT_FILE_EMPTY;
      return mfres;
    }

    // Multi-threaded: a helper thread reads the spectra of an mzML file while this thread works, the
    // first file from here on while the database is prepared and the index built, and, in the
    // unchunked search below, file i + 1 while file i is searched. The OpenMP threads of this thread
    // are started first, as in search(file).
    // Registry order: idXML writes the UserParams of an object in the order in which their names were
    // registered. The mzML reader registers names (CV terms and user parameters of the file, and its
    // own, see MzMLHandler) while it runs. Preparing the database and building the index register
    // none, so with the first file the registry ends up as if the file had been read first. While
    // file i is searched, the search registers the names of what it writes (scan_index, the PSM
    // features, target_decoy, PeptideIndexer:*, spectra_data, ...); if the reader of file i + 1
    // registers names in between, these keep their order relative to each other. ProSE writes none
    // of the reader's names, unless a file carries a user parameter that is named like one of them.
    const Size threads_started = startOpenMPThreads();
    const bool read_in_background = threads_started > 1;
    const auto is_mzml = [&in_spectra_files](Size i) { return FileHandler::getTypeByFileName(in_spectra_files[i]) == FileTypes::MZML; };
    BackgroundSpectra next_spectra; // the spectra of the next file to search, if read in the background
    // A single file is the only search of the context (r3 merge of WP7 and WP10): as in search(file), the index
    // holds only the peptides its spectra can reach if they are read by the time the peptides are generated, or
    // if waiting for them pays off (see MAX_THREADS_WAITING_FOR_SPECTRA); the spectra taken for that are searched.
    const bool restrict_index = in_spectra_files.size() == 1 && restrictIndexToSpectra_();
    // With more threads than MAX_THREADS_WAITING_FOR_SPECTRA, the index build waits for the spectra of such a file
    // (and helps to read them, see BackgroundSpectra::get()).
    const bool wait_for_first_spectra = restrict_index && threads_started > MAX_THREADS_WAITING_FOR_SPECTRA;
    const auto read_next = [&](Size i)
    {
      if (read_in_background && i < in_spectra_files.size() && is_mzml(i))
      {
        // File 0 next to the index build. A later one next to the search of the file before: from its rescoring on if
        // Percolator, which trains on at most 9 threads, leaves the reader's threads (a quarter) idle; otherwise from the
        // start (measured on three Orbitrap files at 16 threads without rescoring: from the end of the scoring on, 3%
        // slower, as more of the read is left for get()). The reader then registers the names of the file's meta values
        // after the protein mapping of the file before (see the registry order above).
        if (i == 0) { next_spectra.start(in_spectra_files[i], wait_for_first_spectra, 0, chimericSearch_()); }
        else if (rescore_ && threads_started >= 12) { next_spectra.startDeferred(in_spectra_files[i], chimericSearch_()); }
        else { next_spectra.start(in_spectra_files[i], false, readerThreads(true, false), chimericSearch_()); }
      }
    };
    read_next(0);
    // Whether the first file names an electron-based activation early on, if build_index_while_first_file_is_read()
    // below asks: looked up meanwhile, from here on.
    std::future<bool> first_file_names_electron_activation;
    if (next_spectra.valid() && ions_by_activation_)
    {
      first_file_names_electron_activation = std::async(std::launch::async | std::launch::deferred,
        [&in_spectra_files]() { return mzMLHeadNamesElectronActivation(in_spectra_files[0]); });
    }
    // A database read here is read while the first file is: neither registers meta value names (see above).
    if (in_db != nullptr) { loadFASTA(*in_db, fasta_db_read); }

    // Resolve decoy handling once from the shared input FASTA; reused for the
    // chunk-major path, the single-context path, and the downstream FDR steps.
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    // Surface the effective marker so a caller doing merged-PSM protein FDR
    // (e.g. ProSE -out_merged) recognises the same decoys that were searched.
    mfres.decoy_string = strategy.decoy_string;
    mfres.decoy_is_prefix = strategy.is_prefix;
    mfres.have_decoys = strategy.have_decoys;

    // -- Shared configuration recap for the end-of-search report (built once) --
    {
      SharedSearchStats& sh = mfres.shared;
      sh.enzyme = enzyme_;
      sh.precursor_tol_lower = precursor_mass_tolerance_lower_;
      sh.precursor_tol_upper = precursor_mass_tolerance_upper_;
      sh.precursor_tol_unit = precursor_mass_tolerance_unit_;
      sh.fragment_tol = fragment_mass_tolerance_;
      sh.fragment_tol_unit = fragment_mass_tolerance_unit_;
      sh.min_charge = static_cast<Int>(precursor_min_charge_);
      sh.max_charge = static_cast<Int>(precursor_max_charge_);
      sh.missed_cleavages = peptide_missed_cleavages_;
      sh.fixed_mods.assign(modifications_fixed_.begin(), modifications_fixed_.end());
      sh.variable_mods.assign(modifications_variable_.begin(), modifications_variable_.end());
      if (add_a_ions_) sh.ion_series.push_back("a");
      if (add_b_ions_) sh.ion_series.push_back("b");
      if (add_c_ions_) sh.ion_series.push_back("c");
      if (add_x_ions_) sh.ion_series.push_back("x");
      if (add_y_ions_) sh.ion_series.push_back("y");
      if (add_z_ions_) sh.ion_series.push_back("z");
      if (add_zp1_ions_) sh.ion_series.push_back("z+1");
      if (ions_by_activation_) sh.ion_series.push_back("+c/z+1 for ETD/ECD/EThcD/ETciD");
      sh.open_search = isOpenSearchMode_();
      sh.calibration_enabled = calibration_enabled_;
      sh.psm_fdr_threshold = fdr_psm_;
      sh.protein_fdr_threshold = fdr_protein_;
      // Decoy handling resolved via the auto/generate/ignore strategy (#9634):
      // target-only when no decoys are present, "generated" when ProSE synthesises
      // them, "external" when pre-existing decoys in the FASTA are reused.
      sh.decoy_mode = !strategy.have_decoys ? "none (target-only)"
                      : (strategy.generate ? "generated" : "external");
    }

    // Decide chunking on the decoy-augmented DB size (#9180): if decoy generation
    // doubles the target DB, a 3000-protein target against chunk_size=5000 should
    // still chunk because the augmented DB is 6000 — otherwise the resulting FI
    // would exceed the user's memory budget by 2×.
    // A database read from *in_db is not used after it is augmented (at most once, here or in prepare_context()
    // below): its entries are taken over instead of being copied.
    const auto augment_db = [&]()
    {
      return in_db == nullptr ? buildDecoyAugmentedDB_(fasta_db, strategy) : buildDecoyAugmentedDB_(std::move(fasta_db_read), strategy);
    };
    std::vector<FASTAFile::FASTAEntry> full_db;
    bool use_chunked = false;
    if (database_chunk_size_ > 0)
    {
      full_db = augment_db();
      use_chunked = (full_db.size() > database_chunk_size_);
    }

    if (use_chunked)
    {
      // ================================================================
      // Chunk-major multi-file path (MSFragger-style): build each chunk's
      // FragmentIndex ONCE and score ALL files against it before moving
      // to the next chunk. This gives C FI builds instead of N×C.
      // No file is searched on its own here, so there is no "Searching <file>" progress step per file as in the
      // non-chunked path below; the steps of the phases and chunks are reported as before.
      // ================================================================
      // full_db already built above.
      if (rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static"))
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode " << reuse_mode_ << " needs an unchunked search: every file is "
                        << "rescored by Percolator trained on it ('per_file')." << std::endl;
      }
      const bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
      const bool open_search_mode = isOpenSearchMode_();

      // Shared report stats: chunk-major builds C indices (one per chunk),
      // reused across all files — count/time them once, not per file.
      mfres.shared.chunked = true;
      // Count by the RESOLVED decoy marker (prefix OR suffix), not the hardcoded decoy_prefix_: otherwise reused
      // external/suffix decoys (decoy_mode "external") would be miscounted as targets. have_decoys is false for
      // target-only (ignore), where decoys are stripped and decoy_string is empty.
      mfres.shared.db_decoy_proteins = strategy.have_decoys ? countDecoyEntries_(full_db, strategy.decoy_string, strategy.is_prefix) : 0;
      mfres.shared.db_target_proteins = full_db.size() - mfres.shared.db_decoy_proteins;

      OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search_mode ? "true" : "false")
                      << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                      << ", +" << precursor_mass_tolerance_upper_ << "] "
                      << precursor_mass_tolerance_unit_ << ")" << std::endl;

      // Phase 1: Load + preprocess all files.
      std::vector<PeakMap> all_spectra(in_spectra_files.size());
      // annotate:self_trained_ion_priors: the compact peak lists of every file, until its PSMs are annotated
      std::vector<FragmentIonLikelihoodModel::PeakLists> all_ion_evidence(in_spectra_files.size());
      const bool retain_evidence = param_.getValue("annotate:local_fragment_evidence").toBool();
      std::vector<PeakMap> all_evidence_spectra(retain_evidence ? in_spectra_files.size() : 0);
      std::vector<PeakMap> all_query_spectra(query_raw_spectrum_ ? in_spectra_files.size() : 0);
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        OPENMS_LOG_INFO << "[ProSE] Loading " << in_spectra_files[i] << std::endl;
        all_spectra[i] = next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra_files[i], readerThreads(false), nullptr, 0, chimericSearch_());
        preprocessSearchedSpectra_(all_spectra[i], in_spectra_files[i], all_ion_evidence[i],
                                   retain_evidence ? &all_evidence_spectra[i] : nullptr,
                                   query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
      }

      // ions:by_activation: the chunk indices, shared by all files, hold c and z+1 ions if any file has
      // electron-activated spectra. Only those spectra are matched against them (see
      // scoreSpectraAgainstIndex_), so the other files get the results they would get alone.
      bool electron_ions = false;
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const Size n_electron_activated = countElectronActivated_(all_spectra[i]);
        if (n_electron_activated == 0) continue;
        electron_ions = true;
        OPENMS_LOG_INFO << "[ProSE] " << in_spectra_files[i] << ": " << n_electron_activated << " of "
                        << all_spectra[i].size() << " spectra are electron-activated (ETD, ECD, EThcD or ETciD)"
                        << " and are also scored with c and z+1 ions." << std::endl;
      }

      // Per-file calibration: build a strided calibration FI once, run calibration per file.
      // Stores per-file effective tolerances (asymmetric lower/upper preserved — see #9180)
      // for use during scoring.
      const Size chunk_size = database_chunk_size_;
      const Size n_chunks = (full_db.size() + chunk_size - 1) / chunk_size;
      struct PerFileCalibration
      {
        double effective_precursor_tol_lower;
        double effective_precursor_tol_upper;
        double effective_fragment_tol;
        double mod_match_tol;  // for open-search mod analysis
      };
      std::vector<PerFileCalibration> per_file_cal(in_spectra_files.size());

      // Default: user-configured tolerances.
      for (auto& cal : per_file_cal)
      {
        cal.effective_precursor_tol_lower = precursor_mass_tolerance_lower_;
        cal.effective_precursor_tol_upper = precursor_mass_tolerance_upper_;
        cal.effective_fragment_tol = fragment_mass_tolerance_;
        cal.mod_match_tol = computeModMatchTolerance_();
      }

      if (calibration_enabled_ && !open_search_mode)
      {
        // Build a strided-sample calibration FI once, reused across files.
        std::vector<FASTAFile::FASTAEntry> cal_db = buildCalibrationSample_(full_db);
        FragmentIndex cal_fi;
        cal_fi.setParameters(fragmentIndexParameters_(electron_ions));
        StopWatch sw_cal_idx; sw_cal_idx.start();
        cal_fi.build(cal_db);
        sw_cal_idx.stop();
        mfres.shared.seconds_index_build += sw_cal_idx.getClockTime();

        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
          OPENMS_LOG_INFO << "[ProSE] Calibration for " << in_spectra_files[i]
                          << " (strided sample, " << cal_db.size() << " proteins)" << std::endl;
          CalibrationResult_ cal = runCalibrationPass_(all_spectra[i], cal_fi, cal_db, query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
          if (cal.success)
          {
            per_file_cal[i].effective_fragment_tol = cal.fragment_tolerance;
            if (!cal.extreme_bias)
            {
              per_file_cal[i].effective_precursor_tol_lower = cal.cal_lower;
              per_file_cal[i].effective_precursor_tol_upper = cal.cal_upper;
              OPENMS_LOG_INFO << "[ProSE] Calibration: shift=" << cal.precursor_shift
                              << " " << precursor_mass_tolerance_unit_
                              << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                              << " fragment=" << cal.fragment_tolerance << std::endl;
            }
            else
            {
              OPENMS_LOG_WARN << "[ProSE] Calibration for " << in_spectra_files[i]
                              << ": extreme bias, precursor calibration discarded. Fragment calibration applied." << std::endl;
            }
            // Recompute mod-match tolerance with calibrated values.
            // Temporarily set member variables, compute, then restore.
            const double orig_lower = precursor_mass_tolerance_lower_;
            const double orig_upper = precursor_mass_tolerance_upper_;
            if (!cal.extreme_bias)
            {
              precursor_mass_tolerance_lower_ = cal.cal_lower;
              precursor_mass_tolerance_upper_ = cal.cal_upper;
            }
            per_file_cal[i].mod_match_tol = computeModMatchTolerance_();
            precursor_mass_tolerance_lower_ = orig_lower;
            precursor_mass_tolerance_upper_ = orig_upper;
          }
          else
          {
            OPENMS_LOG_INFO << "[ProSE] Calibration failed for " << in_spectra_files[i]
                            << ", using configured tolerances." << std::endl;
          }
        }
        // cal_fi freed here.
      }
      else if (calibration_enabled_ && open_search_mode)
      {
        OPENMS_LOG_WARN << "Warning: calibration not applicable in open-search mode." << std::endl;
      }

      // Prepare spectrum generators (once).
      const SpectrumGenerators_ generators = spectrumGenerators_();

      // Per-file hit accumulators.
      std::vector<std::vector<std::vector<AnnotatedHit_>>> per_file_hits(in_spectra_files.size());
      // One pool summary per spectrum per file; accumulates across chunks.
      std::vector<std::vector<CandidatePoolStats_>> per_file_pool_stats(in_spectra_files.size());
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        per_file_hits[i].resize(all_spectra[i].size());
        per_file_pool_stats[i].resize(all_spectra[i].size());
      }
      if (chimeric_enabled_)
      {
        OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
      }
      OPENMS_LOG_INFO << "[ProSE] Chunk-major multi-file: " << full_db.size()
                      << " proteins, " << n_chunks << " chunks, "
                      << in_spectra_files.size() << " files." << std::endl;

      Size chunk_idx = 0;
      for (Size start = 0; start < full_db.size(); start += chunk_size)
      {
        ++chunk_idx;
        const Size end = std::min(start + chunk_size, full_db.size());
        OPENMS_LOG_INFO << "[ProSE] Chunk " << chunk_idx << "/" << n_chunks
                        << " (" << (end - start) << " proteins)" << std::endl;

        std::vector<FASTAFile::FASTAEntry> chunk_db(full_db.begin() + start, full_db.begin() + end);
        FragmentIndex chunk_fi;
        chunk_fi.setParameters(fragmentIndexParameters_(electron_ions));
        StopWatch sw_chunk; sw_chunk.start();
        chunk_fi.build(chunk_db);
        sw_chunk.stop();
        mfres.shared.seconds_index_build += sw_chunk.getClockTime();
        mfres.shared.indexed_peptides += chunk_fi.getPeptides().size();
        mfres.shared.indexed_fragments += chunk_fi.getNumFragments();
        if (chunk_fi.isPeptideMassMode()) { mfres.shared.peptide_mass_mode = true; }

        // Score ALL files against this chunk's index.
        // Each file may have different calibrated tolerances — apply per-file
        // precursor bounds to the FI before scoring, then use per-file
        // fragment tolerance for HyperScore.
        const Param base_fi_params = chunk_fi.getParameters();
        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
          // Apply per-file calibrated precursor bounds to FI query params. Asymmetric
          // lower/upper preserved — collapsing to max() would re-open the tight side
          // of the calibrated window and admit spurious decoy candidates (#9180).
          if (calibration_enabled_ && !open_search_mode)
          {
            Param fi_params = base_fi_params;
            fi_params.setValue("fragment:mass_tolerance", per_file_cal[i].effective_fragment_tol);
            fi_params.setValue("precursor:mass_tolerance_lower", per_file_cal[i].effective_precursor_tol_lower);
            fi_params.setValue("precursor:mass_tolerance_upper", per_file_cal[i].effective_precursor_tol_upper);
            chunk_fi.setParameters(fi_params);
          }
          scoreSpectraAgainstIndex_(all_spectra[i], chunk_fi, chunk_db,
                                    generators, per_file_cal[i].effective_fragment_tol,
                                    fragment_mass_tolerance_unit_ppm, open_search_mode,
                                    per_file_hits[i], per_file_pool_stats[i],
                                    "  file " + StringUtils::toStr(i + 1) + " chunk " + StringUtils::toStr(chunk_idx),
                                    query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
        }
        // Restore base FI params for next chunk (in case calibration modified them).
        if (calibration_enabled_ && !open_search_mode)
          chunk_fi.setParameters(base_fi_params);

        // Per-chunk pruning for each file.
        const Size keep = std::max(searchTopHits_(), Size(2));
        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
#pragma omp parallel for default(none) shared(per_file_hits, i, keep)
          for (SignedSize si = 0; si < (SignedSize)per_file_hits[i].size(); ++si)
          {
            if (per_file_hits[i][si].size() > keep)
            {
              std::partial_sort(per_file_hits[i][si].begin(),
                                per_file_hits[i][si].begin() + keep,
                                per_file_hits[i][si].end(),
                                AnnotatedHit_::hasBetterScore);
              per_file_hits[i][si].resize(keep);
            }
          }
        }
      } // end chunk loop

      // Phase 3: Per-file postprocess with per-file calibrated tolerances.
      for (auto& file_stats : per_file_pool_stats)
      {
        for (auto& stats : file_stats)
        {
          std::unordered_set<std::string>().swap(stats.seen_candidates);
        }
      }
      mfres.per_file.reserve(in_spectra_files.size());

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");
        last_mod_match_tolerance_used_ = per_file_cal[i].mod_match_tol;

        SearchResult result;
        result.is_open_search = open_search_mode;

        if (query_raw_spectrum_) all_query_spectra[i].clear(true);

        postProcessHits_(all_spectra[i], per_file_hits[i], per_file_pool_stats[i],
          result.protein_ids, result.peptide_ids,
          searchTopHits_(), modifications_fixed_, modifications_variable_,
          peptide_missed_cleavages_,
          std::max(per_file_cal[i].effective_precursor_tol_lower,
                   per_file_cal[i].effective_precursor_tol_upper),
          per_file_cal[i].effective_fragment_tol,
          precursor_mass_tolerance_unit_, fragment_mass_tolerance_unit_,
          precursor_min_charge_, precursor_max_charge_, enzyme_, "",
          retain_evidence ? &all_evidence_spectra[i] : nullptr);
        if (retain_evidence) all_evidence_spectra[i].clear(true);

        // Per input file, while each file is still its own identification run. Moving this
        // after the merge below would leave one run to calibrate NCE and RT on, which are
        // per-run quantities; see the note on PeptDeepRescoring::annotate().
        annotatePeptDeepFeatures_(all_spectra[i], result.protein_ids, result.peptide_ids);

        PeptideIndexing indexer;
        Param param_pi = indexer.getParameters();
        param_pi.setValue("decoy_string", strategy.decoy_string.empty() ? decoy_prefix_ : strategy.decoy_string);
        param_pi.setValue("decoy_string_position", strategy.is_prefix ? "prefix" : "suffix");
        param_pi.setValue("enzyme:name", enzyme_);
        param_pi.setValue("enzyme:specificity",
                          EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
        param_pi.setValue("missing_decoy_action", "silent");
        indexer.setParameters(param_pi);
        PeptideIndexing::ExitCodes indexer_exit =
            indexer.run(full_db, result.protein_ids, result.peptide_ids);
        if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
            (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
        {
          OPENMS_LOG_WARN << "[ProSE] PeptideIndexing failed for " << in_spectra
                          << " (exit code " << static_cast<int>(indexer_exit)
                          << "). Skipping this file." << std::endl;
          if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
            result.exit_code = ExitCodes::INPUT_FILE_EMPTY;
          else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
            result.exit_code = ExitCodes::UNEXPECTED_RESULT;
          else
            result.exit_code = ExitCodes::UNKNOWN_ERROR;
          result.stats.input_file = File::basename(in_spectra);
          mfres.per_file.push_back(std::move(result));
          all_ion_evidence[i].clear();
          continue;
        }

        // PSM-level FDR only (matching non-chunked search semantics). Target+decoy PSMs are
        // filtered alike by the q-value threshold; no decoy-specific stripping happens here
        // (decoupled from decoy removal). The decoys that pass are kept because per-file results
        // feed cross-file merging and protein-level FDR. Categorical decoy removal is a
        // protein-FDR finalization step performed by ProSE.cpp. have_decoys is the marker-aware
        // result from resolveDecoyStrategy_ (prefix or suffix, generated or external).
        const bool has_decoys = strategy.have_decoys;
        // Optional per-run ion priors, after PeptideIndexing and before FDR (see search()); each file learns its own.
        if (self_trained_ion_priors_)
        {
          annotateIonPriors_(all_spectra[i], all_ion_evidence[i], result.protein_ids, result.peptide_ids);
          all_ion_evidence[i].clear();
        }
        // Pre-FDR stats (target/decoy counts + HyperScore distribution).
        capturePreFdrStats_(result.peptide_ids, result.stats);
        // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
        reset_run_name_ = in_spectra;
        result.percolator = applyPsmFDR_(result.protein_ids, result.peptide_ids, has_decoys, result.stats,
                                         keepPreRescoringPsms_() ? &result.pre_fdr_peptide_ids : nullptr,
                                         &result.inference_only_peptide_ids);

        result.exit_code = ExitCodes::EXECUTION_OK;

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, all_spectra[i]);

        // In chunk-major mode scoring is shared across files (one pass per chunk),
        // so per-file seconds_search is not separable; it is left 0 and the shared
        // index/total timing carries the cost.
        collectRunStatistics_(all_spectra[i], result.protein_ids, result.peptide_ids, result.stats);
        result.stats.input_file = File::basename(in_spectra);

        // Per-file modification analysis (uses per-file calibrated tolerance).
        if (result.is_open_search)
        {
          OpenSearchModificationAnalysis mod_analyzer;
          std::string output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, per_file_cal[i].mod_match_tol,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
        }

        mfres.per_file.push_back(std::move(result));
      }
    }
    else
    {
      // ================================================================
      // Non-chunked multi-file: shared SearchContext (existing path).
      // ================================================================
      SearchContext ctx;
      bool ctx_built = false;
      // ions:by_activation: the shared index gets c and z+1 ions once a file has electron-activated
      // spectra, by rebuilding it. Only such spectra are matched against these ions (see
      // scoreSpectraAgainstIndex_), so the results of the other files depend neither on whether
      // the index holds them nor on the input order, and no file has to be read in advance.
      // A single file is the only search of the context: its index holds only the peptides that the
      // spectra (read before) can reach.
      auto prepare_context = [&](bool electron_ions, const std::function<const PeakMap*(Size)>& searched_spectra)
      {
        if (ctx_built && (ctx.electron_ions || !electron_ions)) { return; }
        StopWatch sw_idx; sw_idx.start();
        if (ctx_built)
        {
          startProgress(0, 1, "Building fragment index with c and z+1 ions...");
          ctx.fragment_index.clear();
          ctx.fragment_index.setParameters(fragmentIndexParameters_(true));
          ctx.fragment_index.build(ctx.db);
          ctx.electron_ions = true;
          endProgress();
        }
        else
        {
          // As prepareContext(fasta_db, electron_ions), but with the decoy strategy resolved above
          // instead of detecting the decoys of fasta_db a second time. If chunk_size was set but the
          // augmented DB fits in one chunk, full_db holds it already.
          if (full_db.empty())
          {
            startProgress(0, 1, "Generate decoys...");
            full_db = augment_db();
            endProgress();
          }
          ctx.db = std::move(full_db);
          ctx.decoy_string = strategy.decoy_string;
          ctx.decoy_is_prefix = strategy.is_prefix;
          ctx.have_decoys = strategy.have_decoys;
          startProgress(0, 1, "Building fragment index...");
          ctx.fragment_index.setParameters(fragmentIndexParameters_(electron_ions));
          ctx.fragment_index.build(ctx.db, searched_spectra);
          ctx.electron_ions = electron_ions;
          endProgress();
          reportDecoyParity_(ctx);
        }
        sw_idx.stop();

        // Shared report stats: index built once (and rebuilt at most once, with c and z+1 ions)
        // and reused across all files.
        mfres.shared.chunked = false;
        mfres.shared.seconds_index_build += sw_idx.getClockTime();
        mfres.shared.indexed_fragments = ctx.fragment_index.getNumFragments();
        if (ctx_built) { return; } // a rebuild: same database and peptides
        ctx_built = true;
        mfres.shared.indexed_peptides = ctx.fragment_index.getPeptides().size();
        mfres.shared.peptide_mass_mode = ctx.fragment_index.isPeptideMassMode();
        // (counted by the resolved decoy marker, as in the chunk-major path)
        mfres.shared.db_decoy_proteins = strategy.have_decoys ? countDecoyEntries_(ctx.db, strategy.decoy_string, strategy.is_prefix) : 0;
        mfres.shared.db_target_proteins = ctx.db.size() - mfres.shared.db_decoy_proteins;
      };

      mfres.per_file.reserve(in_spectra_files.size());
      // the consumer of the files' results (setFileResultConsumer()); waited for before the results change
      BackgroundConsumer file_results;

      // rescore_reuse joint/static: the searches of the first reuse_files files defer their rescoring and FDR:PSM until
      // the models are trained on them (finishDeferredRescoring_()); with 'static', the later files are rescored with
      // these models in their search.
      const bool reuse_across_files = rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static")
                                      && in_spectra_files.size() > 1 && !isOpenSearchMode_();
      if (rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static") && in_spectra_files.size() > 1 && isOpenSearchMode_())
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode " << reuse_mode_ << " needs a closed search: every file is rescored "
                        << "by Percolator trained on it ('per_file')." << std::endl;
      }
      const Size reuse_files = !reuse_across_files ? 0
                               : reuse_mode_ == "joint" ? in_spectra_files.size()
                               : std::min(reuse_static_files_, in_spectra_files.size());
      std::vector<Size> reuse_deferred;
      reuse_model_.reset();
      // the state of this search is not left behind, whichever way it ends
      const OnScopeExit reset_reuse([this]() { reuse_defer_ = false; reuse_model_.reset(); reuse_file_index_ = 0; });

      // While the first file is read: build the index, as search(file) does without c and z+1 ions
      // unless the spectra are read already or the file names an electron-based activation early on.
      // prepare_context() below adds them if the spectra need them after all.
      // With restrict_index (see above), the spectra taken for the index are searched.
      PeakMap first_spectra;
      bool first_spectra_taken = false;
      // Called for the first file, inside its progress step (the build overlaps with reading it), if its spectra are
      // not read yet.
      const auto build_index_while_first_file_is_read = [&]()
      {
        if (next_spectra.valid() && !next_spectra.ready()
            && !(ions_by_activation_ && first_file_names_electron_activation.get()))
        {
          std::function<const PeakMap*(Size)> searched_spectra;
          if (restrict_index)
          {
            const UInt64 spectra_bytes = File::fileSize(in_spectra_files[0]); // UInt64(-1) if unknown: no waiting
            searched_spectra = [&, spectra_bytes](Size peptides) -> const PeakMap*
            {
              if (!first_spectra_taken)
              {
                const bool wait = wait_for_first_spectra || (threads_started <= MAX_THREADS_WAITING_FOR_SPECTRA
                                  && spectra_bytes <= MAX_SPECTRA_BYTES_PER_PEPTIDE * peptides / threads_started);
                if (!wait && !next_spectra.ready())
                {
                  return nullptr;
                }
                first_spectra = next_spectra.get();
                first_spectra_taken = true;
              }
              return &first_spectra;
            };
          }
          prepare_context(false, searched_spectra);
        }
      };

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");

        // One progress step per input file, identified by its label, so that a progress display can tell which file
        // is searched. It ends with this iteration, whichever way it is left: the end, the continue of a failed
        // search, an exception. The steps of the search of the file nest below it (depth 1) only with one thread:
        // the background mzML reader shares the static step depth of the ProgressLogger.
        if (reuse_across_files && i == reuse_files)
        {
          finishDeferredRescoring_(mfres, reuse_deferred); // 'static': the models for this file and the later ones
          reuse_deferred.clear();
        }
        reuse_defer_ = reuse_across_files && i < reuse_files;
        reuse_file_index_ = i;

        startProgress(0, 1, "Searching " + File::basename(in_spectra) + " (" + StringUtils::toStr(i + 1) + " of "
                                + StringUtils::toStr(in_spectra_files.size()) + ")...");
        const OnScopeExit end_file_progress([this]() { endProgress(); });

        OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                        << "] Searching " << in_spectra << std::endl;

        if (i == 0) { build_index_while_first_file_is_read(); }

        PeakMap spectra = first_spectra_taken ? std::move(first_spectra)
                          : next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra, readerThreads(false), nullptr, 0, chimericSearch_());
        first_spectra_taken = false;
        read_next(i + 1);
        std::function<const PeakMap*(Size)> searched_spectra;
        if (restrict_index) { searched_spectra = [&spectra](Size) { return &spectra; }; }
        prepare_context(countElectronActivated_(spectra) > 0, searched_spectra);

        SearchResult result;
        result.is_open_search = isOpenSearchMode_();
        // No file after this one: search() releases the index right after scoring (in parallel, and
        // before the post-processing allocates) instead of the context's destructor at the end.
        ctx.release_fragment_index_after_scoring = i + 1 == in_spectra_files.size();
        reset_run_name_ = in_spectra_files[i]; // names the RESET stem table of the run (reset:stem_table a directory)
        result.exit_code = search(spectra, ctx, result.protein_ids, result.peptide_ids);

        if (result.exit_code != ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_WARN << "[ProSE] Search failed for " << in_spectra
                          << " (exit code " << static_cast<int>(result.exit_code) << "). Continuing." << std::endl;
          result.stats = last_run_stats_;
          result.percolator = last_percolator_;
          result.stats.input_file = File::basename(in_spectra);
          mfres.per_file.push_back(std::move(result));
          continue;
        }

        result.stats = last_run_stats_;
        result.percolator = last_percolator_;
        result.stats.input_file = File::basename(in_spectra);
        if (keepPreRescoringPsms_()) { result.pre_fdr_peptide_ids = std::move(pre_fdr_peptide_ids_); }
        pre_fdr_peptide_ids_ = PeptideIdentificationList();
        result.inference_only_peptide_ids = std::move(inference_only_peptide_ids_);
        inference_only_peptide_ids_ = PeptideIdentificationList();

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);

        if (result.is_open_search)
        {
          OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for " << in_spectra << std::endl;
          OpenSearchModificationAnalysis mod_analyzer;
          std::string output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, last_mod_match_tolerance_used_,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
        }
        else
        {
          OPENMS_LOG_INFO << "[ProSE] Closed search mode - per-file modification analysis skipped" << std::endl;
        }
        // After the last file, its spectra are released by a helper thread, which the caller then writes the results
        // next to (0.5 s for the 250,000 spectra of a timsTOF run at 64 threads, most of it in the kernel returning their
        // memory); the holder joins it like the hits that the rescoring dropped (see released_hits_). Otherwise here.
        bool released_later = false;
#if !defined(OPENMS_WINDOWSPLATFORM) && defined(_OPENMP) // Windows: the holder may go under the loader lock (a static object of a DLL), where joining deadlocks
        if (i + 1 == in_spectra_files.size() && spectra.size() >= 1024 && omp_get_max_threads() > 1)
        {
          try
          {
            auto release = std::make_shared<SpectraRelease>();
            release->before = released_hits_;
            // with a quarter of the threads of this thread, as the dropped hits (the release of the database and the
            // output after this use parallel regions of their own)
            release->thread = std::thread([to_release = std::move(spectra), threads = std::max(1, omp_get_max_threads() / 4)]() mutable
            {
              omp_set_num_threads(threads);
              freeInParallel(to_release);
            });
            released_hits_ = std::move(release);
            released_later = true;
          }
          catch (...) {} // the spectra (if not moved into the thread's function, and released with it) are released here
        }
#endif
        if (!released_later) { freeInParallel(spectra); }
        // after the last file, the database the index was built from
        if (i + 1 == in_spectra_files.size()) { freeInParallel(ctx.db); }

        if (reuse_defer_) { reuse_deferred.push_back(mfres.per_file.size()); }
        mfres.per_file.push_back(std::move(result));
        // (per_file holds all files' places, reserved above: the consumer's result stays where it is)
        // (with pool_rescoring_, a file that Percolator did not rescore on its own may still change: rescoreFilesTogether_())
        if (file_result_consumer_ && !reuse_defer_ && i + 1 < in_spectra_files.size()
            && mfres.per_file.back().exit_code == ExitCodes::EXECUTION_OK
            && (!pool_rescoring_ || mfres.per_file.back().percolator.applied))
        {
          file_results.start(file_result_consumer_, i, mfres.per_file.back());
        }
      }
      file_results.finish();
      finishDeferredRescoring_(mfres, reuse_deferred); // 'joint': all files (nothing without rescore_reuse)
    }

    // 'rescore': the files that Percolator could not rescore on their own are rescored with all others, with one model.
    if (pool_rescoring_) { rescoreFilesTogether_(mfres, in_spectra_files); }

    // Build the aggregate result by pooling per-file PSMs.
    // If every per-file run failed, propagate the first non-OK exit code into the
    // aggregate (rather than overloading UNEXPECTED_RESULT, which has a specific
    // meaning - PeptideIndexing failure).
    bool any_ok = false;
    for (const auto& pf : mfres.per_file)
    {
      if (pf.exit_code == ExitCodes::EXECUTION_OK) { any_ok = true; break; }
    }

    if (!any_ok)
    {
      for (const auto& pf : mfres.per_file)
      {
        if (pf.exit_code != ExitCodes::EXECUTION_OK)
        {
          mfres.aggregate.exit_code = pf.exit_code;
          break;
        }
      }
      return mfres;
    }

    // Single-file fast path: the aggregate would just duplicate the only per-file
    // result and re-run modification analysis on the same PSMs. Skip the pooling
    // and analysis entirely. The aggregate is left with only @c is_open_search
    // and @c exit_code populated; callers should use @c per_file[0] for the
    // actual identifications. This is documented on @c MultiFileSearchResult.
    if (mfres.per_file.size() == 1 && mfres.per_file[0].exit_code == ExitCodes::EXECUTION_OK)
    {
      return mfres;
    }

    // Pooling copies every per-file PSM a second time. Skip it entirely when the caller
    // does not want the pooled set AND nothing here needs it: only the open-search
    // aggregate modification analysis below consumes it, so in closed search a
    // build_pooled_aggregate == false caller gets the same almost-empty aggregate as the
    // single-file fast path above. In open search the pooled set is built, analyzed and
    // then released again (see the release block after the analysis).
    const bool need_pooled = build_pooled_aggregate || mfres.aggregate.is_open_search;

    if (need_pooled)
    {
      // Merge per-file identifications into a single aggregate using
      // IDMergerAlgorithm — the canonical OpenMS pattern for cross-file protein
      // inference (mirrors ConsensusMapMergerAlgorithm::mergeAllIDRuns in
      // ProteomicsLFQ). This deduplicates ProteinHits by accession (union) and
      // remaps all PeptideIdentification identifiers to the merged run. No
      // PeptideIndexing re-run needed: per-file searches already linked
      // peptides to proteins.
      // Copying is mandatory here: per_file is returned to the caller and must
      // stay intact, so nothing may be moved out of it.
      poolSearchedFiles_(mfres.per_file, /*move=*/false, in_spectra_files,
                         mfres.aggregate.protein_ids, mfres.aggregate.peptide_ids);
    }

    // Note: protein inference + FDR on the aggregate is left to the caller
    // (e.g., the TOPP tool's -out_merged option) so that per-file outputs
    // retain run-level information without being overwritten by aggregate
    // protein lists. The aggregate here contains the merged (unfiltered)
    // proteins + pooled PSMs for downstream use.

    // Aggregate modification analysis on the pooled PSM set.
    if (mfres.aggregate.is_open_search && !mfres.aggregate.peptide_ids.empty())
    {
      OPENMS_LOG_INFO << "[ProSE] Running aggregate modification analysis on "
                      << mfres.aggregate.peptide_ids.size() << " pooled PSM(s) from "
                      << in_spectra_files.size() << " input file(s)." << std::endl;

      OpenSearchModificationAnalysis mod_analyzer;
      std::string agg_output_file;
      if (!aggregate_base_name.empty())
      {
        agg_output_file = aggregate_base_name + "_ModificationAnalysis.idXML";
      }

      mfres.aggregate.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
        mfres.aggregate.peptide_ids,
        computeModMatchTolerance_(),
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        agg_output_file
      );
    }

    // The pooled PSMs were only needed as input to the analysis above; release the second
    // copy again so a caller that asked for no pooled aggregate does not pay for it. The
    // analysis result stays. Leaves the aggregate exactly as in the single-file fast path
    // (only is_open_search / exit_code populated), as documented on MultiFileSearchResult.
    if (!build_pooled_aggregate)
    {
      mfres.aggregate.peptide_ids.clear();
      mfres.aggregate.peptide_ids.shrink_to_fit();
      mfres.aggregate.protein_ids.clear();
      mfres.aggregate.protein_ids.shrink_to_fit();
    }

    return mfres;
  }

  // =====================================================================
  // Multi-file searchWithModificationAnalysis (FASTA file path).
  //
  // Loads the FASTA from disk once, delegates to the in-memory multi-file
  // overload, and patches the database file path on each per-file (and the
  // aggregate) ProteinIdentification.
  // =====================================================================
  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      const std::vector<std::string>& in_spectra_files,
      const std::string& in_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    // load FASTA once (while the first spectrum file is read in the background)
    MultiFileSearchResult mfres = searchWithModificationAnalysis_(
      in_spectra_files, {}, &in_db, output_base_names, aggregate_base_name, build_pooled_aggregate);

    mfres.shared.database_file = in_db;

    // Patch database file path into each per-file and aggregate ProteinIdentification.
    for (auto& pf : mfres.per_file)
    {
      if (pf.exit_code == ExitCodes::EXECUTION_OK && !pf.protein_ids.empty())
      {
        pf.protein_ids[0].getSearchParameters().db = in_db;
      }
    }
    if (mfres.aggregate.exit_code == ExitCodes::EXECUTION_OK && !mfres.aggregate.protein_ids.empty())
    {
      mfres.aggregate.protein_ids[0].getSearchParameters().db = in_db;
    }

    return mfres;
  }

  // =====================================================================
  // Pooling of the searched files of a multi-file search (IDMergerAlgorithm), shared
  // by the search's pooled aggregate and mergeFileResults().
  // =====================================================================
  // static
  void ProSEAlgorithm::poolSearchedFiles_(std::vector<SearchResult>& per_file, bool move_psms,
                                          const StringList& in_spectra_files,
                                          std::vector<ProteinIdentification>& protein_ids,
                                          PeptideIdentificationList& peptide_ids)
  {
    // IDMergerAlgorithm merges only runs of one search engine (and version). With 'rescore', a run that Percolator did
    // not rescore keeps ProSE as its search engine, while the rescored ones report Percolator. The copies of the runs
    // that are pooled all get the engine of the files with PSMs: Percolator's if it rescored all of them (a file without
    // PSMs is not rescored), else ProSE's (the pooled evaluation then scores all files by ProSE's score, see
    // mergeFileResults()); the runs of the files stay as they are.
    // If no file has a PSM left (e.g. FDR:PSM removed all), the engine of the searched files by the same rule:
    // IDMergerAlgorithm requires one engine of runs without PSMs as well.
    const SearchResult* engine_of = nullptr; // the first file that Percolator did not rescore, else the first
    for (const bool require_psms : {true, false})
    {
      for (const SearchResult& pf : per_file)
      {
        if (pf.exit_code != ExitCodes::EXECUTION_OK || pf.protein_ids.empty()) { continue; }
        if (require_psms && !hasPSMs_(pf.peptide_ids)) { continue; }
        if (engine_of == nullptr || (engine_of->percolator.applied && !pf.percolator.applied)) { engine_of = &pf; }
      }
      if (engine_of != nullptr) { break; }
    }
    std::string engine, engine_version;
    if (engine_of != nullptr)
    {
      engine = engine_of->protein_ids[0].getSearchEngine();
      engine_version = engine_of->protein_ids[0].getSearchEngineVersion();
    }

    IDMergerAlgorithm merger;
    std::vector<Size> pooled_files; // input positions of the pooled files, in insertion order
    for (Size i = 0; i < per_file.size(); ++i)
    {
      SearchResult& pf = per_file[i];
      if (pf.exit_code != ExitCodes::EXECUTION_OK) { continue; }
      if (!pf.protein_ids.empty()) { pooled_files.push_back(i); } // insertRuns() skips a file without a run
      // Moving: the merger takes (and clears) the PSMs and a copy of the runs, which stay with the files -- their search
      // parameters, meta values (e.g. PeptDeepRescoring's NCE) and run paths are still the caller's.
      std::vector<ProteinIdentification> runs = pf.protein_ids;
      if (!engine.empty())
      {
        for (ProteinIdentification& run : runs)
        {
          run.setSearchEngine(engine);
          run.setSearchEngineVersion(engine_version);
        }
      }
      if (move_psms) { merger.insertRuns(std::move(runs), std::move(pf.peptide_ids)); }
      else { merger.insertRuns(std::move(runs), PeptideIdentificationList(pf.peptide_ids)); }
    }
    ProteinIdentification merged_proteins;
    merger.returnResultsAndClear(merged_proteins, peptide_ids);
    protein_ids.clear();
    protein_ids.emplace_back(std::move(merged_proteins));
    protein_ids[0].setPrimaryMSRunPath(in_spectra_files);

    // IDMergerAlgorithm numbers the primary MS run paths of the runs ('id_merge_index') in the order it meets them --
    // one per pooled file, as the search gives each run the path of its input file -- so it counts only the pooled
    // files, while the run lists all input files: after a failed file, index k would name the wrong file. Map k to the
    // input position of the k-th pooled file (nothing to do if no file failed).
    if (!pooled_files.empty() && pooled_files.back() + 1 != pooled_files.size())
    {
      for (PeptideIdentification& id : peptide_ids)
      {
        if (!id.metaValueExists(Constants::UserParam::ID_MERGE_INDEX)) { continue; }
        const Size k = static_cast<Size>(static_cast<int>(id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)));
        if (k < pooled_files.size()) { id.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, pooled_files[k]); }
      }
    }
  }

  // =====================================================================
  // 'rescore' of a multi-file search: the files that Percolator could not rescore on their own (e.g. fewer than 100
  // PSMs) get the scores of a model trained on the files that it could rescore, or, if it could rescore none, all
  // files are rescored together.
  // =====================================================================
  // static
  std::vector<Size> ProSEAlgorithm::stratifiedQuotas_(const std::vector<Size>& sizes, Size cap)
  {
    std::vector<Size> quotas = sizes;
    Size total = std::accumulate(sizes.begin(), sizes.end(), Size(0));
    if (cap == 0 || total <= cap) { return quotas; }
    // water filling: equal shares, the share a small stratum leaves goes to the others
    std::vector<Size> order(sizes.size());
    std::iota(order.begin(), order.end(), Size(0));
    std::sort(order.begin(), order.end(), [&](Size a, Size b) { return sizes[a] < sizes[b]; });
    Size left = cap;
    for (Size k = 0; k < order.size(); ++k)
    {
      const Size share = left / (order.size() - k);
      quotas[order[k]] = std::min(sizes[order[k]], share);
      left -= quotas[order[k]];
    }
    return quotas;
  }

  void ProSEAlgorithm::rescoreFilesTogether_(MultiFileSearchResult& mfres, const std::vector<std::string>& in_spectra_files) const
  {
    // Each file ends as a search without pooling would leave it: its proteins those of its PSMs (and of the PSMs only
    // for the protein inference), no copy of its PSMs before the rescoring unless asked for (setKeepPreFdrPsms()).
    const auto finish = [this](SearchResult& pf)
    {
      removeUnreportedProteins(pf.protein_ids, pf.peptide_ids, pf.inference_only_peptide_ids);
      if (!keep_pre_fdr_psms_) { pf.pre_fdr_peptide_ids = PeptideIdentificationList(); }
      else if (rescore_top_hits_ > 0) { keepTopHits_(pf.pre_fdr_peptide_ids, report_top_hits_); } // (see applyPsmFDR_())
    };
    std::vector<Size> searched, own, alone; // files searched; rescored on their own; with PSMs, but not rescored
    StringList alone_names;
    for (Size i = 0; i < mfres.per_file.size(); ++i)
    {
      const SearchResult& pf = mfres.per_file[i];
      if (pf.exit_code != ExitCodes::EXECUTION_OK || pf.protein_ids.empty()) { continue; }
      searched.push_back(i);
      if (pf.percolator.applied) { own.push_back(i); }
      else if (hasPSMs_(pf.pre_fdr_peptide_ids))
      {
        alone.push_back(i);
        alone_names.push_back(File::basename(in_spectra_files[i]) + " (" + pf.percolator.skipped_reason + ")");
      }
    }
    const auto finish_all = [&]() { for (Size i : searched) { finish(mfres.per_file[i]); } };
    if (alone.empty())
    {
      finish_all();
      return;
    }
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    // Percolator as for the files rescored on their own (rescoring:..., the training cap percolator:max_train_psms)
    const Param pp = rescoringPercolatorParameters(Percolator().getDefaults(), threads, percolator_training_, Param());

    // A file back from the rescoring (its PSMs before the rescoring, now with Percolator's score, q-value and PEP):
    // q-values per file (FDR per run, as for a file rescored on its own), by target-decoy competition of the file's
    // PSMs on the SVM score (MS:1001492), (D + 1) / T as Percolator's own (tdcQValues_(), also 1/T for a file without
    // a decoy PSM); FDR:PSM filters on them.
    const std::string svm_score = "MS:1001492";
    const auto adopt = [&](Size i, PeptideIdentificationList&& ids, const std::string& warning)
    {
      SearchResult& pf = mfres.per_file[i];
      for (PeptideIdentification& id : ids)
      {
        for (PeptideHit& hit : id.getHits()) { hit.setScore(static_cast<double>(hit.getMetaValue(svm_score))); }
        id.setScoreType(svm_score);
        id.setHigherScoreBetter(true);
      }
      if (!ids.empty()) { tdcQValues_(ids, true); }
      for (PeptideIdentification& id : ids)
      {
        for (PeptideHit& hit : id.getHits())
        {
          hit.removeMetaValue(svm_score + "_score"); // (FalseDiscoveryRate's copy of the SVM score)
          hit.setMetaValue("MS:1001491", hit.getScore());
          hit.setMetaValue("percolator_q_value", hit.getScore());
        }
      }
      // rescoring:stratum_q_values: as for a file rescored on its own (applyRescoredPsmFDR_()), before FDR:PSM
      if (stratum_q_values_ && !ids.empty())
      {
        annotateStratumQValues(ids);
        for (ProteinIdentification& run : pf.protein_ids) run.getSearchParameters().setMetaValue("rescoring:stratum_q_values", "true");
      }
      pf.inference_only_peptide_ids = PeptideIdentificationList();
      bool filtered = false;
      if (fdr_psm_ > 0.0)
      {
        if (fdr_protein_ > 0.0) { collectInferenceOnlyPsms_(ids, fdr_psm_, pf.inference_only_peptide_ids); }
        IDFilter::filterHitsByScore(ids, fdr_psm_);
        filtered = true;
      }
      pf.peptide_ids = std::move(ids);
      markRescoredByPercolator_(pf.protein_ids, pp, percolator_training_);
      pf.percolator.applied = true;
      pf.percolator.skipped_reason.clear();
      pf.percolator.warning = warning;
      updateFinalStats(pf.stats, pf.peptide_ids, enzyme_, filtered);
      // (an open search keeps the modification analysis of the PSMs of the file's own scoring)
    };

    // Percolator gets the first rescoring:top_hits hits of a spectrum, as for a file rescored on its own (the copies
    // before the rescoring may hold more, for a .pin with report:top_hits above it; see applyPsmFDR_())
    const Size rescored_hits = rescore_top_hits_ > 0 ? rescoreTopHits_() : std::numeric_limits<Size>::max();
    const auto hits_rescored = [rescored_hits](const PeptideIdentification& id) { return std::min(id.getHits().size(), rescored_hits); };
    const auto cut_to_rescored = [&](PeptideIdentificationList& ids) { if (rescore_top_hits_ > 0) keepTopHits_(ids, rescored_hits); };

    // The precursor charges of the charge features: of all PSMs the model sees, for one feature set
    int min_charge = std::numeric_limits<int>::max(), max_charge = std::numeric_limits<int>::min();
    for (Size i : searched)
    {
      for (const PeptideIdentification& id : mfres.per_file[i].pre_fdr_peptide_ids)
      {
        for (Size h = 0; h < hits_rescored(id); ++h)
        {
          min_charge = std::min(min_charge, id.getHits()[h].getCharge());
          max_charge = std::max(max_charge, id.getHits()[h].getCharge());
        }
      }
    }
    const ProteinIdentification& first_run = mfres.per_file[searched[0]].protein_ids[0];
    const StringList feature_set = rescoringFeatureSet(first_run.getSearchParameters(), min_charge, max_charge);
    const std::string enzyme = rescoringEnzyme(first_run.getSearchParameters().digestion_enzyme.getName(), percolator_training_);

    // rescore_method 'reset' has no transfer of a model (Percolator::trainPSMs()/scorePSMs() are plain Percolator): the
    // files too small for a rescoring of their own keep their search engine scores, as before.
    if (!own.empty() && percolator_extra_.exists("algorithm"))
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator could not rescore " << alone.size() << " of " << searched.size()
                      << " input files on their own: " << ListUtils::concatenate(alone_names, ", ")
                      << ". rescore_method 'reset' cannot score them with a model of the other files; they keep their "
                      << "search engine scores." << std::endl;
      finish_all();
      return;
    }
    if (!own.empty())
    {
      // Transfer: the files rescored on their own keep their results. A model trained on their PSMs (before the
      // rescoring; each file an equal share of at most percolator:max_train_psms, so that none dominates) scores the
      // others, on the scale of its training PSMs: a PSM gets the score it would get among them, and the PEP of the
      // training PSMs at that score (Percolator::scorePSMs()), not one estimated on a file too small for that.
      // 'id_merge_index' keeps the spectra of the training files apart (spectrum file numbers).
      // The cap counts PSMs (Percolator's rows: every hit of a spectrum), the sample takes whole spectra.
      std::vector<std::vector<Size>> rows(own.size()); // of each file, those with PSMs
      std::vector<Size> sizes; // PSMs of each file
      for (Size k = 0; k < own.size(); ++k)
      {
        const PeptideIdentificationList& ids = mfres.per_file[own[k]].pre_fdr_peptide_ids;
        Size psms = 0;
        for (Size r = 0; r < ids.size(); ++r)
        {
          if (!ids[r].getHits().empty()) { rows[k].push_back(r); psms += hits_rescored(ids[r]); }
        }
        sizes.push_back(psms);
      }
      const std::vector<Size> quotas = stratifiedQuotas_(sizes, max_train_psms_);
      PeptideIdentificationList training;
      for (Size k = 0; k < own.size(); ++k)
      {
        const PeptideIdentificationList& ids = mfres.per_file[own[k]].pre_fdr_peptide_ids;
        if (quotas[k] < sizes[k])
        {
          // reproducible, and the same on every platform (std::shuffle is not); whole spectra up to the file's share
          Math::RandomShuffler(static_cast<int>(own[k] + 1)).portable_random_shuffle(rows[k].begin(), rows[k].end());
          std::vector<Size> kept;
          Size psms = 0;
          for (Size r : rows[k])
          {
            const Size n = hits_rescored(ids[r]);
            if (psms + n > quotas[k]) { continue; }
            kept.push_back(r);
            psms += n;
          }
          rows[k] = std::move(kept);
          std::sort(rows[k].begin(), rows[k].end());
        }
        for (Size r : rows[k])
        {
          PeptideIdentification& copy = training.emplace_back(ids[r]);
          copy.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<int>(own[k]));
        }
      }
      cut_to_rescored(training);
      const Size available = std::accumulate(sizes.begin(), sizes.end(), Size(0));
      Size training_psms = 0;
      for (const PeptideIdentification& id : training) { training_psms += id.getHits().size(); }
      OPENMS_LOG_INFO << "[ProSE] Percolator could not rescore " << alone.size() << " of " << searched.size()
                      << " input files on their own: " << ListUtils::concatenate(alone_names, ", ")
                      << ". Scoring them with a model trained on " << training_psms << " PSMs of the " << own.size()
                      << " files it rescored on their own";
      if (training_psms < available)
      {
        OPENMS_LOG_INFO << " (sampled from " << available << ", an equal share of each file, Search:percolator:max_train_psms "
                        << max_train_psms_ << ")";
      }
      OPENMS_LOG_INFO << "." << std::endl;
      try
      {
        Percolator perc;
        perc.setParameters(pp);
        const PercolatorModel model = perc.trainPSMs(training, feature_set, enzyme, min_charge, max_charge);
        const std::string transferred = "its PSMs are scored with a model trained on the PSMs of the "
          + StringUtils::toStr(own.size()) + " input(s) that Percolator rescored on their own, on the scale of that "
          "model's training PSMs (their PEP at the same score); the q-values are those of target-decoy competition of "
          "this input's PSMs on that score." + (perc.usedDefaultDirection() ? " Percolator's sanity check rejected "
          "the learned model, so the PSMs are scored with the best single feature." : "");
        // The PSMs of all these files in one call: it scores every PSM on its own (no target-decoy competition among
        // them), and prepares and scores the training PSMs (the reference) once.
        PeptideIdentificationList ids;
        std::vector<Size> file_of; // per element of ids
        for (Size i : alone)
        {
          PeptideIdentificationList& pre_fdr = mfres.per_file[i].pre_fdr_peptide_ids;
          for (PeptideIdentification& id : pre_fdr)
          {
            if (keep_pre_fdr_psms_) { ids.push_back(id); }
            else { ids.push_back(std::move(id)); }
            file_of.push_back(i);
          }
          if (!keep_pre_fdr_psms_) { pre_fdr = PeptideIdentificationList(); }
        }
        // the best hit of each spectrum by the SVM score (the q-values of the training PSMs at a score have plateaus,
        // which would leave ties to the search engine's order); target-decoy competition follows in adopt()
        cut_to_rescored(ids);
        perc.scorePSMs(ids, model, training, feature_set, enzyme, min_charge, max_charge, "svm");
        IDFilter::keepNBestHits(ids, 1);
        std::vector<PeptideIdentificationList> parts(mfres.per_file.size());
        for (Size k = 0; k < ids.size(); ++k) { parts[file_of[k]].push_back(std::move(ids[k])); }
        ids = PeptideIdentificationList();
        for (Size i : alone)
        {
          adopt(i, std::move(parts[i]), "Percolator could not rescore this input on its own ("
                                        + mfres.per_file[i].percolator.skipped_reason + "): " + transferred);
        }
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "[ProSE] Percolator cannot score the files too small for a model of their own with a model "
                        << "of the other files either (" << e.what() << "); they keep their search engine scores."
                        << std::endl;
      }
      finish_all();
      return;
    }

    // No file could be rescored on its own: the PSMs of all files together, with one model, cross-validated as for one
    // file (Percolator samples at most percolator:max_train_psms PSMs to train on; every file has too few PSMs for a
    // model of its own, so none can dominate). 'id_merge_index' keeps the files apart (spectrum file numbers, for the
    // folds and the target-decoy competition of each spectrum).
    PeptideIdentificationList pooled;
    std::vector<Size> file_of; // per element of pooled
    std::vector<ProteinIdentification> runs = mfres.per_file[searched[0]].protein_ids; // the search parameters
    for (Size i : searched)
    {
      for (const PeptideIdentification& id : mfres.per_file[i].pre_fdr_peptide_ids)
      {
        PeptideIdentification& copy = pooled.emplace_back(id);
        copy.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<int>(i));
        file_of.push_back(i);
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Percolator could not rescore any of the " << searched.size()
                    << " input files on its own: " << ListUtils::concatenate(alone_names, ", ")
                    << ". Rescoring the " << pooled.size() << " PSMs of all files together, with one model." << std::endl;
    // (rescore_method 'reset': without a stem table, which names one run)
    Param extra = percolator_extra_;
    if (extra.exists("reset:stem_table")) { extra.setValue("reset:stem_table", ""); }
    cut_to_rescored(pooled);
    const PercolatorRescoring outcome = rescoreWithPercolator_(runs, pooled, threads, nullptr, percolator_training_, extra);
    if (!outcome.applied)
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator cannot rescore the PSMs of all input files together either ("
                      << outcome.skipped_reason << "); they keep their search engine scores." << std::endl;
      finish_all();
      return;
    }
    std::vector<PeptideIdentificationList> parts(mfres.per_file.size());
    for (Size k = 0; k < pooled.size(); ++k)
    {
      pooled[k].removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
      parts[file_of[k]].push_back(std::move(pooled[k]));
    }
    pooled = PeptideIdentificationList();
    const std::string warning = "Percolator could not rescore any input on its own (too few PSMs or decoys); the PSMs of "
      "all " + StringUtils::toStr(searched.size()) + " inputs are rescored together, with one model; the q-values are "
      "those of target-decoy competition of this input's PSMs on that model's score." + (outcome.warning.empty() ? std::string() : " " + outcome.warning);
    for (Size i : searched) { adopt(i, std::move(parts[i]), warning); }
    finish_all();
  }

  // =====================================================================
  // The pooled evaluation of a multi-file search (the ProSE tool's -out_merged):
  // pooling, cross-file protein inference and the optional picked-protein FDR.
  // =====================================================================
  ProSEAlgorithm::MergedResult
  ProSEAlgorithm::mergeFileResults(MultiFileSearchResult& result, const StringList& in_spectra_files,
                                   bool remove_decoys) const
  {
    // The pooled run lists the input files, and 'id_merge_index' indexes them: one per file of the search.
    if (in_spectra_files.size() != result.per_file.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "mergeFileResults(): " + StringUtils::toStr(in_spectra_files.size()) + " input file(s) given for the result of a "
        "search of " + StringUtils::toStr(result.per_file.size()) + " file(s); pass the input files of the search, in its "
        "order.");
    }

    // Percolator ('rescore') rescores each file on its own; a file it could not rescore (e.g. one with fewer than 100
    // PSMs) keeps ProSE's scores. The PSMs of each file passed FDR:PSM on the q-values of their file, Percolator's or
    // target-decoy q-values of ProSE's score, which are valid per file and stay. Only the protein inference compares
    // PSMs across files, by the best PEP of a protein's PSMs: the PSMs that Percolator did not rescore have none. A
    // mixed set is therefore pooled on ProSE's own score, which every PSM has: the inference scores proteins by the
    // best 'ln(hyperscore)' (as without rescoring, see inferProteins()), and with FDR:PSM 0, where the main score is
    // not a q-value, every pooled PSM gets 'ln(hyperscore)' as main score (Percolator's scores stay meta values).
    Size n_searched = 0, n_rescored = 0;
    std::vector<std::string> not_rescored;
    for (Size i = 0; i < result.per_file.size(); ++i)
    {
      const SearchResult& pf = result.per_file[i];
      if (pf.exit_code != ExitCodes::EXECUTION_OK || !hasPSMs_(pf.peptide_ids)) { continue; } // (no PSM to compare)
      ++n_searched;
      if (pf.percolator.applied) { ++n_rescored; }
      else { not_rescored.push_back(File::basename(in_spectra_files[i]) + ": " + pf.percolator.skipped_reason); }
    }
    const bool rescoring_mixed = n_rescored > 0 && n_rescored < n_searched;
    if (rescoring_mixed)
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator rescored " << n_rescored << " of " << n_searched << " input files; not "
                      << ListUtils::concatenate(not_rescored, " ") << " The pooled evaluation scores the proteins of "
                      << "all files by ProSE's score ln(hyperscore) instead of Percolator's PEP, which only the rescored "
                      << "PSMs have; each PSM keeps the q-value of its file." << std::endl;
    }

    MergedResult merged;
    if (n_searched == 0) { return merged; } // no file was searched: nothing to pool

    // Use the effective decoy marker/position the search resolved (which may be an auto-detected external marker,
    // prefix or suffix), not the raw decoy_prefix parameter.
    const std::string& decoy_string = result.decoy_string;
    const bool decoy_is_prefix = result.decoy_is_prefix;

    // The multi-file search already produced exactly this merge: it pools the per-file results with
    // poolSearchedFiles_(), as below -- so re-merging here would only build a third copy of every PSM. Consume the
    // aggregate by move instead. Without an aggregate (build_pooled_aggregate false), pool the runs of the files here.
    if (!result.aggregate.protein_ids.empty())
    {
      merged.protein_ids.emplace_back(std::move(result.aggregate.protein_ids[0]));
      merged.peptide_ids = std::move(result.aggregate.peptide_ids);
      result.aggregate.protein_ids.clear();
      result.aggregate.peptide_ids.clear();
    }
    else
    {
      // The PSMs are moved out of per_file, as documented; its runs are copied and stay, as do the statistics, exit
      // codes, Percolator outcomes and modification analyses.
      poolSearchedFiles_(result.per_file, /*move_psms=*/true, in_spectra_files, merged.protein_ids, merged.peptide_ids);
    }

    std::vector<ProteinIdentification>& merged_protein_ids = merged.protein_ids;
    PeptideIdentificationList& merged_peptides = merged.peptide_ids;
    merged_protein_ids[0].setPrimaryMSRunPath(in_spectra_files);
    // IDMergerAlgorithm carries over search engine/params but not a run date; set one so the
    // merged idXML is schema-valid (matches the single-file path).
    merged_protein_ids[0].setDateTime(DateTime::now());

    // Without FDR:PSM the main scores of a mixed set differ: Percolator's q-values and ProSE's scores. One score for all.
    if (rescoring_mixed && fdr_psm_ <= 0.0) { useSearchEngineScore_(merged_peptides); }

    // Cross-file protein inference (inferProteins(): the best PEP, or search engine score, per peptide per protein) and
    // the optional picked-protein FDR. Picked-protein FDR needs identified decoy proteins: gate on DB-level decoys here,
    // applyCompleteSetProteinFDR() checks the result of the inference (it could be target-only even with a decoy
    // database if no decoy survived in any file). It infers the proteins on the pooled PSMs together with the PSMs of
    // the files between FDR:PSM and FDR:protein_inference_PSM, and finalizes the result: with remove_decoys it strips
    // the decoys, and in any case the proteins without a pooled PSM and the references to removed proteins (peptide
    // evidences, indistinguishable proteins and protein groups), so that the result can still be stored as idXML.
    // Without the protein FDR the merged output is a pooled intermediate that retains the target+decoy evidence for
    // a downstream FDR.
    if (fdr_protein_ > 0.0 && result.have_decoys && !decoy_string.empty())
    {
      PeptideIdentificationList inference_only;
      for (SearchResult& pf : result.per_file)
      {
        if (pf.exit_code != ExitCodes::EXECUTION_OK) { continue; }
        inference_only.insert(inference_only.end(), std::make_move_iterator(pf.inference_only_peptide_ids.begin()),
                              std::make_move_iterator(pf.inference_only_peptide_ids.end()));
        pf.inference_only_peptide_ids = PeptideIdentificationList();
      }
      merged.protein_fdr_applied = applyCompleteSetProteinFDR(merged_protein_ids, merged_peptides, decoy_string,
                                                              decoy_is_prefix, fdr_protein_, remove_decoys, inference_only);
    }
    else
    {
      IDFilter::keepNBestHits(merged_peptides, 1); // as the protein FDR keeps them
      inferProteins(merged_protein_ids, merged_peptides);
    }

    return merged;
  }

  // =====================================================================
  // File-based single-file searchWithModificationAnalysis: thin wrapper around
  // the multi-file overload using a single-element list.
  // =====================================================================
  ProSEAlgorithm::SearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(const std::string& in_spectra,
                                                                  const std::string& in_db,
                                                                  const std::string& output_base_name) const
  {
    std::vector<std::string> in_files{in_spectra};
    std::vector<std::string> base_names;
    if (!output_base_name.empty()) { base_names.push_back(output_base_name); }

    MultiFileSearchResult mfres = searchWithModificationAnalysis(in_files, in_db, base_names, "");

    if (mfres.per_file.empty())
    {
      SearchResult empty_result;
      empty_result.exit_code = ExitCodes::INPUT_FILE_EMPTY;
      return empty_result;
    }

    // A single input is a complete experiment: protein FDR as in search(file).
    SearchResult& result = mfres.per_file[0];
    if (fdr_protein_ > 0.0 && mfres.have_decoys && result.exit_code == ExitCodes::EXECUTION_OK
        && !result.protein_ids.empty() && !result.peptide_ids.empty())
    {
      applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids, mfres.decoy_string, mfres.decoy_is_prefix, fdr_protein_,
                                 true, result.inference_only_peptide_ids);
      updateFinalStats(result.stats, result.peptide_ids, enzyme_, result.stats.fdr_applied);
    }
    result.inference_only_peptide_ids = PeptideIdentificationList(); // used by the protein FDR only
    return std::move(result);
  }

  // =====================================================================
  // Helper: decoy-marker membership test (prefix/suffix), pure std::string.
  // =====================================================================
  bool ProSEAlgorithm::accessionHasDecoyMarker_(const std::string& accession,
                                                const std::string& marker, bool is_prefix)
  {
    if (marker.empty() || accession.size() < marker.size()) { return false; }
    return is_prefix ? (accession.compare(0, marker.size(), marker) == 0)
                     : (accession.compare(accession.size() - marker.size(), marker.size(), marker) == 0);
  }

  bool ProSEAlgorithm::useSearchEngineScore_(PeptideIdentificationList& ids)
  {
    // The name of the score of each PSM: none for those that have it as main score already
    const std::string main_score = SEARCH_ENGINE_SCORE, meta_score = main_score + "_score";
    std::vector<const std::string*> source(ids.size(), nullptr);
    for (Size i = 0; i < ids.size(); ++i)
    {
      const PeptideIdentification& id = ids[i];
      if (id.getHits().empty() || id.getScoreType() == main_score) { continue; }
      const PeptideHit& hit = id.getHits()[0];
      if (hit.metaValueExists(main_score)) { source[i] = &main_score; }
      else if (hit.metaValueExists(meta_score)) { source[i] = &meta_score; }
      else { return false; }
    }
    // one switcher per meta value name
    const auto switcher = [&main_score](const std::string& name)
    {
      IDScoreSwitcherAlgorithm s;
      Param p = s.getParameters();
      p.setValue("new_score", name);
      p.setValue("new_score_type", main_score);
      p.setValue("new_score_orientation", "higher_better");
      s.setParameters(p);
      return s;
    };
    IDScoreSwitcherAlgorithm from_main = switcher(main_score), from_meta = switcher(meta_score);
    Size counter = 0;
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (source[i] != nullptr) { (source[i] == &main_score ? from_main : from_meta).switchScores(ids[i], counter); }
    }
    return true;
  }

  Size ProSEAlgorithm::countDecoyEntries_(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& marker, bool is_prefix)
  {
    return static_cast<Size>(std::count_if(db.begin(), db.end(), [&](const FASTAFile::FASTAEntry& e)
                                           { return accessionHasDecoyMarker_(e.identifier, marker, is_prefix); }));
  }

  // =====================================================================
  // Helper: capture pre-FDR statistics (target/decoy counts + HyperScore
  // distribution). Must run BEFORE FalseDiscoveryRate, which overwrites the
  // hit score with the q-value and may drop decoy hits.
  // =====================================================================
  void ProSEAlgorithm::capturePreFdrStats_(const PeptideIdentificationList& peptide_ids,
                                           RunStatistics& stats)
  {
    stats.target_psms = 0;
    stats.decoy_psms = 0;
    // The top hits of the spectra are read in parallel (the target/decoy value by name index: looking the name up takes
    // a lock on the registry, which the threads would contend for; by name if it is not registered), then their scores
    // collected and counted in spectrum order, as before (also where an error stops it).
    enum TopHit : char { NO_HIT, OTHER, TARGET, DECOY };
    const UInt td_index = MetaInfoInterface::metaRegistry().getIndex(Constants::UserParam::TARGET_DECOY);
    const SignedSize num_ids = static_cast<SignedSize>(peptide_ids.size());
    std::vector<TopHit> top_hits(peptide_ids.size(), NO_HIT);
    std::vector<double> top_scores(peptide_ids.size());
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(static, 1024)
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      try
      {
        const PeptideIdentification& pid = peptide_ids[i];
        if (pid.getHits().empty()) { continue; }
        const PeptideHit& top = pid.getHits().front();
        top_hits[i] = OTHER;
        top_scores[i] = top.getScore();

        const bool by_name = td_index == UInt(-1);
        if (by_name ? !top.metaValueExists(Constants::UserParam::TARGET_DECOY) : !top.metaValueExists(td_index)) { continue; }
        const std::string td = (by_name ? top.getMetaValue(Constants::UserParam::TARGET_DECOY) : top.getMetaValue(td_index)).toString();
        // "target+decoy" counts as target (OpenMS FDR semantics).
        if (td == "decoy") { top_hits[i] = DECOY; }
        else if (!td.empty()) { top_hits[i] = TARGET; }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    std::vector<double> hyperscores;
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      if (top_hits[i] == NO_HIT) { continue; }
      hyperscores.push_back(top_scores[i]);
      if (errors[i]) { std::rethrow_exception(errors[i]); }
      if (top_hits[i] == DECOY) { ++stats.decoy_psms; }
      else if (top_hits[i] == TARGET) { ++stats.target_psms; }
    }

    if (!hyperscores.empty())
    {
      auto [mn, mx] = std::minmax_element(hyperscores.begin(), hyperscores.end());
      stats.hyperscore_min = *mn;
      stats.hyperscore_max = *mx;
      stats.hyperscore_median = Math::median(hyperscores.begin(), hyperscores.end());
      stats.score_stats_valid = true;
    }
  }

  namespace
  {
    // What collectRunStatistics_() and updateFinalStats() count over the top hits
    struct TopHitSummary
    {
      Size matched = 0, n_target = 0, n_decoy = 0;
      std::map<Int, Size> charges;
      std::map<Size, Size> missed_cleavages;
      Size unique_peptides = 0, unique_proteins = 0; // distinct sequence strings and protein accessions
      std::vector<double> precursor_errors, fragment_errors; // top hits in input order (with_errors only)
      Size fragment_errors_psms = 0; // the PSMs with any fragment error (with_errors only)
    };

    // One pass over the top hits, in parallel: each thread takes a contiguous range of the identifications, and the
    // ranges are merged in order, so the error lists come out in input order as from a single pass. The regular
    // expression of the missed-cleavage count and the sequence strings made this a serial step of 20 - 40 ms per file.
    // The distinct peptides and proteins are counted in parallel, too: each range notes its sequence strings (in one
    // buffer) and accessions in buckets by hash, and the buckets are counted one per thread. Merging sets of strings
    // per range and freeing them afterwards took 0.6 s on one thread for the 250,000 top hits of a timsTOF run.
    //
    // The mass errors of the tolerance estimate (with_errors), of the top hits of the target PSMs (a decoy is a wrong
    // match: its errors are noise), signed (observed minus theoretical, ppm): of the precursor, corrected by the isotope
    // error of the hit (as the search matched it; @p isotope_observed_minus_theoretical: the convention of the reported
    // 'isotope_error'), without the hits with a delta mass of an open search; and of the annotated b and y fragment peaks.
    TopHitSummary summarizeTopHits(const PeptideIdentificationList& peptide_ids, const std::string& enzyme, bool with_errors,
                                   bool isotope_observed_minus_theoretical = true)
    {
      const DigestionEnzymeProtein* protease = ProteaseDB::getInstance()->getEnzyme(enzyme);
      const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      const UInt target_decoy_index = registry.getIndex(Constants::UserParam::TARGET_DECOY); // UInt(-1): exists nowhere
      const UInt delta_mass_index = registry.getIndex("DeltaMass");
      const UInt isotope_error_index = registry.getIndex(Constants::UserParam::ISOTOPE_ERROR);
      static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
      const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
      // Declared without OpenMP, too: with the -fopenmp-simd fallback the compiler still parses the clauses below.
      int num_parts = 1;
#ifdef _OPENMP
      num_parts = static_cast<int>(std::clamp<SignedSize>(n / 512, 1, omp_get_max_threads()));
#endif
      const Size num_buckets = static_cast<Size>(num_parts);
      struct Part
      {
        TopHitSummary summary;
        std::string sequences;                                             // the sequence strings, one after the other
        std::vector<std::vector<std::pair<Size, Size>>> peptide_buckets;   // (offset, length) in sequences, by hash
        std::vector<std::vector<const std::string*>> protein_buckets;      // accessions of the peptide evidences, by hash
      };
      std::vector<Part> parts(num_parts);
      std::vector<std::exception_ptr> errors(num_parts); // an exception must not leave the parallel region
#pragma omp parallel for schedule(static, 1) num_threads(num_parts)
      for (int part = 0; part < num_parts; ++part)
      try
      {
        Part& p = parts[part];
        TopHitSummary& s = p.summary;
        p.peptide_buckets.resize(num_buckets);
        p.protein_buckets.resize(num_buckets);
        const std::hash<std::string_view> hash;
        EnzymaticDigestion digestor;
        digestor.setEnzyme(protease);
        std::unordered_map<std::string, Size> missed_cleavages_of; // by unmodified sequence
        std::vector<double> residue_masses; // of a hit, cumulative (prefix sums, for the b and y ions)
        for (SignedSize i = n * part / num_parts; i < n * (part + 1) / num_parts; ++i)
        {
          const PeptideIdentification& pid = peptide_ids[i];
          if (pid.getHits().empty()) { continue; }
          ++s.matched;
          const PeptideHit& top = pid.getHits().front();
          const std::string sequence = top.getSequence().toString();
          p.peptide_buckets[hash(sequence) % num_buckets].emplace_back(p.sequences.size(), sequence.size());
          p.sequences += sequence;
          for (const auto& ev : top.getPeptideEvidences())
          {
            const std::string& accession = ev.getProteinAccession();
            p.protein_buckets[hash(accession) % num_buckets].push_back(&accession);
          }
          if (top.getCharge() > 0) { ++s.charges[top.getCharge()]; }
          // the target/decoy type as PeptideHit::getTargetDecoyType() (by the index of the meta value; case-insensitive;
          // another value is neither, where getTargetDecoyType() throws)
          std::string td;
          if (top.metaValueExists(target_decoy_index)) { td = StringUtils::toLowered(top.getMetaValue(target_decoy_index).toString()); }
          const bool is_target = td == "target" || td == "target+decoy";
          if (with_errors)
          {
            // Precursor error: skip hits with unresolved charge (0) - getMZ() would throw - and, in an open search, those
            // with a delta mass (a modification the sequence lacks, not a mass error; as OpenSearchModificationAnalysis
            // tells unmodified PSMs).
            const bool shifted = top.metaValueExists(delta_mass_index)
                                 && !OpenSearchModificationAnalysis::isUnmodified(static_cast<double>(top.getMetaValue(delta_mass_index)));
            if (is_target && !shifted && top.getCharge() > 0)
            {
              // the isotope offset the search matched at (theoretical minus observed)
              int offset = 0;
              if (top.metaValueExists(isotope_error_index))
              {
                const int reported = static_cast<int>(top.getMetaValue(isotope_error_index));
                offset = isotope_observed_minus_theoretical ? -reported : reported;
              }
              const double corrected_mz = monoisotopicPrecursorMZ(pid.getMZ(), offset, static_cast<double>(top.getCharge()));
              s.precursor_errors.push_back(Math::getPPM(corrected_mz, top.getSequence().getMZ(top.getCharge())));
            }
            // Fragment errors: of the annotated b and y peaks ("b3+", "y12++"; not those with a neutral loss or of other
            // ions), against their m/z from the hit's sequence: (prefix or suffix residue masses + terminal modification
            // [+ H2O] + charge protons) / charge, the masses TheoreticalSpectrumGenerator gives the annotated ions
            if (is_target && !top.getPeakAnnotations().empty())
            {
              const AASequence& sequence = top.getSequence();
              const Size length = sequence.size();
              residue_masses.assign(length + 1, 0.0);
              for (Size k = 0; k < length; ++k)
              {
                residue_masses[k + 1] = residue_masses[k] + sequence[k].getMonoWeight(Residue::Internal);
              }
              const double n_term = sequence.hasNTerminalModification() ? sequence.getNTerminalModification()->getDiffMonoMass() : 0.0;
              const double c_term = sequence.hasCTerminalModification() ? sequence.getCTerminalModification()->getDiffMonoMass() : 0.0;
              bool any = false;
              for (const PeptideHit::PeakAnnotation& pa : top.getPeakAnnotations())
              {
                const std::string& name = pa.annotation;
                bool prefix = false;
                Size ordinal = 0;
                if (name.empty() || (name[0] != 'b' && name[0] != 'y') || pa.charge < 1
                    || !FragmentIonLikelihoodModel::parseIonName(name, prefix, ordinal) || ordinal >= length
                    || name.find_first_not_of('+', name.find_first_not_of("0123456789", 1)) != std::string::npos)
                {
                  continue; // another ion, or one with a neutral loss ("y3-H2O1+")
                }
                const double neutral = prefix ? residue_masses[ordinal] + n_term
                                              : residue_masses[length] - residue_masses[length - ordinal] + c_term + water_mass;
                const double theoretical_mz = (neutral + pa.charge * Constants::PROTON_MASS_U) / pa.charge;
                s.fragment_errors.push_back(Math::getPPM(pa.mz, theoretical_mz));
                any = true;
              }
              if (any) { ++s.fragment_errors_psms; }
            }
          }
          // counted once per sequence of the range: the regular expression takes a block from boost's process-wide
          // memory cache per search, which the threads contended for (93 ms at 16 and 64 threads for a timsTOF run)
          std::string unmodified = top.getSequence().toUnmodifiedString();
          auto cleavages = missed_cleavages_of.find(unmodified);
          if (cleavages == missed_cleavages_of.end())
          {
            const Size count = digestor.countInternalCleavageSites(unmodified);
            cleavages = missed_cleavages_of.emplace(std::move(unmodified), count).first;
          }
          ++s.missed_cleavages[cleavages->second];
          if (td == "decoy") { ++s.n_decoy; }
          else if (is_target) { ++s.n_target; }
        }
      }
      catch (...)
      {
        errors[part] = std::current_exception();
      }
      for (const std::exception_ptr& error : errors)
      {
        if (error) { std::rethrow_exception(error); }
      }

      // The distinct strings of each bucket (equal strings have equal hashes, so they share a bucket), and the buffers of
      // each range freed by a thread of the team
      std::vector<Size> unique_peptides(num_buckets, 0), unique_proteins(num_buckets, 0);
      std::vector<std::exception_ptr> bucket_errors(num_buckets);
#pragma omp parallel num_threads(num_parts)
      {
        std::vector<std::string_view> strings;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(num_buckets); ++sb)
        try
        {
          const Size b = static_cast<Size>(sb);
          strings.clear();
          for (const Part& p : parts)
          {
            for (const auto& [offset, length] : p.peptide_buckets[b]) { strings.emplace_back(p.sequences.data() + offset, length); }
          }
          std::sort(strings.begin(), strings.end());
          unique_peptides[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
          strings.clear();
          for (const Part& p : parts)
          {
            for (const std::string* accession : p.protein_buckets[b]) { strings.emplace_back(*accession); }
          }
          std::sort(strings.begin(), strings.end());
          unique_proteins[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
        }
        catch (...)
        {
          bucket_errors[sb] = std::current_exception();
        }
#pragma omp for schedule(static, 1)
        for (int part = 0; part < num_parts; ++part)
        {
          std::string().swap(parts[part].sequences);
          decltype(Part::peptide_buckets)().swap(parts[part].peptide_buckets);
          decltype(Part::protein_buckets)().swap(parts[part].protein_buckets);
        }
      }

      for (const std::exception_ptr& error : bucket_errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
      Size num_precursor_errors = 0, num_fragment_errors = 0;
      for (const Part& p : parts)
      {
        num_precursor_errors += p.summary.precursor_errors.size();
        num_fragment_errors += p.summary.fragment_errors.size();
      }
      TopHitSummary all = std::move(parts.front().summary);
      all.unique_peptides = std::accumulate(unique_peptides.begin(), unique_peptides.end(), Size(0));
      all.unique_proteins = std::accumulate(unique_proteins.begin(), unique_proteins.end(), Size(0));
      all.precursor_errors.reserve(num_precursor_errors);
      all.fragment_errors.reserve(num_fragment_errors);
      for (int part = 1; part < num_parts; ++part)
      {
        TopHitSummary& s = parts[part].summary;
        all.matched += s.matched;
        all.n_target += s.n_target;
        all.n_decoy += s.n_decoy;
        all.fragment_errors_psms += s.fragment_errors_psms;
        for (const auto& [charge, count] : s.charges) { all.charges[charge] += count; }
        for (const auto& [cleavages, count] : s.missed_cleavages) { all.missed_cleavages[cleavages] += count; }
        all.precursor_errors.insert(all.precursor_errors.end(), s.precursor_errors.begin(), s.precursor_errors.end());
        all.fragment_errors.insert(all.fragment_errors.end(), s.fragment_errors.begin(), s.fragment_errors.end());
      }
      return all;
    }
    // Math::median() and Math::MAD() of @p values (reordered), by selection instead of sorts: the same order statistics,
    // so the same values (the sorts of the 250,000 errors of a timsTOF run took 20-30 ms each). Unless a value is not
    // finite (a NaN is not ordered; an infinite median gives NaN deviations) or a negative zero (equal to 0, but written
    // differently), where the result may depend on the algorithm or on the order the sort left: then Math itself.
    std::pair<double, double> medianAndMADBySelection(std::vector<double>& values)
    {
      const auto ambiguous = [](double v) { return !std::isfinite(v) || (v == 0.0 && std::signbit(v)); };
      if (values.empty() || std::any_of(values.begin(), values.end(), ambiguous))
      {
        const double median = Math::median(values.begin(), values.end());
        return {median, Math::MAD(values.begin(), values.end(), median)};
      }
      // the median of finite values without a negative zero; their deviations are not negative zeros and not NaN
      // Large lists (the b/y fragment errors of a timsTOF run: millions) with the OpenMP threads: the selection by
      // libstdc++'s parallel nth_element (the same k-th smallest value), the deviations in parallel.
      constexpr Size parallel_from = 65536;
      const auto median_of = [](std::vector<double>& v)
      {
        const auto upper = v.begin() + static_cast<SignedSize>(v.size() / 2);
#ifdef PROSE_GNU_PARALLEL
        if (v.size() >= parallel_from && omp_get_max_threads() > 1) { __gnu_parallel::nth_element(v.begin(), upper, v.end()); }
        else
#endif
        { std::nth_element(v.begin(), upper, v.end()); }
        if (v.size() % 2 == 1) { return *upper; }
        return (*std::max_element(v.begin(), upper) + *upper) / 2.0;
      };
      const double median = median_of(values);
      std::vector<double> diffs(values.size());
#pragma omp parallel for schedule(static) if (values.size() >= parallel_from)
      for (SignedSize i = 0; i < static_cast<SignedSize>(values.size()); ++i) { diffs[i] = fabs(values[i] - median); }
      return {median, median_of(diffs)};
    }
  }

  // =====================================================================
  // Recompute result-level stats from a FINAL (post-rescoring/post-FDR) PSM list.
  // =====================================================================
  namespace
  {
    // Spectra with at least one PSM. The PSMs of co-isolated precursors (chimeric:enabled) belong to the spectrum
    // (scan_index) of the selected precursor, which counts once.
    Size countMatchedSpectra(const PeptideIdentificationList& peptide_ids)
    {
      Size matched = 0;
      std::set<unsigned int> selected_scans, coisolated_scans;
      for (const PeptideIdentification& pid : peptide_ids)
      {
        if (pid.getHits().empty()) continue;
        const PeptideHit& top = pid.getHits().front();
        const bool coisolated = top.metaValueExists(ProSEAlgorithm::CHIMERIC_RANK)
                                && static_cast<int>(top.getMetaValue(ProSEAlgorithm::CHIMERIC_RANK)) == 2;
        if (!pid.metaValueExists("scan_index"))
        {
          if (!coisolated) ++matched;
          continue;
        }
        const unsigned int scan = static_cast<unsigned int>(pid.getMetaValue("scan_index"));
        (coisolated ? coisolated_scans : selected_scans).insert(scan);
        if (!coisolated) ++matched;
      }
      for (unsigned int scan : coisolated_scans)
      {
        if (selected_scans.count(scan) == 0) ++matched;
      }
      return matched;
    }

    // countMatchedSpectra(), given @p with_hits, the number of identifications with hits: the same unless a top hit is
    // of a co-isolated precursor (chimeric:enabled), which is checked first (in parallel: no set of scans per search)
    Size matchedSpectra(const PeptideIdentificationList& peptide_ids, Size with_hits)
    {
      const UInt rank_index = MetaInfoInterface::metaRegistry().getIndex(ProSEAlgorithm::CHIMERIC_RANK);
      if (rank_index == UInt(-1)) { return with_hits; } // (a name that was never registered is on no hit)
      bool any_coisolated = false;
#pragma omp parallel for schedule(static) reduction(||: any_coisolated)
      for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
      {
        const auto& hits = peptide_ids[i].getHits();
        if (!hits.empty() && hits.front().metaValueExists(rank_index)
            && static_cast<int>(hits.front().getMetaValue(rank_index)) == 2) { any_coisolated = true; }
      }
      return any_coisolated ? countMatchedSpectra(peptide_ids) : with_hits;
    }
  }

  void ProSEAlgorithm::updateFinalStats(RunStatistics& stats,
                                        const PeptideIdentificationList& peptide_ids,
                                        const std::string& enzyme,
                                        bool fdr_applied)
  {
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    const TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme, false);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    stats.charge_histogram = summary.charges;
    stats.missed_cleavage_histogram = summary.missed_cleavages;
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    stats.fdr_applied = fdr_applied;
    stats.achieved_psm_fdr = fdr_applied ? maxRetainedScore_(peptide_ids) : -1.0;
  }

  void ProSEAlgorithm::releaseInParallel(SearchResult& result)
  {
    for (auto& protein_id : result.protein_ids) { freeInParallel(protein_id.getHits()); }
    std::vector<ProteinIdentification>().swap(result.protein_ids);
    freeInParallel(result.peptide_ids);
    freeInParallel(result.pre_fdr_peptide_ids);
  }

  // =====================================================================
  // Helper: maximum top-hit score among retained PSMs (== achieved FDR after
  // q-value filtering). Returns -1.0 for an empty set.
  // =====================================================================
  double ProSEAlgorithm::maxRetainedScore_(const PeptideIdentificationList& peptide_ids)
  {
    double max_score = -1.0;
    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      max_score = std::max(max_score, pid.getHits().front().getScore());
    }
    return max_score;
  }

  // =====================================================================
  // Helper: fill per-run identification statistics (silent — no logging).
  // Computes target/decoy counts from the (post-FDR) hits it is given, so
  // SearchResult::stats is consistent for any caller; leaves achieved FDR and
  // timing fields untouched (captured at well-defined points in the search paths).
  // =====================================================================
  ProSEAlgorithm::SpectraCounts_ ProSEAlgorithm::countSpectra_(const PeakMap& spectra) const
  {
    // in parallel: one pass over the spectra of a timsTOF run (600,000) took 0.13 s on one thread (memory bound)
    const bool search_outside = param_.getValue("precursor:recorded_charge_outside_range").toString() == "search";
    const int min_charge = static_cast<int>(precursor_min_charge_), max_charge = static_cast<int>(precursor_max_charge_);
    Size ms2 = 0, profile = 0, out_of_range = 0;
#pragma omp parallel reduction(+: ms2, profile, out_of_range) if (spectra.size() >= 16384)
    {
      std::vector<uint16_t> charges;
#pragma omp for schedule(static)
      for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
      {
        const MSSpectrum& s = spectra[i];
        if (s.getMSLevel() != 2) continue;
        ++ms2;
        if (s.getType() == SpectrumSettings::SpectrumType::PROFILE) ++profile;
        if (s.getPrecursors().empty()) continue;
        FragmentIndex::precursorCharges(s.getPrecursors()[0], min_charge, max_charge, charges, search_outside);
        if (charges.empty()) ++out_of_range;
      }
    }
    SpectraCounts_ counts;
    counts.ms2 = ms2;
    counts.profile = profile;
    counts.charge_out_of_range = out_of_range;
    return counts;
  }

  void ProSEAlgorithm::preprocessSearchedSpectra_(PeakMap& spectra, const std::string& input,
                                                   FragmentIonLikelihoodModel::PeakLists& ion_evidence,
                                                   PeakMap* evidence_spectra, PeakMap* query_spectra) const
  {
    const SpectraCounts_ counts = countSpectra_(spectra);
    const std::string of = input.empty() ? std::string() : " of '" + input + "'";
    if (counts.ms2 == 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << (input.empty() ? std::string("the input") : "'" + input + "'")
                      << " holds no MS2 spectra (e.g. an MS1-only file): "
                      << "nothing is searched and no peptides are identified." << std::endl;
    }
    if (counts.profile > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.profile << " of " << counts.ms2 << " MS2 spectra" << of
                      << " are recorded as profile spectra. ProSE expects centroided spectra and searches profile "
                      << "spectra as they are, which identifies few peptides, if any. Centroid them first (e.g. with "
                      << "PeakPickerHiRes, or peak picking when converting the raw files)." << std::endl;
    }
    if (counts.charge_out_of_range > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.charge_out_of_range << " of " << counts.ms2
                      << " MS2 spectra" << of << " have a precursor charge outside of precursor:min_charge.."
                      << "max_charge (" << precursor_min_charge_ << ".." << precursor_max_charge_
                      << ") and are not searched." << std::endl;
    }
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm", deisotope_requested_,
                       peaks_keep_n_, peaks_window_top_, peaks_window_type_, deisotoping_,
                       self_trained_ion_priors_ ? &ion_evidence : nullptr, ion_prior_scored_peaks_, evidence_spectra,
                       query_spectra);
  }

  void ProSEAlgorithm::collectRunStatistics_(
      const PeakMap& spectra,
      const std::vector<ProteinIdentification>& /*protein_ids*/,
      const PeptideIdentificationList& peptide_ids,
      RunStatistics& stats) const
  {
    const SpectraCounts_ counts = countSpectra_(spectra);
    stats.ms2_spectra = counts.ms2;
    stats.profile_spectra = counts.profile;
    stats.charge_out_of_range_spectra = counts.charge_out_of_range;
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    // Target/decoy of the FINAL (post-FDR) hits, so SearchResult::stats stays consistent
    // for library callers that read it without calling updateFinalStats().
    // Mass errors for the tolerance estimate: see summarizeTopHits() (signed, isotope-corrected, of the target PSMs).
    TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme_, true, isotope_error_observed_minus_theoretical_);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    for (const auto& [charge, count] : summary.charges) { stats.charge_histogram[charge] += count; }
    for (const auto& [cleavages, count] : summary.missed_cleavages) { stats.missed_cleavage_histogram[cleavages] += count; }
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    vector<double>& precursor_errors = summary.precursor_errors;
    vector<double>& fragment_errors = summary.fragment_errors;
    const Size fragment_errors_psms = summary.fragment_errors_psms;

    // -- Per-run tolerance estimate: the median of the signed errors (the bias) and their MAD (the spread; unscaled).
    //    The recommended (symmetric) tolerance covers the bias plus three standard deviations of a normal distribution,
    //    estimated robustly as 1.4826 * MAD: ceil(|median| + 3 * 1.4826 * MAD). Of the fragments it is an estimate
    //    within the searched tolerance: peaks beyond it were not matched. --
    const Size min_psms_for_estimation = 10;
    const auto recommended_tolerance = [](double median, double mad) { return std::ceil(std::fabs(median) + 3.0 * 1.4826 * mad); };
    if (precursor_errors.size() >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(precursor_errors); // = Math::median() / Math::MAD()
      stats.prec_err_median = med;
      stats.prec_err_mad = mad;
      stats.prec_err_recommended = recommended_tolerance(med, mad);
      stats.prec_tol_valid = true;
    }
    if (fragment_errors_psms >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(fragment_errors);
      stats.frag_err_median = med;
      stats.frag_err_mad = mad;
      stats.frag_err_recommended = recommended_tolerance(med, mad);
      stats.frag_tol_valid = true;
    }
  }

  // =====================================================================
  // Self-trained ion priors: learn presence, intensity rank and mass error
  // of fragment ions from the confident PSMs of this run, one model per
  // half of the spectra, and annotate every hit with the other half's model.
  // =====================================================================
  // static
  AASequence ProSEAlgorithm::reversedNoiseSequence_(const AASequence& sequence)
  {
    // Reverse all but the C-terminal residue and keep every residue's modification, so the noise
    // hypothesis shares composition, mass and enzymatic C-terminus with the peptide.
    const Size n = sequence.size();
    if (n < 3) return sequence;
    AASequence reversed;
    for (Size i = n - 1; i-- > 0;) { reversed += &sequence[i]; }
    reversed += &sequence[n - 1];
    if (sequence.hasNTerminalModification()) { reversed.setNTerminalModification(sequence.getNTerminalModification()); }
    if (sequence.hasCTerminalModification()) { reversed.setCTerminalModification(sequence.getCTerminalModification()); }
    return reversed;
  }

  void ProSEAlgorithm::annotateIonPriors_(const PeakMap& spectra,
                                          const FragmentIonLikelihoodModel::PeakLists& evidence,
                                          std::vector<ProteinIdentification>& protein_ids,
                                          PeptideIdentificationList& peptide_ids) const
  {
    if (!self_trained_ion_priors_) return;
    StopWatch sw_train, sw_annotate;
    sw_train.start();
    const bool ppm = fragment_mass_tolerance_unit_ == "ppm";
    const double tolerance = fragment_mass_tolerance_;
    const SpectrumGenerators_ generators = spectrumGenerators_();
    const bool rich = ion_prior_model_ == FragmentIonLikelihoodModel::ContextSet::RICH;
    // Fragment charges of the theoretical ions: b/y ions of a z+ precursor carry at most z - 1 charges.
    const int max_fragment_charge = ion_prior_max_fragment_charge_;
    const auto fragment_charges = [max_fragment_charge](int precursor_charge) {
      return std::max(1, std::min(max_fragment_charge, precursor_charge - 1));
    };
    // The ions of @p sequence (charge 1 to fragment_charges(charge)) matched against peak list @p scan: from the
    // series, ordinals and charges of getIonSeriesSpectrum() if the generator supports it (no ion names to write and
    // parse back: the same ions), else from getSpectrum() with the names.
    struct IonScratch
    {
      PeakSpectrum theo;
      std::vector<uint8_t> series;
      std::vector<UInt> ordinals;
      std::vector<Int> charges;
    };
    const auto match_ions = [&](const FragmentIonLikelihoodModel& model, const TheoreticalSpectrumGenerator& tsg, Size scan,
                                const AASequence& sequence, int charge, IonScratch& scratch, std::vector<FragmentIonLikelihoodModel::Ion>& ions)
    {
      if (tsg.getIonSeriesSpectrum(scratch.theo, scratch.series, scratch.ordinals, scratch.charges, sequence, 1, fragment_charges(charge)))
      {
        model.matchIons(evidence, scan, scratch.theo, scratch.series, scratch.ordinals, scratch.charges, sequence, charge, tolerance, ppm, ions);
        return;
      }
      scratch.theo.clear(true);
      tsg.getSpectrum(scratch.theo, sequence, 1, fragment_charges(charge));
      model.matchIons(evidence, scan, scratch.theo, sequence, charge, tolerance, ppm, ions);
    };

    // 1. The spectrum of every PSM (scan_index of postProcessHits_); -1 if it has no peaks to match. In parallel, with
    //    the meta values read by registry index (a lookup by name takes the registry's process-wide lock).
    const Size n_ids = peptide_ids.size();
    std::vector<SignedSize> scans(n_ids, -1);
    const UInt mv_scan_index = MetaInfoInterface::metaRegistry().getIndex("scan_index"); // UInt(-1): never registered
    const UInt mv_target_decoy = MetaInfoInterface::metaRegistry().getIndex(Constants::UserParam::TARGET_DECOY);
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      const PeptideIdentification& pi = peptide_ids[i];
      if (pi.getHits().empty() || !pi.metaValueExists(mv_scan_index)) continue;
      const SignedSize scan = static_cast<int>(pi.getMetaValue(mv_scan_index));
      if (scan >= 0 && static_cast<Size>(scan) < spectra.size() && static_cast<Size>(scan) < evidence.size()
          && evidence.peaks(static_cast<Size>(scan)) > 0)
      {
        scans[i] = scan;
      }
    }
    // Best hit by the identification's score orientation (postProcessHits_ sorts the hits, but the
    // annotation must not depend on it).
    const auto best_hit = [](const PeptideIdentification& pi) -> const PeptideHit& {
      const std::vector<PeptideHit>& hits = pi.getHits();
      const bool higher_better = pi.isHigherScoreBetter();
      Size best = 0;
      for (Size h = 1; h < hits.size(); ++h)
      {
        if (higher_better ? hits[h].getScore() > hits[best].getScore() : hits[h].getScore() < hits[best].getScore()) best = h;
      }
      return hits[best];
    };

    // 2. Training PSMs of each half (fold = scan parity): its best target hits at TDC q <= threshold over the
    //    native score, computed within the half, so the other half's labels are not used. Decoys win exact
    //    ties, as in the native yield evaluation, which keeps the selection conservative. The search has decoys
    //    (target-only searches never get here), so (D + 1) / T is an estimate even in a half without a decoy hit.
    //    The candidate rows of the PSMs are made in parallel (label: -1 none, 0 target, 1 decoy) and collected in PSM
    //    order; each half is sorted by a total order (the index decides the last ties), so the order does not depend
    //    on the sort, and the two halves are sorted at the same time.
    struct Row { double score; bool decoy; Size index; };
    std::vector<Row> candidate(n_ids);
    std::vector<signed char> label(n_ids, -1);
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      if (scans[i] < 0) continue;
      const PeptideHit& top = best_hit(peptide_ids[i]);
      if (!top.metaValueExists(mv_target_decoy)) continue;
      const bool decoy = top.getMetaValue(mv_target_decoy).toString() == "decoy";
      candidate[i] = {peptide_ids[i].isHigherScoreBetter() ? top.getScore() : -top.getScore(), decoy, static_cast<Size>(i)};
      label[i] = decoy ? 1 : 0;
    }
    std::array<std::vector<Row>, 2> rows;
    for (Size i = 0; i < n_ids; ++i)
    {
      if (label[i] >= 0) rows[static_cast<Size>(scans[i]) % 2].push_back(candidate[i]);
    }
    std::vector<Row>().swap(candidate);
    std::array<std::vector<Size>, 2> training;
    // (two threads, if this thread has them; declared without OpenMP, too: -fopenmp-simd parses the clause)
    [[maybe_unused]] int sort_threads = 1;
#ifdef _OPENMP
    sort_threads = std::min(2, omp_get_max_threads());
#endif
#pragma omp parallel for schedule(static, 1) num_threads(sort_threads)
    for (int fold = 0; fold < 2; ++fold)
    {
      std::vector<Row>& fold_rows = rows[fold];
      std::sort(fold_rows.begin(), fold_rows.end(), [](const Row& a, const Row& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.decoy != b.decoy) return a.decoy;
        return a.index < b.index;
      });
    }
    for (Size fold = 0; fold < 2; ++fold)
    {
      std::vector<Row>& fold_rows = rows[fold];
      std::vector<double> q(fold_rows.size(), 1.0);
      Size targets = 0, decoys = 0;
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (fold_rows[r].decoy) { ++decoys; }
        else { ++targets; }
        q[r] = static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1));
      }
      for (Size r = fold_rows.size(); r-- > 1;) { q[r - 1] = std::min(q[r - 1], q[r]); }
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (!fold_rows[r].decoy && q[r] <= ion_prior_train_fdr_) training[fold].push_back(fold_rows[r].index);
      }
    }
    const bool trained = training[0].size() >= ion_prior_min_psms_ && training[1].size() >= ion_prior_min_psms_
                         && std::isfinite(tolerance) && tolerance > 0.0;

    // 3. One model per half; its signal: the ions of the training PSMs, its noise: those of their reversed sequences
    //    on the same spectra. Ions are matched in parallel; the integer counts are added up in training order.
    std::array<FragmentIonLikelihoodModel, 2> models{FragmentIonLikelihoodModel(ion_prior_model_),
                                                     FragmentIonLikelihoodModel(ion_prior_model_)};
    if (trained)
    {
      for (Size fold = 0; fold < 2; ++fold)
      {
        const std::vector<Size>& psms = training[fold];
        FragmentIonLikelihoodModel& model = models[fold];
        std::vector<std::vector<FragmentIonLikelihoodModel::Ion>> signal(psms.size()), noise(psms.size());
#pragma omp parallel for schedule(dynamic, 16)
        for (SignedSize t = 0; t < static_cast<SignedSize>(psms.size()); ++t)
        {
          const Size scan = static_cast<Size>(scans[psms[t]]);
          const PeptideHit& hit = best_hit(peptide_ids[psms[t]]);
          const int charge = hit.getCharge();
          const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spectra[scan]);
          IonScratch scratch;
          match_ions(model, tsg, scan, hit.getSequence(), charge, scratch, signal[t]);
          const AASequence reversed = reversedNoiseSequence_(hit.getSequence());
          if (reversed == hit.getSequence()) continue; // no noise hypothesis (a palindrome)
          match_ions(model, tsg, scan, reversed, charge, scratch, noise[t]);
        }
        for (Size t = 0; t < psms.size(); ++t)
        {
          model.addObservations(signal[t], false);
          if (!noise[t].empty()) model.addObservations(noise[t], true);
        }
        model.finalize();
      }
    }
    else
    {
      OPENMS_LOG_WARN << "[ProSE] Ion priors: " << training[0].size() << " and " << training[1].size()
                      << " confident training PSMs in the two halves of the spectra, fewer than " << ion_prior_min_psms_
                      << "; the ion_prior_* features are 0 for this file." << std::endl;
    }
    sw_train.stop();

    // 4. Annotate every hit, scoring it with the model of the other half (cross-fitting). Untrained runs get zeros
    //    so the feature columns stay complete. The names of annotate:ion_prior_features are registered here, always in
    //    this order, before the parallel loop.
    sw_annotate.start();
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const bool write_llr = ion_prior_feature_llr_, write_explained = ion_prior_feature_explained_, write_topk = ion_prior_feature_topk_;
    const UInt mv_llr = write_llr ? registry.registerName(Constants::UserParam::ION_PRIOR_LLR) : 0;
    const UInt mv_explained = write_explained ? registry.registerName(Constants::UserParam::ION_PRIOR_EXPLAINED) : 0;
    const UInt mv_topk = write_topk ? registry.registerName(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED) : 0;
    Size annotated = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : annotated)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      PeptideIdentification& pi = peptide_ids[i];
      const SignedSize scan = scans[i];
      const bool scorable = trained && scan >= 0;
      const FragmentIonLikelihoodModel& model = models[scorable ? 1 - static_cast<Size>(scan) % 2 : 0];
      IonScratch scratch;
      std::vector<FragmentIonLikelihoodModel::Ion> ions;
      for (PeptideHit& hit : pi.getHits())
      {
        FragmentIonLikelihoodModel::Features features;
        if (scorable)
        {
          const int charge = hit.getCharge();
          match_ions(model, generators.forSpectrum(spectra[scan]), static_cast<Size>(scan), hit.getSequence(), charge, scratch, ions);
          features = model.score(ions);
        }
        if (write_llr) hit.setMetaValue(mv_llr, features.log_likelihood_ratio);
        if (write_explained) hit.setMetaValue(mv_explained, features.explained_presence);
        if (write_topk) hit.setMetaValue(mv_topk, features.top_predicted_observed);
        ++annotated;
      }
    }
    sw_annotate.stop();
    OPENMS_LOG_INFO << "[ProSE] Ion priors (" << (rich ? "rich" : "basic") << " model, "
                    << (ion_prior_scored_peaks_ ? "scored" : "all") << " peaks): "
                    << (trained ? "trained on " : "not trained; ") << training[0].size() << " + " << training[1].size()
                    << " PSMs at q <= " << ion_prior_train_fdr_ << " (halves by scan parity, each scored by the other's model) in "
                    << sw_train.getClockTime() << " s; " << annotated << " hits annotated in " << sw_annotate.getClockTime()
                    << " s; peak lists " << evidence.totalPeaks() << " peaks of " << evidence.size() << " spectra, "
                    << static_cast<double>(evidence.memoryUsage()) / (1024.0 * 1024.0) << " MiB." << std::endl;

    // 5. Percolator sees the new columns; the search parameters record every setting that determines the features
    //    (model, peaks, fragment charges, training threshold and minimum) and the training set.
    if (protein_ids.empty()) return;
    ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
    StringList features;
    if (search_parameters.metaValueExists("extra_features"))
    {
      const std::string existing = search_parameters.getMetaValue("extra_features").toString();
      if (!existing.empty()) features = ListUtils::create<std::string>(existing);
    }
    StringList written;
    if (write_llr) written.push_back(Constants::UserParam::ION_PRIOR_LLR);
    if (write_explained) written.push_back(Constants::UserParam::ION_PRIOR_EXPLAINED);
    if (write_topk) written.push_back(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED);
    for (const std::string& name : written)
    {
      if (std::find(features.begin(), features.end(), name) == features.end()) features.push_back(name);
    }
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
    search_parameters.setMetaValue("annotate:self_trained_ion_priors", "true");
    search_parameters.setMetaValue("annotate:ion_prior_features", ListUtils::concatenate(written, ","));
    search_parameters.setMetaValue("annotate:ion_prior_model", rich ? "rich" : "basic");
    search_parameters.setMetaValue("annotate:ion_prior_peaks", ion_prior_scored_peaks_ ? "scored" : "all");
    search_parameters.setMetaValue("annotate:ion_prior_max_fragment_charge", ion_prior_max_fragment_charge_);
    search_parameters.setMetaValue("annotate:ion_prior_train_fdr", ion_prior_train_fdr_);
    search_parameters.setMetaValue("annotate:ion_prior_min_psms", static_cast<int>(ion_prior_min_psms_));
    search_parameters.setMetaValue("ion_prior:trained", trained ? "true" : "false");
    search_parameters.setMetaValue("ion_prior:training_psms", static_cast<int>(training[0].size() + training[1].size()));
    search_parameters.setMetaValue("ion_prior:fold_training_psms",
                                   IntList{static_cast<int>(training[0].size()), static_cast<int>(training[1].size())});
  }

  // =====================================================================
  // Helper: run calibration pass on a subset of spectra
  // =====================================================================
  ProSEAlgorithm::CalibrationResult_
  ProSEAlgorithm::runCalibrationPass_(
      PeakMap& spectra,
      FragmentIndex& fragment_index,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const PeakMap* query_spectra) const
  {
    CalibrationResult_ result;

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    // Select subset by TIC (highest-quality spectra first)
    vector<pair<double, Size>> tic_index;
    tic_index.reserve(spectra.size());
    for (Size i = 0; i < spectra.size(); ++i)
    {
      double tic = 0;
      for (const auto& p : spectra[i]) { tic += p.getIntensity(); }
      tic_index.emplace_back(tic, i);
    }
    std::sort(tic_index.rbegin(), tic_index.rend());
    Size subset_size = std::max<Size>(1, static_cast<Size>(spectra.size() * calibration_subset_ratio_));
    if (subset_size > tic_index.size()) subset_size = tic_index.size();

    OPENMS_LOG_INFO << "[ProSE] Calibration: scoring " << subset_size << " / " << spectra.size()
                    << " spectra (top TIC)..." << std::endl;

    // Score subset and collect errors from the best hit per spectrum, with the ion series
    // the fragment index was built from (e.g. c/z+1 for ETD, where b/y would match nothing)
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // Collect per-spectrum best hits with scores and errors
    struct CalHit { double score; double prec_error; double frag_error; };
    vector<CalHit> cal_hits;
    // At most one hit per spectrum of the subset: the push_back() in the critical section below never reallocates, so
    // it cannot throw there (an exception must not leave an OpenMP critical section; the catch of the loop body is
    // outside it).
    cal_hits.reserve(subset_size);

    // Parallelize over the calibration subset, mirroring the main scoring loop
    // (scoreSpectraAgainstIndex_). Each iteration is independent: querySpectrum and
    // the shared TheoreticalSpectrumGenerators expose const, thread-safe methods (the
    // main loop already calls them concurrently), and every working variable below is
    // loop-local. The only cross-thread write is the push into cal_hits, guarded by a
    // critical section. Without this the calibration pass ran single-threaded: on a
    // many-core machine the 10% subset took about as long as the entire parallel
    // search, roughly doubling wall time for no extra useful work. Downstream is
    // order-independent — cal_hits is sorted by score and the error vectors are sorted
    // before quantiles — so parallel insertion order does not change the result.
    // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
    // leave the parallel region, which would terminate the process: the first one is rethrown after the loop.
    std::exception_ptr calibration_error;
    std::atomic<bool> calibration_failed{false};
#pragma omp parallel for schedule(dynamic)
    for (SignedSize si = 0; si < (SignedSize)subset_size; ++si)
    {
      if (calibration_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      const Size scan_idx = tic_index[si].second;
      const MSSpectrum& spec = spectra[scan_idx];
      const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);

      FragmentIndex::SpectrumMatchesTopN top_sms;
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_idx] : spec;
      fragment_index.querySpectrum(query, db, top_sms, generators.electronIons(spec));

      // Find the best-scoring hit for this spectrum
      double best_score = 0;
      AASequence best_seq;
      int best_isotope_error = 0;
      uint16_t best_charge = 0;
      float best_mean_error = 0;

      // Reused across this spectrum's candidates — same rationale as the main
      // scoring loop: a fresh PeakSpectrum per candidate churns its DataArrays
      // (ion names / charges filled by add_metainfo) on the heap.
      PeakSpectrum theo;

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& peptide = fragment_index.getPeptides()[sms.peptide_idx_];
        AASequence seq = fragment_index.isPeptideMassMode()
          ? fragment_index.reconstructRealizedSubSequence(peptide, db, sms.realized_length_, sms.subset_bitmask_)
          : fragment_index.reconstructModifiedSequence(peptide, db);
        // Clear peaks + data arrays before refilling; getSpectrum appends to
        // whatever is there. Its output is already sorted with add_metainfo=true.
        theo.clear(true);
        tsg.getSpectrum(theo, seq, 1, scoringMaxCharge_(sms.precursor_charge_));

        // The calibration PSMs are selected with the search's own score (upstream 825c33bb). The mass-accuracy score
        // favours PSMs with small fragment errors, so it narrows the fragment tolerance estimated from them below
        // (2-10% narrower than with HyperScore on 12 of 14 ppm runs); selecting by HyperScore instead did not
        // change the yield measurably.
        HyperScore::PSMDetail detail;
        double score = mass_accuracy_score_
          ? HyperScore::computeMassAccuracy(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, spec, theo, mass_error_sd_ppm_, detail)
          : HyperScore::computeWithDetail(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, spec, theo, detail);

        if (score > best_score)
        {
          best_score = score;
          best_seq = std::move(seq);
          best_isotope_error = sms.isotope_error_;
          best_charge = sms.precursor_charge_;
          best_mean_error = static_cast<float>(detail.mean_error);
        }
      }

      if (best_score == 0 || best_seq.empty()) continue;

      // Skip PSMs matched at a non-zero isotope offset. Their precursor m/z carries
      // an extra source of uncertainty (the +1/+2 isotope peak can be ambiguous in
      // the MS1 precursor picking), which inflates the variance of the calibration
      // quantile estimate. iso_err=0 PSMs are the gold-standard "true monoisotopic
      // peak picked" subset — the right population for estimating instrument bias.
      if (best_isotope_error != 0) continue;

      // Compute precursor error (signed), isotope-corrected. FragmentIndex searches
      // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs +
      // isotope_error * C13C12; the observed-to-monoiso m/z correction is
      //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
      // Matches the sign used by postProcessHits_'s PRECURSOR_ERROR_PPM annotation.
      double exp_mz = spec.getPrecursors()[0].getMZ();
      double theo_mz = best_seq.getMZ(best_charge);
      double corrected_exp_mz = monoisotopicPrecursorMZ(exp_mz, best_isotope_error, static_cast<double>(best_charge));
      double prec_err = (precursor_mass_tolerance_unit_ == "ppm")
                          ? Math::getPPM(corrected_exp_mz, theo_mz)
                          : (corrected_exp_mz - theo_mz);

#pragma omp critical (prose_calibration_hits)
      cal_hits.push_back({best_score, prec_err, static_cast<double>(best_mean_error)});
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_calibration_error)
        {
          if (!calibration_error) calibration_error = std::current_exception();
        }
        calibration_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (calibration_error) std::rethrow_exception(calibration_error);


    // Filter to high-confidence PSMs: keep only top 50% by score (robust against
    // random matches inflating the error distribution tails). Only crop if the
    // top 50% still meets the minimum PSM threshold; otherwise leave all hits
    // and let the calibration_min_psms_ check below disable calibration.
    Size cal_hits_total = cal_hits.size();
    std::sort(cal_hits.begin(), cal_hits.end(),
              [](const CalHit& a, const CalHit& b) { return a.score > b.score; });
    Size keep = cal_hits.size() / 2;
    if (keep >= calibration_min_psms_ && keep < cal_hits.size())
    {
      cal_hits.resize(keep);
    }

    // Extract error vectors from filtered hits.
    // Additionally discard PSMs whose signed precursor error falls outside the band the
    // FragmentIndex actually delivers. The mass-space window is [observed - lower, observed + upper];
    // rewriting in terms of the signed error e = observed - theoretical:
    //   theoretical in [observed - lower, observed + upper]
    //     <=>  e in [-upper, +lower]
    // So legitimate hits from FragmentIndex have e in [-upper, +lower], NOT [-lower, +upper].
    // The bounds flip for asymmetric windows; matching on the wrong band silently drops
    // real matches when the user has compensated an instrument bias.
    vector<double> precursor_errors;
    vector<double> fragment_errors_abs;
    for (const auto& h : cal_hits)
    {
      if (h.prec_error < -precursor_mass_tolerance_upper_ ||
          h.prec_error > precursor_mass_tolerance_lower_) continue; // wrong match
      precursor_errors.push_back(h.prec_error);
      if (h.frag_error > 0) fragment_errors_abs.push_back(h.frag_error);
    }

    if (precursor_errors.size() < calibration_min_psms_)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: insufficient PSMs (" << precursor_errors.size()
                      << " < " << calibration_min_psms_
                      << "), using configured tolerances." << std::endl;
      return result; // success=false
    }

    // --- Compute calibrated tolerances ---
    const double min_tolerance = 1e-6; // avoid non-positive tolerances

    // Signed precursor errors: 0.5% and 99.5% empirical quantiles give a distribution-free
    // 99% window (same approach as NuXL autotune, OpenNuXL.cpp:4939-4941). Asymmetric by
    // construction, no Gaussian assumption required — heavy tails and biased distributions
    // are handled correctly. We previously used residual_median + 3*residual_MAD which
    // captures only ~94% coverage on a Gaussian and less on heavy-tailed data; the
    // quantile method consistently produces more realistic windows.
    std::vector<double> sorted_errs = precursor_errors;
    std::sort(sorted_errs.begin(), sorted_errs.end());
    const size_t n = sorted_errs.size();
    // Bias estimate — median of the signed errors. Kept as a diagnostic field; the
    // quantile bounds below don't use it directly (they already encode the bias via
    // the asymmetry of [lo, hi]), but callers log it and tests assert on it.
    result.precursor_shift = Math::median(sorted_errs.begin(), sorted_errs.end(), /*sorted=*/true);
    const double lo = sorted_errs[static_cast<size_t>(n * 0.005)];                   // ~most negative
    const double hi = sorted_errs[std::min(n - 1, static_cast<size_t>(n * 0.995))];  // ~most positive
    // Convention (see lines 2193-2197): signed error e = observed - theoretical lies in
    // [-cal_upper, +cal_lower]. Matching the observed [lo, hi] against that gives:
    //   -cal_upper = lo  →  cal_upper = -lo
    //   +cal_lower = hi  →  cal_lower =  hi
    const double cal_lower_raw = std::max(min_tolerance, hi);
    const double cal_upper_raw = std::max(min_tolerance, -lo);
    // "Extreme" guards the degenerate case where the signed-error distribution
    // has essentially zero spread — e.g. a pure uniform shift fixture or a
    // single-peptide corner case — and the quantile bounds can't usefully inform
    // a calibrated window. Precursor calibration is then discarded; fragment
    // calibration still applies. A distribution strictly on one side of zero is
    // NOT extreme by this definition: we emit a legal half-line window
    // (cal_upper or cal_lower clamped to min_tolerance).
    //
    // Unit-aware threshold: 1 ppm is the right floor for ppm tolerances —
    // realistic proteomics-scale fixtures never trip it. In Da mode, 1.0 Da
    // would flag every realistic calibration (spreads are <0.05 Da), so fall
    // back to min_tolerance.
    const double extreme_bias_threshold =
        (precursor_mass_tolerance_unit_ == "ppm") ? 1.0 : min_tolerance;
    result.extreme_bias = (hi - lo) < extreme_bias_threshold;
    result.precursor_spread = std::max(cal_lower_raw, cal_upper_raw);  // diagnostic only
    if (!result.extreme_bias)
    {
      // Only tighten — cap against user-configured bounds.
      result.cal_lower = std::min(cal_lower_raw, precursor_mass_tolerance_lower_);
      result.cal_upper = std::min(cal_upper_raw, precursor_mass_tolerance_upper_);
    }
    // else: cal_lower/cal_upper stay at 0; writeback block skips the calibration result.

    // Fragment: 4 × 68th percentile (following NuXL / identipy convention).
    std::sort(fragment_errors_abs.begin(), fragment_errors_abs.end());
    if (!fragment_errors_abs.empty())
    {
      double frag_68 = fragment_errors_abs[static_cast<Size>((fragment_errors_abs.size() - 1) * 0.68)];
      result.fragment_tolerance = 4.0 * frag_68;
      if (result.fragment_tolerance < min_tolerance) result.fragment_tolerance = min_tolerance;
    }
    else
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    result.fragment_shift = 0.0;

    // Don't widen beyond configured fragment tolerance (only tighten).
    if (result.fragment_tolerance > fragment_mass_tolerance_)
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    result.success = true;

    OPENMS_LOG_INFO << "[ProSE] Calibration: " << precursor_errors.size() << " PSMs used (top "
                    << static_cast<int>(100.0 * precursor_errors.size() / cal_hits_total) << "% by score)" << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor signed: shift=" << std::fixed << std::setprecision(2) << result.precursor_shift
                    << " " << precursor_mass_tolerance_unit_ << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor spread: -> " << result.precursor_spread
                    << " " << precursor_mass_tolerance_unit_
                    << (result.extreme_bias ? " (extreme bias, discarded)" : "") << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Fragment tolerance:  " << fragment_mass_tolerance_
                    << " -> " << result.fragment_tolerance << " " << fragment_mass_tolerance_unit_ << std::endl;

    return result;
  }

  // =====================================================================
  // Serialize the end-of-search report to a YAML string. Hand-rolled (no YAML
  // library), so neither the TOPP tool nor the library needs an extra dependency
  // for this. Every string scalar is double-quoted (q()) so values containing ':'
  // (e.g. Windows paths) or other metacharacters cannot break the structure, and
  // non-finite numbers are emitted as null (num()).
  // =====================================================================
  std::string ProSEAlgorithm::renderRunSummaryYaml(
      const MultiFileSearchResult& mfres,
      const std::vector<std::pair<std::string, std::vector<std::string>>>& manifest,
      Size files_failed,
      Size files_total)
  {
    // Double-quoted YAML scalar: escape the two structural chars, the common
    // whitespace escapes, and every other C0 control + DEL as \xNN. A literal
    // control character is ill-formed in a double-quoted scalar and would corrupt
    // or invalidate the document — reachable via arbitrary input file paths.
    auto q = [](const std::string& s) {
      static const char* const hex = "0123456789ABCDEF";
      std::string out;
      out.reserve(s.size() + 2);
      out += '"';
      for (char ch : s)
      {
        const unsigned char c = static_cast<unsigned char>(ch);
        if (c == '\\' || c == '"') { out += '\\'; out += static_cast<char>(c); }
        else if (c == '\n') { out += "\\n"; }
        else if (c == '\t') { out += "\\t"; }
        else if (c == '\r') { out += "\\r"; }
        else if (c < 0x20 || c == 0x7F) { out += "\\x"; out += hex[c >> 4]; out += hex[c & 0x0F]; }
        else { out += static_cast<char>(c); }
      }
      out += '"';
      return out;
    };
    // Locale-independent number; non-finite -> null. YAML 1.1 reads an exponent
    // without a mantissa dot ("8e-05") as a string, so ensure a '.' before 'e'.
    auto num = [](double v) -> std::string {
      if (!std::isfinite(v)) { return "null"; }
      std::ostringstream o;
      o.imbue(std::locale::classic());
      o << v;
      std::string s = o.str();
      const auto e = s.find_first_of("eE");
      if (e != std::string::npos && s.find('.') == std::string::npos) { s.insert(e, ".0"); }
      return s;
    };
    auto bol = [](bool b) -> std::string { return b ? "true" : "false"; };

    std::ostringstream y;
    const SharedSearchStats& sh = mfres.shared;

    y << "shared:\n";
    y << "  database_file: " << q(sh.database_file) << "\n";
    y << "  enzyme: " << q(sh.enzyme) << "\n";
    y << "  precursor_tol_lower: " << num(sh.precursor_tol_lower) << "\n";
    y << "  precursor_tol_upper: " << num(sh.precursor_tol_upper) << "\n";
    y << "  precursor_tol_unit: " << q(sh.precursor_tol_unit) << "\n";
    y << "  fragment_tol: " << num(sh.fragment_tol) << "\n";
    y << "  fragment_tol_unit: " << q(sh.fragment_tol_unit) << "\n";
    y << "  min_charge: " << sh.min_charge << "\n";
    y << "  max_charge: " << sh.max_charge << "\n";
    y << "  missed_cleavages: " << sh.missed_cleavages << "\n";
    auto str_list = [&](const char* key, const std::vector<std::string>& items) {
      y << "  " << key << ":";
      if (items.empty()) { y << " []\n"; return; }
      y << "\n";
      for (const auto& it : items) { y << "    - " << q(it) << "\n"; }
    };
    str_list("fixed_mods", std::vector<std::string>(sh.fixed_mods.begin(), sh.fixed_mods.end()));
    str_list("variable_mods", std::vector<std::string>(sh.variable_mods.begin(), sh.variable_mods.end()));
    str_list("ion_series", std::vector<std::string>(sh.ion_series.begin(), sh.ion_series.end()));
    y << "  open_search: " << bol(sh.open_search) << "\n";
    y << "  calibration_enabled: " << bol(sh.calibration_enabled) << "\n";
    y << "  peptide_mass_mode: " << bol(sh.peptide_mass_mode) << "\n";
    y << "  chunked: " << bol(sh.chunked) << "\n";
    y << "  decoy_mode: " << q(sh.decoy_mode) << "\n";
    y << "  psm_fdr_threshold: " << num(sh.psm_fdr_threshold) << "\n";
    y << "  protein_fdr_threshold: " << num(sh.protein_fdr_threshold) << "\n";
    y << "  db_target_proteins: " << sh.db_target_proteins << "\n";
    y << "  db_decoy_proteins: " << sh.db_decoy_proteins << "\n";
    y << "  indexed_peptides: " << sh.indexed_peptides << "\n";
    y << "  indexed_fragments: " << sh.indexed_fragments << "\n";
    y << "  seconds_index_build: " << num(sh.seconds_index_build) << "\n";
    y << "  seconds_total: " << num(sh.seconds_total) << "\n";
    y << "  seconds_output: " << num(sh.seconds_output) << "\n";
    y << "  seconds_qpx_extended: " << num(sh.seconds_qpx_extended) << "\n";

    y << "per_file:";
    if (mfres.per_file.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& pf : mfres.per_file)
      {
        const RunStatistics& st = pf.stats;
        // The first key of each list item carries the "- " marker; the rest indent to 4.
        y << "  - input_file: " << q(st.input_file) << "\n";
        y << "    exit_code: " << static_cast<int>(pf.exit_code) << "\n";
        y << "    ms2_spectra: " << st.ms2_spectra << "\n";
        y << "    profile_spectra: " << st.profile_spectra << "\n";
        y << "    charge_out_of_range_spectra: " << st.charge_out_of_range_spectra << "\n";
        y << "    matched_spectra: " << st.matched_spectra << "\n";
        y << "    target_psms: " << st.target_psms << "\n";
        y << "    decoy_psms: " << st.decoy_psms << "\n";
        y << "    fdr_applied: " << bol(st.fdr_applied) << "\n";
        y << "    achieved_psm_fdr: " << num(st.achieved_psm_fdr) << "\n";
        y << "    unique_peptides: " << st.unique_peptides << "\n";
        y << "    unique_proteins: " << st.unique_proteins << "\n";
        y << "    hyperscore:\n";
        y << "      valid: " << bol(st.score_stats_valid) << "\n";
        y << "      min: " << num(st.hyperscore_min) << "\n";
        y << "      median: " << num(st.hyperscore_median) << "\n";
        y << "      max: " << num(st.hyperscore_max) << "\n";
        y << "    charge_histogram:";
        if (st.charge_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [z, c] : st.charge_histogram) { y << "      " << q(std::to_string(z)) << ": " << c << "\n"; } }
        y << "    missed_cleavage_histogram:";
        if (st.missed_cleavage_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [m, c] : st.missed_cleavage_histogram) { y << "      " << q(std::to_string(m)) << ": " << c << "\n"; } }
        y << "    precursor_error:\n";
        y << "      valid: " << bol(st.prec_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.prec_err_median) << "\n";
        y << "      mad_ppm: " << num(st.prec_err_mad) << "\n";
        y << "      recommended_ppm: " << num(st.prec_err_recommended) << "\n";
        y << "    fragment_error:\n";
        y << "      valid: " << bol(st.frag_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.frag_err_median) << "\n";
        y << "      mad_ppm: " << num(st.frag_err_mad) << "\n";
        y << "      recommended_ppm: " << num(st.frag_err_recommended) << "\n";
        y << "    timing_seconds:\n";
        y << "      calibration: " << num(st.seconds_calibration) << "\n";
        y << "      search: " << num(st.seconds_search) << "\n";
        y << "      fdr: " << num(st.seconds_fdr) << "\n";
      }
    }

    y << "outputs:";
    if (manifest.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& [label, paths] : manifest)
      {
        y << "  - type: " << q(label) << "\n";
        y << "    paths:";
        if (paths.empty()) { y << " []\n"; }
        else { y << "\n"; for (const auto& p : paths) { y << "      - " << q(p) << "\n"; } }
      }
    }
    y << "files_failed: " << files_failed << "\n";
    y << "files_total: " << files_total << "\n";

    return y.str();
  }

  // =====================================================================
  // Render the modification-discovery section (open search) to an ostream.
  // =====================================================================
  void ProSEAlgorithm::renderModificationSummary(
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      std::ostream& os)
  {
    const auto& dm_stats = mod_analysis.delta_mass_stats;
    const auto& ptm_stats = mod_analysis.ptm_stats;

    os << "[ProSE] ---------------- Modification discovery ---------------\n";
    os << "[ProSE] Delta mass : " << dm_stats.modified_psms << " modified / "
       << dm_stats.total_psms << " PSMs ("
       << std::fixed << std::setprecision(1)
       << (dm_stats.total_psms > 0 ? (100.0 * dm_stats.modified_psms / dm_stats.total_psms) : 0.0)
       << "%), median Δ=" << std::setprecision(4) << dm_stats.median_delta_mass
       << " Da, " << dm_stats.entries.size() << " bins\n";
    os << "[ProSE] PTMs      : " << ptm_stats.total_modified_psms << " PSMs with known PTMs, "
       << ptm_stats.unknown_modification_psms << " unknown, "
       << ptm_stats.num_unique_modifications << " unique PTMs\n";

    if (!ptm_stats.entries.empty())
    {
      os << "[ProSE]   Top PTMs (rank | name | count | % | mass Da):\n";
      size_t rank = 1;
      for (const auto& ptm : ptm_stats.entries)
      {
        if (rank > 15) { break; }
        std::string name = ptm.name;
        if (name.size() > 30) { name = StringUtils::substr(name, 0, 27) + "..."; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | "
           << std::setw(31) << std::left << name << std::right << " | "
           << std::setw(6) << ptm.count << " | "
           << std::setw(5) << std::fixed << std::setprecision(1) << ptm.percentage << " | "
           << std::setw(9) << std::fixed << std::setprecision(4) << ptm.theoretical_mass << "\n";
      }
    }

    std::vector<OpenSearchModificationAnalysis::DeltaMassEntry> unknown_dm;
    for (const auto& entry : dm_stats.entries)
    {
      if (!entry.is_known_modification && entry.count >= 5) { unknown_dm.push_back(entry); }
    }
    if (!unknown_dm.empty())
    {
      os << "[ProSE]   Top unknown delta masses (potential novel PTMs):\n";
      size_t rank = 1;
      for (const auto& dm : unknown_dm)
      {
        if (rank > 10) { break; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | Δ="
           << std::setw(11) << std::fixed << std::setprecision(4) << dm.delta_mass << " Da | "
           << std::setw(6) << dm.count << " PSMs | "
           << dm.unique_peptides << " peptides\n";
      }
    }
  }

  // =====================================================================
  // Render a human-readable single-run summary block.
  // =====================================================================
  std::string ProSEAlgorithm::renderTotalTime(const SharedSearchStats& sh)
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << "total " << sh.seconds_total << "s | search of the file(s) "
       << std::max(0.0, sh.seconds_total - sh.seconds_output) << "s | outputs " << sh.seconds_output << "s";
    if (sh.seconds_qpx_extended > 0.0) { os << " (of which -qpx_extended " << sh.seconds_qpx_extended << "s)"; }
    return os.str();
  }

  void ProSEAlgorithm::renderRunSummary(
      const RunStatistics& s,
      const SharedSearchStats& sh,
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      bool is_open_search,
      std::ostream& os)
  {
    auto join = [](const std::vector<std::string>& v) -> std::string {
      if (v.empty()) { return "(none)"; }
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? ", " : "") + v[i]; }
      return out;
    };
    auto join_int = [](const std::vector<std::string>& v) -> std::string {
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? "," : "") + v[i]; }
      return out.empty() ? std::string("(none)") : out;
    };

    if (!s.input_file.empty())
    {
      os << "[ProSE] Input        : " << s.input_file << "  (" << s.ms2_spectra << " MS2 spectra)\n";
    }
    if (s.ms2_spectra == 0)
    {
      os << "[ProSE] Warning      : no MS2 spectra in the input (e.g. an MS1-only file): nothing was searched\n";
    }
    if (s.profile_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.profile_spectra << " MS2 spectra are profile spectra; centroid them first\n";
    }
    if (s.charge_out_of_range_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.charge_out_of_range_spectra << " MS2 spectra not searched: precursor charge outside of "
         << sh.min_charge << "-" << sh.max_charge << "\n";
    }
    // -- Configuration recap --
    os << "[ProSE] Config       : " << sh.enzyme << ", " << sh.missed_cleavages << " MC | prec [-"
       << sh.precursor_tol_lower << ", +" << sh.precursor_tol_upper << "] " << sh.precursor_tol_unit
       << " | frag " << sh.fragment_tol << " " << sh.fragment_tol_unit
       << " | z " << sh.min_charge << "-" << sh.max_charge << "\n";
    os << "[ProSE]                fixed: " << join(sh.fixed_mods)
       << " | variable: " << join(sh.variable_mods) << "\n";
    os << "[ProSE]                ions: " << join_int(sh.ion_series)
       << " | calibration: " << (sh.calibration_enabled ? "on" : "off")
       << " | mode: " << (sh.open_search ? "open" : "closed")
       << (sh.peptide_mass_mode ? " | low_memory" : "")
       << (sh.chunked ? " | chunked" : "") << "\n";

    // -- Database / fragment index --
    os << "[ProSE] Database     : " << (sh.database_file.empty() ? std::string("(in-memory)") : std::string(sh.database_file))
       << "  (" << sh.db_target_proteins << " target";
    if (sh.db_decoy_proteins > 0) { os << " + " << sh.db_decoy_proteins << " decoy"; }
    os << " proteins, decoys: " << (sh.decoy_mode.empty() ? std::string("n/a") : std::string(sh.decoy_mode)) << ")\n";
    os << "[ProSE] Fragment idx : " << sh.indexed_peptides << " peptides / "
       << sh.indexed_fragments << " fragments  (build "
       << std::fixed << std::setprecision(1) << sh.seconds_index_build << " s)\n";

    // -- Results --
    os << "[ProSE] --------------------------- Results ---------------------------\n";
    os << "[ProSE] Matched      : " << s.matched_spectra << " / " << s.ms2_spectra << " spectra";
    if (s.ms2_spectra > 0)
    {
      os << "   (ID rate " << std::fixed << std::setprecision(1)
         << (100.0 * s.matched_spectra / s.ms2_spectra) << "%)";
    }
    os << "\n";
    os << "[ProSE] Peptides/prot: " << s.unique_peptides << " unique peptides | "
       << s.unique_proteins << " unique proteins\n";

    // FDR / target-decoy. Decide target-only strictly from the database, not from
    // a zero PSM count (a decoy DB with no decoy PSMs is still FDR-capable).
    if (sh.decoy_mode == "none (target-only)")
    {
      os << "[ProSE] FDR          : n/a (target-only database)\n";
    }
    else
    {
      os << "[ProSE] Target/decoy : " << s.target_psms << " target / " << s.decoy_psms << " decoy PSMs";
      if (s.fdr_applied)
      {
        os << "  | PSM FDR <= " << std::fixed << std::setprecision(1) << (sh.psm_fdr_threshold * 100.0) << "%";
        if (s.achieved_psm_fdr >= 0.0)
        {
          os << ", achieved " << std::setprecision(2) << (s.achieved_psm_fdr * 100.0) << "%";
        }
        else
        {
          os << " (0 PSMs retained)";
        }
      }
      else
      {
        os << "  | PSM FDR filtering: off";
      }
      os << "\n";
    }

    // Charge distribution
    if (!s.charge_histogram.empty())
    {
      os << "[ProSE] Charges      :";
      for (const auto& [z, c] : s.charge_histogram) { os << " " << z << ":" << c; }
      os << "\n";
    }
    // Missed cleavages
    if (!s.missed_cleavage_histogram.empty())
    {
      os << "[ProSE] Missed clv   :";
      for (const auto& [mc, c] : s.missed_cleavage_histogram) { os << " " << mc << ":" << c; }
      os << "\n";
    }
    // Score distribution
    if (s.score_stats_valid)
    {
      os << "[ProSE] HyperScore   : min " << std::fixed << std::setprecision(1) << s.hyperscore_min
         << "  median " << s.hyperscore_median << "  max " << s.hyperscore_max << "\n";
    }

    // -- Tolerance estimation --
    if (s.prec_tol_valid || s.frag_tol_valid)
    {
      os << "[ProSE] --------------------- Tolerance estimate ---------------------\n";
      if (s.prec_tol_valid)
      {
        os << "[ProSE] Precursor    : median " << std::fixed << std::setprecision(2) << s.prec_err_median
           << " ppm, MAD " << s.prec_err_mad << "  -> recommend " << static_cast<int>(s.prec_err_recommended) << " ppm\n";
      }
      if (s.frag_tol_valid)
      {
        os << "[ProSE] Fragment     : median " << std::fixed << std::setprecision(2) << s.frag_err_median
           << " ppm, MAD " << s.frag_err_mad << "  -> recommend " << static_cast<int>(s.frag_err_recommended) << " ppm\n";
      }
    }

    // -- Timing --
    os << "[ProSE] ------------------------- Timing -----------------------------\n";
    os << "[ProSE] " << std::fixed << std::setprecision(1);
    if (s.seconds_calibration > 0.0) { os << "calib " << s.seconds_calibration << "s | "; }
    if (s.seconds_search > 0.0)      { os << "scoring " << s.seconds_search << "s | "; }
    if (s.seconds_fdr > 0.0)         { os << "fdr " << s.seconds_fdr << "s | "; }
    os << "index(shared) " << sh.seconds_index_build << "s\n";
    // the whole run (set by the caller): the search of the file (reading it, rescoring, FDR, ...) and the outputs
    if (sh.seconds_total > 0.0) { os << "[ProSE] " << renderTotalTime(sh) << "\n"; }

    // -- Modification discovery (open search only) --
    if (is_open_search)
    {
      renderModificationSummary(mod_analysis, os);
    }
  }

} // namespace OpenMS
