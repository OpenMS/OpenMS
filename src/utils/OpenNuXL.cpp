// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

    // TODO: 
    // - precursor correction

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>

#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>
#include <OpenMS/ANALYSIS/RNPXL/MorpheusScore.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlFragmentAnnotationHelper.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlParameterParsing.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/regex.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>

#include <map>
#include <algorithm>
#include <iterator>

#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

#include <boost/accumulators/statistics/p_square_quantile.hpp>
using namespace boost::accumulators;

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//#define DEBUG_OpenNuXL 1
//#define OPENNUXL_SEPARATE_FEATURES 1
//#define DANGEROUS_FEAUTURES increases FDR on entrappment dataset
//#define FILTER_BAD_SCORES_ID_TAGS filter out some good hits
//#define FILTER_AMBIGIOUS_PEAKS  so far only worse results
//#define FILTER_NO_ARBITRARY_TAG_PRESENT
#define CALCULATE_LONGEST_TAG
//#define MODDS_ON_ABY_IONS_ONLY
//#define SVM_RECALIBRATE 1
//#define FILTER_RANKS 1
//#define CALCULATE_NUCLEOTIDE_TAGS 1
#define DONT_ACCUMULATE_PARTIAL_ION_SCORES 1

#define CONSIDER_AA_LOSSES 1

#define ANNOTATED_QUANTILES 1

#ifdef ANNOTATED_QUANTILES
typedef accumulator_set<double, stats<tag::p_square_quantile> > quantile_accu_t;

#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>

typedef accumulator_set<double, stats<tag::extended_p_square_quantile(quadratic)> > accumulator_t_quadratic;

struct SpectrumLevelScoreQuantiles
{
  SpectrumLevelScoreQuantiles():
    acc_(extended_p_square_probabilities = std::vector<double>{ 0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999, 1.00 })
  {   
  }
  void insert(double v) { acc_(v); }
  
  double quantileOfValue(double v)
  {
    double l = 0.0;
    double h = 1.0;
    double mid = l + (h - l) / 2.0; // we start with the median (0.5)
    double p_value = quantile(acc_, quantile_probability = mid); // value of quantile p

    size_t iter(0);
    while (fabs(p_value - v) > 0.01 && iter < 100)
    {
      mid = l + (h - l) / 2.0;
      p_value = quantile(acc_, quantile_probability = mid); // value of quantile p (e.g., 1234.56)
      if (p_value > v) // if the current quantile value (e.g., of the median) has a value larger than our value of interest
      {
        h = mid;  // then we need to search in the lower quantile range
      }
      else 
      {
        l = mid;
      }
      ++iter;
    }
    return l;
  }
  private:
    accumulator_t_quadratic acc_;
};

#endif

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;
// stores which residues (known to give rise to immonium ions) are in the sequence
struct ImmoniumIonsInPeptide
{
  explicit ImmoniumIonsInPeptide(const String& s)
  {
    for (const char & c : s)
    {
      switch (c)
      {
        case 'Y': Y = true; break;
        case 'W': W = true; break;
        case 'F': F = true; break;
        case 'H': H = true; break;
        case 'C': C = true; break;
        case 'P': P = true; break;
        case 'I':
        case 'L': L = true; break;
        case 'K': K = true; break;
        case 'M': M = true; break;
        case 'Q': Q = true; break;
        case 'E': E = true; break;
        default: break;
      }   
    } 
  }

  void swap()
  { // swap entries for decoy generation
    std::swap(Y,E);
    std::swap(W,Q);
    std::swap(F,M);
    std::swap(H,K);
    std::swap(C,L);
    std::swap(P,L);
  }
 
  bool Y = false;
  bool W = false;
  bool F = false;
  bool H = false;
  bool C = false;
  bool P = false;
  bool L = false;
  bool K = false;
  bool M = false;
  bool Q = false;
  bool E = false;
}; 

/// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
/// floats need to be initialized to zero as default
class AnnotatedHit
{
  public:
  /*
     Slim indices/views to lookup the actual sequence
   */
  StringView sequence;
  SignedSize peptide_mod_index = 0; // enumeration index of the non-NA peptide modification
  Size NA_mod_index = 0; // index of the NA modification

  int isotope_error = 0; // wheter the hit has been matched with isotopic misassignment

  static constexpr const char UNKNOWN_NUCLEOTIDE = '?';
  char cross_linked_nucleotide = UNKNOWN_NUCLEOTIDE;

  /**
       The main score (score) is a linear combination of (weighted) subscores

       For the fast score (ignoring all shifted peaks) we calculate: 
         score = 1.0 * total_loss_score 
               + 1.0 * total_MIC 
               + 0.333 * mass_error_p;
 
       For the all-ion score we calculate:
         peptides:
	      score = -6.486416409280039 
               + 4.059968526608637   * ah.total_MIC         
               + 0.5842539236790404  * ah.modds
               + 0.21721652155697285 * ah.total_loss_score
               + 1.9988345415208777  * ah.mass_error_p;
         XLs:
	       score = -6.648631037190969
               + 0.4688059636415974  * ah.Morph
               + 4.0386886051238     * ah.MIC         
               + 0.5446999629799386  * ah.modds
               + 0.25318342707227187 * ah.total_loss_score
               + 0.12472562244230834 * ah.partial_loss_score
               + 1.2107674392113372  * ah.mass_error_p
               + 2.3319284783288805  * ah.pl_MIC;
  */
  float score = 0;

  /**
       Normalized precursor mass error score.
       Mass error is assumed normally distributed with:
         - mean = 0
         - sd = sqrt(precursor_mass_tolerance) => variance = precursor tolerance
  */
  float mass_error_p = 0;

  //
  // Scores exclusively calculated from peaks without nucleotide shifts:
  //
  
  /**
      The total loss score is the X!Tandem HyperScore calculated from b-,y-ions 
      without any nucleotide shift.
  */
  float total_loss_score = 0;

  /**
      The matched ion current in immonium (immonium_score) and precursor ions (precursor_score) 
      without any nucleotide shift.

      see DOI: 10.1021/pr3007045 A Systematic Investigation into the Nature of Tryptic HCD Spectra
      imY = EmpiricalFormula("C8H10NO").getMonoWeight(); // 85%
      imW = EmpiricalFormula("C10H11N2").getMonoWeight(); // 84%
      imF = EmpiricalFormula("C8H10N").getMonoWeight(); // 84%
      imL = EmpiricalFormula("C5H12N").getMonoWeight(); // I/L 76%
      imH = EmpiricalFormula("C5H8N3").getMonoWeight(); // 70%
      imC = EmpiricalFormula("C2H6NS").getMonoWeight(); // CaC 61%
      imK1 = EmpiricalFormula("C5H13N2").getMonoWeight(); // 2%
      imP = EmpiricalFormula("C4H8N").getMonoWeight(); //?
      imQ = 101.0715; // 52%
      imE = 102.0555; // 37%
      imM = 104.0534; // 3%
  */
  float immonium_score = 0;
  float precursor_score = 0;

  /**
      The matched ion current (MIC), average fragment error (err), and morpheus score (Morph) are calculated 
      for b-,y-,a-ions without nucleotide shift. Morph is just the number of matched peaks + the fraction of MIC
  */
  float MIC = 0;
  float err = 0;
  float Morph = 0;

  /**
     The match odds (-log10) of observing this number of b-,a-, and y-ions assuming a uniform distribution of noise peaks.      
  */
  float modds = 0;

  //
  // Scores exclusively calculated from nucleotide shifted peaks:
  //
 
  /**
      The partial loss score is the X!Tandem HyperScore calculated from b-,a-, and y-ions 
      with nucleotide shifts. Matches from b- and a-ions are combined, i.e. a matching a_n-ion is counted as b_n-ion.
      For a precursor with charge N, all fragment ion charges up to N-1 are considered.

      Calculation of HyperScore:
      yFact = logfactorial_(y_ion_count);
      bFact = logfactorial_(b_ion_count);
      hyperScore = log1p(dot_product) + yFact + bFact;
  */
  float partial_loss_score = 0;

  /**
      The matched ion current (pl_MIC) of ladder ions, average fragment error (pl_err), and morpheus score (pl_Morph) are calculated 
      from b-,y-,a-ions with nucleotide shift.
      Morph: number of matched peaks + the fraction of MIC
  */
  float pl_MIC = 0;
  float pl_err = 0;
  float pl_Morph = 0;

  /*
     The match odds (-log10) of observing this number of b-,a-, and y-ions with nucleotide shifts assuming a uniform distribution of noise peaks.      
  */
  float pl_modds = 0;

  /*
     The MIC of precursor with all nucleotide shifts.
     Three variants: No additional loss, loss of water, and loss ammonia.
     Charge states considered: 1..N (precursor charge)
  */
  float pl_pc_MIC = 0;

  /**
      The matched ion current calculated from immonium ions with nucleotide shifts.
      Only singly charged immonium ions are considered.

      imY = EmpiricalFormula("C8H10NO").getMonoWeight();
      imW = EmpiricalFormula("C10H11N2").getMonoWeight();
      imF = EmpiricalFormula("C8H10N").getMonoWeight();
      imH = EmpiricalFormula("C5H8N3").getMonoWeight();
      imC = EmpiricalFormula("C2H6NS").getMonoWeight();
      imP = EmpiricalFormula("C4H8N").getMonoWeight();
      imL = EmpiricalFormula("C5H12N").getMonoWeight();
      imK1 = EmpiricalFormula("C5H13N2").getMonoWeight();
      imK2 = EmpiricalFormula("C5H10N1").getMonoWeight();
      imK3 = EmpiricalFormula("C6H13N2O").getMonoWeight();
      imQ = 101.0715;
      imE = 102.0555;
      imM = 104.0534;
   */
  float pl_im_MIC = 0;

  //
  // Scores calculated from peaks with AND without nucleotide shifts:
  //
  
  /**
       The complete TIC fraction of explained peaks (total_MIC) (excludes marker ions)
       For peptides: total_MIC = MIC + im_MIC + pc_MIC (b-,a-,y-ions, immonium ions, precursor ions)
       For XLs:      total_MIC = MIC + im_MIC + pc_MIC + pl_MIC + pl_pc_MIC + pl_im_MIC (same as above + shifted versions)
  */
  float total_MIC = 0;

  /**
       The matched ion current in marker ions (marker_ions_score) is not considered in scoring.
  */  
  float marker_ions_score = 0;

  /**
       Coverage of peptide by prefix or suffix ions (fraction)
       For example: PEPTIDER
                    01000100 (two of eight ions observed -> 2/8)       
       Shifted and non-shifted are combined to determine coverage.
  */
  float ladder_score = 0;
  /**
       Longest sequence covered in peptide by prefix or suffix ions (fraction).
       Coverage of peptide by prefix or suffix ions (fraction)
       For example: PEPTIDER
                    01110001 (three ions form the longest sequence -> 3/8)
       Shifted and non-shifted are combined to determine coverage.
  */
  float sequence_score = 0;

  float best_localization_score = 0;
  String localization_scores = 0;
  String best_localization;
  int best_localization_position = -1; // UNKNOWN
  std::vector<PeptideHit::PeakAnnotation> fragment_annotations;

  size_t tag_unshifted = 0;
  size_t tag_shifted = 0;
  size_t tag_XLed = 0;  // tag that contains the transition from unshifted to shifted

  double explained_peak_fraction = 0;
  double wTop50 = 0;

  size_t n_theoretical_peaks = 0;

  static bool hasBetterScore(const AnnotatedHit& a, const AnnotatedHit& b)
  {
    return a.score > b.score;
  }
};

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_OpenNuXL OpenNuXL 

    @brief Annotate NA to peptide crosslinks in MS/MS spectra.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_OpenNuXL.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_OpenNuXL.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class OpenNuXL :
  public TOPPBase
{
  bool fast_scoring_ = true; // fast or all fragment adduct scoring mode
  set<char> can_xl_; ///< nucleotides that can form cross-links

public:
  OpenNuXL() :
    TOPPBase("OpenNuXL", "Annotate RNA/DNA-peptide cross-links in MS/MS spectra.", false)
  {
  }

  static constexpr double MIN_HYPERSCORE = 0.1; // hit's with lower score than this will be neglected (usually 1 or 0 matches)
  static constexpr double MIN_TOTAL_LOSS_IONS = 1; // minimum number of matches to unshifted ions
  static constexpr double MIN_SHIFTED_IONS = 1; // minimum number of matches to shifted ions (applies to XLs only)
  static constexpr Size IA_CHARGE_INDEX = 0;
  static constexpr Size IA_RANK_INDEX = 1;
  static constexpr Size IA_DENOVO_TAG_INDEX = 2;
protected:
  /// percolator feature set
  StringList feature_set_;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML,raw"));
    registerInputFile_("NET_executable", "<executable>", "", "The .NET framework executable. Only required on linux and mac.", false, true, ListUtils::create<String>("skipexists"));
    registerInputFile_("ThermoRaw_executable", "<file>", "ThermoRawFileParser.exe", "The ThermoRawFileParser executable.", false, true, ListUtils::create<String>("skipexists"));

    registerInputFile_("database", "<file>", "", "input file ");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_tsv", "<file>", "", "tsv output file", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 6.0, "Precursor mass tolerance (+/- around precursor m/z).", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, false);

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0};
    registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak).", false, false);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 20.0, "Fragment mass tolerance (+/- around fragment m/z).", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment mass tolerance.", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'.", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>("Oxidation (M)"), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide.", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 6, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:max_size", "<num>", 1e6, "Maximum size a peptide may have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 2, "Number of missed cleavages.", false, false);

    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin/P", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);
    registerDoubleOption_("report:peptideFDR", "<num>", 0.0, "Maximum q-value of non-cross-linked peptides. (0 = disabled).", false, true);
    registerDoubleList_("report:xlFDR", "<num>", { 0.0 }, "Maximum q-value of cross-linked peptides. (0 = disabled). If multiple values are provided, multiple output files will be created.", false, true);

    registerInputFile_("percolator_executable", "<executable>", 
 // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "Percolator executable of the installation e.g. 'percolator.exe'", false, false, ListUtils::create<String>("skipexists"));


    // RNPxl specific
    registerTOPPSubsection_("RNPxl", "RNPxl Options");
    registerIntOption_("RNPxl:length", "", 2, "Oligonucleotide maximum length. 0 = disable search for NA variants.", false);

    registerStringOption_("RNPxl:sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence).", false);

    registerStringList_("RNPxl:target_nucleotides", 
                        "", 
                        {"A=C10H14N5O7P", "C=C9H14N3O8P", "G=C10H14N5O8P", "U=C9H13N2O9P"}, 
                        "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG.", 
                        false, 
                        false);

    registerStringList_("RNPxl:nt_groups",
        "",
        {},
	"Restrict which nucleotides can cooccur in a precursor adduct to be able to search both RNA and DNA (format: 'AU CG').",
        false,
        false);

    registerStringList_("RNPxl:mapping", "", {"A->A", "C->C", "G->G", "U->U"}, "format: source->target e.g. A->A, ..., U->U, U->X.", false, false);

    // define if nucleotide can cross-link (produce y,b,a,immonium-ion shifts) in addition to marker ions
    registerStringOption_("RNPxl:can_cross_link", 
                        "<option>", 
                        "U", 
                        "format: 'U' if only U forms cross-links. 'CATG' if C, A, G, and T form cross-links.", 
                        false, 
                        false);

    StringList modifications;
    modifications.emplace_back("U:");
    modifications.emplace_back("U:-H2O");
    modifications.emplace_back("U:-H2O-HPO3");
    modifications.emplace_back("U:-HPO3");

    // fragment adducts that may occur for every precursor adduct (if chemically feasible in terms of elements may not be negative)
    StringList fragment_adducts = {"U:C9H10N2O5;U-H3PO4", 
                                   "U:C4H4N2O2;U'", 
                                   "U:C4H2N2O1;U'-H2O",
                                   "U:C3O;C3O",
                                   "U:C9H13N2O9P1;U",
                                   "U:C9H11N2O8P1;U-H2O",
                                   "U:C9H12N2O6;U-HPO3"
                                  };    

    registerStringList_("RNPxl:fragment_adducts", 
                        "", 
                        fragment_adducts, 
                        "format: [target nucleotide]:[formula] or [precursor adduct]->[fragment adduct formula];[name]: e.g., 'U:C9H10N2O5;U-H3PO4' or 'U:U-H2O->C9H11N2O8P1;U-H2O'.", 
                        false, 
                        false);

    registerStringList_("RNPxl:modifications", "", modifications, "format: empirical formula e.g U:  U:-H2O, ..., U:H2O+PO3.", false, false);

    registerStringOption_("RNPxl:scoring", "<method>", "fast", "Scoring algorithm used in prescoring (fast: total-loss, slow: all losses).", false, false);
    setValidStrings_("RNPxl:scoring", {"fast", "slow"});

    registerStringOption_("RNPxl:decoys", "<bool>", "true", "Generate decoys internally (recommended).", false, false);
    setValidStrings_("RNPxl:decoys", {"true", "false"});

    registerFlag_("RNPxl:CysteineAdduct", "Use this flag if the +152 adduct is expected.", true);
    registerFlag_("RNPxl:filter_fractional_mass", "Use this flag to filter non-crosslinks by fractional mass.", true);
    registerFlag_("RNPxl:carbon_labeled_fragments", "Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).", true);
    registerFlag_("RNPxl:only_xl", "Only search cross-links and ignore non-cross-linked peptides.", true);

    registerDoubleOption_("RNPxl:filter_small_peptide_mass", "<threshold>", 600.0, "Filter precursor that can only correspond to non-crosslinks by mass.", false, true);
    registerDoubleOption_("RNPxl:marker_ions_tolerance", "<tolerance>", 0.03, "Tolerance used to determine marker ions (Da).", false, true);
  
    registerStringList_("filter", "<list>", {"filter_pc_mass_error", "autotune", "idfilter"}, "Filtering steps applied to results.", false, true);
    setValidStrings_("filter", {"filter_pc_mass_error", "impute_decoy_medians", "filter_bad_partial_loss_scores", "autotune", "idfilter", "spectrumclusterfilter", "pcrecalibration"}); 


    registerDoubleOption_("window_size", "<number>", 75.0, "Peak window for spectra precprocessing.", false, true);
    registerIntOption_("peak_count", "<number>", 20, "Retained peaks in peak window.", false, true);
  }



  // bad score or less then two peaks matching and less than 1% explained signal
  static bool badTotalLossScore(float hyperScore, float tlss_Morph, float tlss_modds, float tlss_total_MIC)
  {
    return (hyperScore < MIN_HYPERSCORE 
      || tlss_Morph < MIN_TOTAL_LOSS_IONS + 1.0
      || tlss_total_MIC < 0.01); 
  }


  static bool badPartialLossScore(float tlss_Morph, float plss_Morph, float plss_MIC, float plss_im_MIC, float plss_pc_MIC, float marker_ions_score)
  {
#if !defined DONT_ACCUMULATE_PARTIAL_ION_SCORES
    // if partial loss scores accumulate on the total loss scores, we first need to calculate the individual components
    plss_Morph -= tlss_Morph;
    float tlss_MIC = tlss_Morph - static_cast<int>(tlss_Morph);
    plss_MIC -= tlss_MIC;
#endif

    if (plss_Morph + tlss_Morph < 5.03) return true; // less than 5 peaks? 3% TIC

    if (plss_MIC + plss_im_MIC + plss_pc_MIC + marker_ions_score < 0.03) return true;

    // if we don't see shifted ladder ions, we need at least some signal in the shifted immonium ions
    return (plss_Morph < MIN_SHIFTED_IONS && plss_im_MIC < 0.03);
  }


/*
  // less strict filtering
  static bool badPartialLossScore(float tlss_Morph, float plss_Morph, float plss_MIC, float plss_im_MIC, float plss_pc_MIC, float marker_ions_score)
  {
    if (plss_Morph + tlss_Morph < 3.0) return true; // less than 3 matched peaks?

    // no signal in NA specific peaks?
    if (plss_MIC + plss_im_MIC + plss_pc_MIC + marker_ions_score == 0.0) return true;

    // if we don't see shifted ladder ions, we need at least some signal in the shifted immonium ions
    return false; 
  }
*/

  /*
     @param N number of theoretical peaks
     @param peak_in_spectrum number of experimental peaks
     @param matched_size number of matched theoretical peaks
  */ 
  static double matchOddsScore_(
    const Size& N,
    const float fragment_mass_tolerance_Da,
    const Size peaks_in_spectrum,
    const float mass_range_Da,
    const Size matched_size)
  {    
    if (matched_size < 1 || N < 1) { return 0; }

    // Nd/w (number of peaks in spectrum * fragment mass tolerance in Da / MS/MS mass range in Da - see phoshoRS)
    //const double p = peaks_in_spectrum * fragment_mass_tolerance_Da / mass_range_Da;

    const double p = 20.0 / 100.0; // level 20.0 / mz 100.0 (see WindowMower)
    const double pscore = boost::math::ibeta(matched_size + 1, N - matched_size, p);
    if (pscore <= std::numeric_limits<double>::min()) return -log10(std::numeric_limits<double>::min());
    const double minusLog10p1pscore = -log10(pscore);
    return minusLog10p1pscore;
  } 

  static void generateTheoreticalMZsZ1_(const AASequence& peptide, 
    const Residue::ResidueType& res_type, 
    std::vector<double>& mzs)
  {    
    const Size N = peptide.size();
    mzs.resize(N-1);
    double mono_weight(Constants::PROTON_MASS_U);
    if (res_type == Residue::BIon || res_type == Residue::AIon || res_type == Residue::CIon)
    {
      if (peptide.hasNTerminalModification())
      {
        mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
      }

      switch (res_type)
      {
        case Residue::AIon: mono_weight += Residue::getInternalToAIon().getMonoWeight(); break;
        case Residue::BIon: mono_weight += Residue::getInternalToBIon().getMonoWeight(); break;
        case Residue::CIon: mono_weight += Residue::getInternalToCIon().getMonoWeight(); break;
        default: break;
      }

      for (Size i = 0; i < N-1; ++i) // < N-1: don't add last residue as it is part of precursor
      {
        mono_weight += peptide[i].getMonoWeight(Residue::Internal);
        mzs[i] = mono_weight;
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if (peptide.hasCTerminalModification())
      {
        mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
      }

      switch (res_type)
      {
        case Residue::XIon: mono_weight += Residue::getInternalToXIon().getMonoWeight(); break;
        case Residue::YIon: mono_weight += Residue::getInternalToYIon().getMonoWeight(); break;
        case Residue::ZIon: mono_weight += Residue::getInternalToZIon().getMonoWeight(); break;
        default: break;
      }

      for (Size i = N-1; i > 0; --i) // > 0: don't add last residue (part of precursor ion)
      {
        mono_weight += peptide[i].getMonoWeight(Residue::Internal);
        mzs[N-1-i] = mono_weight;
      } 
    } 
  } 

  static double logfactorial_(UInt x)
  {
    if (x < 2) { return 0; }
    double z(0);
    for (double y = 2; y <= static_cast<double>(x); ++y) { z += log(static_cast<double>(y)); }
      return z;
    }

  // score ions without nucleotide shift
  static void scorePeptideIons_(
      const PeakSpectrum& exp_spectrum,
      const DataArrays::IntegerDataArray& exp_charges,
      const vector<double>& total_loss_template_z1_b_ions,
      const vector<double>& total_loss_template_z1_y_ions,
      const double peptide_mass_without_NA,
      const unsigned int pc_charge,
      const ImmoniumIonsInPeptide& iip, 
      const double fragment_mass_tolerance, 
      const bool fragment_mass_tolerance_unit_ppm,
      std::vector<double>& intensity_sum,
      std::vector<double>& b_ions,
      std::vector<double>& y_ions,
      vector<bool>& peak_matched, 
      float& hyperScore,
      float& MIC,
      float& Morph,
      float& modds,
      float& err,
      float& pc_MIC,
      float& im_MIC,
      size_t& n_theoretical_peaks)
  {
    OPENMS_PRECONDITION(exp_spectrum.size() >= 1, "Experimental spectrum empty.");
    OPENMS_PRECONDITION(exp_charges.size() == exp_spectrum.size(), "Error: HyperScore: #charges != #peaks in experimental spectrum.");
    OPENMS_PRECONDITION(total_loss_template_z1_b_ions.size() == total_loss_template_z1_y_ions.size(), "b- and y-ion arrays must have same size.");
    OPENMS_PRECONDITION(total_loss_template_z1_b_ions.size() > 0, "b- and y-ion arrays must not be empty.");
    OPENMS_PRECONDITION(intensity_sum.size() == total_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(peak_matched.size() == exp_spectrum.size(), "Peak matched needs to be of same size as experimental spectrum");
    OPENMS_PRECONDITION(std::count_if(peak_matched.begin(), peak_matched.end(), [](bool b){return b == true;}) == 0, "Peak matched must be initialized to false");

    double dot_product(0.0), b_mean_err(0.0), y_mean_err(0.0);
    const Size N = intensity_sum.size();

    size_t matches(0);

    // maximum charge considered
    const unsigned int max_z = std::min(2U, static_cast<unsigned int>(pc_charge - 1));

    // match b- and a-ions (we record a-ions as b-ions)
    for (double diff2b : {0.0, -27.994915} ) // b-ion and a-ion ('CO' mass diff from b- to a-ion)
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (Size i = 0; i < total_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (total_loss_template_z1_b_ions[i] + diff2b 
            + (z-1) * Constants::PROTON_MASS_U) / z;
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              const double intensity = exp_spectrum[index].getIntensity();
              dot_product += intensity;
              b_mean_err += intensity * std::abs(theo_mz - exp_mz);
              b_ions[i] += intensity;
              ++matches;
              peak_matched[index] = true;
            }
          }
          ++n_theoretical_peaks;
        }
      }
    }

    // match y-ions
    for (Size z = 1; z <= max_z; ++z)
    {
      for (Size i = 0; i < total_loss_template_z1_y_ions.size(); ++i)
      {
        const double theo_mz = (total_loss_template_z1_y_ions[i] + (z-1) * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size index = exp_spectrum.findNearest(theo_mz);

        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
        {
          const double intensity = exp_spectrum[index].getIntensity();
          y_mean_err += intensity * std::abs(theo_mz - exp_mz);
          dot_product += intensity;                  
          y_ions[N-1 - i] += intensity;      

          if (!peak_matched[index])
          {
            const double intensity = exp_spectrum[index].getIntensity();
            y_mean_err += intensity * std::abs(theo_mz - exp_mz);
            dot_product += intensity;                  
            y_ions[N-1 - i] += intensity;      
            ++matches;
            peak_matched[index] = true;
          }
        }
        ++n_theoretical_peaks;
      }
    }

#ifdef CONSIDER_AA_LOSSES
    // block peaks matching to AA related neutral losses so they get counted for explained peak fraction calculation
    // b-H2O
    for (double diff2b : { -18.010565 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (Size i = 0; i < total_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (total_loss_template_z1_b_ions[i] + diff2b 
            + (z-1) * Constants::PROTON_MASS_U) / z;
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];
          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              peak_matched[index] = true;
            }
          }
        }
      }
    }

    // y-H2O and y-NH3
    for (double diff2b : { -18.010565, -17.026549 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (Size i = 0; i < total_loss_template_z1_y_ions.size(); ++i)
        {
          const double theo_mz = (total_loss_template_z1_y_ions[i] + diff2b
             + (z-1) * Constants::PROTON_MASS_U) / z;
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];
          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              peak_matched[index] = true;
            }
          }
        }
      }
    }
#endif

    // determine b+a and y-ion count 
    UInt y_ion_count(0), b_ion_count(0);
    double b_sum(0.0);
    for (Size i = 0; i != b_ions.size(); ++i) 
    {
      if (b_ions[i] > 0) 
      {
        intensity_sum[i] += b_ions[i];
        b_sum += b_ions[i];
        ++b_ion_count;
      }       
    } 

    double y_sum(0.0);
    for (Size i = 0; i != y_ions.size(); ++i) 
    {
      if (y_ions[i] > 0) 
      {
        intensity_sum[i] += y_ions[i];
        y_sum += y_ions[i];
        ++y_ion_count;
      }       
    }

    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].getName() == "TIC", "No TIC stored in spectrum meta data.");
    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].size() == 1, "Exactly one TIC expected.");

    const double& TIC = exp_spectrum.getFloatDataArrays()[0][0];

    if (y_ion_count == 0 && b_ion_count == 0) 
    {
      hyperScore = 0;
      MIC = 0;
      Morph = 0;
      err = fragment_mass_tolerance_unit_ppm ? 2000.0 * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
    }
    else
    {
      const double yFact = logfactorial_(y_ion_count);
      const double bFact = logfactorial_(b_ion_count);
      hyperScore = log1p(dot_product) + yFact + bFact;
      MIC = std::accumulate(intensity_sum.begin(), intensity_sum.end(), 0.0);
      for (auto& i : intensity_sum) { i /= TIC; } // scale intensity sum
      MIC /= TIC;
      Morph = b_ion_count + y_ion_count + MIC;
      err = (y_mean_err + b_mean_err)/(b_sum + y_sum);
    }

    // match precusor ions z = 1..pc_charge
    double pc_match_count(0);
    for (double pc_loss : {0.0, -18.010565, -17.026548} ) // normal, loss of water, loss of ammonia
    { 
      for (Size z = 1; z <= pc_charge; ++z)
      {
        const double theo_mz = (peptide_mass_without_NA + pc_loss + z * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
        Size index = exp_spectrum.findNearest(theo_mz);
        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
        {
          if (!peak_matched[index])
          {
            const double intensity = exp_spectrum[index].getIntensity();
            pc_MIC += intensity;
            pc_match_count += 1.0;
            ++matches;
            peak_matched[index] = true;
          }
        }
        ++n_theoretical_peaks;
      }      
    }
    pc_MIC /= TIC;
    pc_MIC += pc_match_count;  // Morpheus score 

    // shifted immonium ions

    // lambda to match one peak and sum up in score
    auto match_one_peak_z1 = [&](const double& theo_mz, float& score)
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;      
        auto index = exp_spectrum.findNearest(theo_mz);
        if (exp_charges[index] == 1 && 
          std::abs(theo_mz - exp_spectrum[index].getMZ()) < max_dist_dalton) // found peak match
        {
          if (!peak_matched[index])
          {
            score += exp_spectrum[index].getIntensity();      
            ++matches;
            peak_matched[index] = true;
          }
        } 
        ++n_theoretical_peaks;
      };

    // see DOI: 10.1021/pr3007045 A Systematic Investigation into the Nature of Tryptic HCD Spectra
    static const double imY = EmpiricalFormula("C8H10NO").getMonoWeight(); // 85%
    static const double imW = EmpiricalFormula("C10H11N2").getMonoWeight(); // 84%
    static const double imF = EmpiricalFormula("C8H10N").getMonoWeight(); // 84%
    static const double imL = EmpiricalFormula("C5H12N").getMonoWeight(); // I/L 76%
    static const double imH = EmpiricalFormula("C5H8N3").getMonoWeight(); // 70%
    static const double imC = EmpiricalFormula("C2H6NS").getMonoWeight(); // CaC 61%
    static const double imK1 = EmpiricalFormula("C5H13N2").getMonoWeight(); // 2%
    static const double imP = EmpiricalFormula("C4H8N").getMonoWeight(); //?
    static const double imQ = 101.0715; // 52%
    static const double imE = 102.0555; // 37%
    static const double imM = 104.0534; // 3%
//  static const double imN = 87.05584; // 11%
//  static const double imD = 88.03986; // 4%

    if (iip.Y) 
    {
      match_one_peak_z1(imY, im_MIC);
    }
    if (iip.W) 
    {
      match_one_peak_z1(imW, im_MIC);
    }
    if (iip.F) 
    {
      match_one_peak_z1(imF, im_MIC);
    }
    if (iip.H) 
    {
      match_one_peak_z1(imH, im_MIC);
    }
    if (iip.C) 
    {
      match_one_peak_z1(imC, im_MIC);
    }
    if (iip.P) 
    {
      match_one_peak_z1(imP, im_MIC);
    }
    if (iip.L) 
    {
      match_one_peak_z1(imL, im_MIC);
    }
    if (iip.K) 
    {
      match_one_peak_z1(imK1, im_MIC);
    }
    if (iip.M) 
    {
      match_one_peak_z1(imM, im_MIC);
    }
    if (iip.Q) 
    {
      match_one_peak_z1(imQ, im_MIC);
    }
    if (iip.E) 
    {
      match_one_peak_z1(imE, im_MIC);
    }
    im_MIC /= TIC;

    // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
    err = Morph > 2 ? err : 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;

    const float fragment_mass_tolerance_Da = 2.0 * fragment_mass_tolerance * 1e-6 * 1000;

#ifdef MODDS_ON_ABY_IONS_ONLY
    modds = matchOddsScore_(total_loss_template_z1_b_ions.size() + total_loss_template_z1_y_ions.size(),
     fragment_mass_tolerance_Da,
     exp_spectrum.size(),
     exp_spectrum.back().getMZ(),
     (int)Morph);
#else
    modds = matchOddsScore_(n_theoretical_peaks,
     fragment_mass_tolerance_Da,
     exp_spectrum.size(),
     exp_spectrum.back().getMZ(),
     matches);
#endif
  }

  static void scoreShiftedLadderIons_(
                        const vector<RNPxlFragmentAdductDefinition>& partial_loss_modification,
                        const vector<double>& partial_loss_template_z1_b_ions,
                        const vector<double>& partial_loss_template_z1_y_ions,
                        const double peptide_mass_without_NA,
                        const unsigned int pc_charge,
                        const ImmoniumIonsInPeptide& iip,
                        const double fragment_mass_tolerance, 
                        const bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        std::vector<double>& intensity_sum,
                        std::vector<double>& b_ions,
                        std::vector<double>& y_ions,
                        std::vector<bool>& peak_matched,
                        float& plss_hyperScore,
                        float& plss_MIC,
                        float& plss_Morph,
                        float& plss_err,
                        float& plss_modds,
                        float& plss_pc_MIC,
                        float& plss_im_MIC,
                        size_t& n_theoretical_peaks,
                        bool is_decoy)
  {
    OPENMS_PRECONDITION(exp_spectrum.size() >= 1, "Experimental spectrum empty.");
    OPENMS_PRECONDITION(exp_charges.size() == exp_spectrum.size(), "Error: HyperScore: #charges != #peaks in experimental spectrum.");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(partial_loss_template_z1_b_ions.size() == partial_loss_template_z1_y_ions.size(), "b- and y-ion arrays must have same size.");
    OPENMS_PRECONDITION(partial_loss_template_z1_b_ions.size() > 0, "b- and y-ion arrays must not be empty.");

#ifdef FILTER_AMBIGIOUS_PEAKS
    // is the mass shift ambigious and common in non-XL data (e.g.: freq > 10%)
    auto ambigious_match = [&](const double& previous_theo_mass, const double& current_theo_mass)->bool
      {
        const double dist = current_theo_mass - previous_theo_mass;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? current_theo_mass * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
        auto low_it = mass2high_frequency_.lower_bound(dist - max_dist_dalton);
        if (low_it != mass2high_frequency_.end() && fabs(low_it->first - dist) < max_dist_dalton) // in tolerance window?
        {
          return true; // ambigious match
        }
        return false;
      };
#endif

    double dot_product(0.0), b_mean_err(0.0), y_mean_err(0.0);
    const Size N = intensity_sum.size(); // number of bonds = length of peptide - 1

    size_t matches(0);

    // maximum charge considered
    const unsigned int max_z = std::min(2U, static_cast<unsigned int>(pc_charge - 1));

    // match b- and a-ions (we record a-ions as b-ions)
    for (double diff2b : {0.0, -27.994915} ) // b-ion and a-ion ('CO' mass diff from b- to a-ion)
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (const RNPxlFragmentAdductDefinition & fa : partial_loss_modification)
        {
          for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)

          {
            const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass + diff2b 
              + (z-1) * Constants::PROTON_MASS_U) / z;

            const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

            // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
            Size index = exp_spectrum.findNearest(theo_mz);

            const double exp_mz = exp_spectrum[index].getMZ();
            const Size exp_z = exp_charges[index];

            // found peak match
            if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
            {
              if (!peak_matched[index])
              {
#ifdef FILTER_AMBIGIOUS_PEAKS
                // skip ambigious matches
                if (i >= 1)
                {
                  if (ambigious_match(partial_loss_template_z1_b_ions[i - 1], partial_loss_template_z1_b_ions[i] + fa.mass)) continue;
                }
#endif
                const double intensity = exp_spectrum[index].getIntensity();
                b_mean_err += intensity * std::abs(theo_mz - exp_mz);
                dot_product += intensity;
                b_ions[i] += intensity;            
                peak_matched[index] = true;
                matches++;
              }
            }
            ++n_theoretical_peaks;
          }
        } 
      }
    } 

    // match y-ions
    for (Size z = 1; z <= max_z; ++z)
    {
      for (const RNPxlFragmentAdductDefinition  & fa : partial_loss_modification)
      {
        for (Size i = 1; i < partial_loss_template_z1_y_ions.size(); ++i)  // Note that we start at (i=1 -> y2) as trypsin would otherwise not cut at cross-linking site
        {
          const double theo_mz = (partial_loss_template_z1_y_ions[i] + fa.mass 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
#ifdef FILTER_AMBIGIOUS_PEAKS
              // skip ambigious matches
              if (i >= 1)
              {
                  if (ambigious_match(partial_loss_template_z1_y_ions[i - 1], partial_loss_template_z1_y_ions[i] + fa.mass)) continue;
              }
#endif
              const double intensity = exp_spectrum[index].getIntensity();
              y_mean_err += intensity * std::abs(theo_mz - exp_mz);
              dot_product += intensity;                  
              y_ions[N-1 - i] += intensity;      
              peak_matched[index] = true;
              matches++;
            }
          }
          ++n_theoretical_peaks;
        }
      }  
    }


#ifdef CONSIDER_AA_LOSSES
    // block peaks matching to AA related neutral losses so they don't get matched to NA shifts
    // b-H2O
    for (double diff2b : { -18.010565 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (const RNPxlFragmentAdductDefinition & fa : partial_loss_modification)
        {
          for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
          {
            const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass + diff2b 
              + (z-1) * Constants::PROTON_MASS_U) / z;

            const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
            Size index = exp_spectrum.findNearest(theo_mz);
            const double exp_mz = exp_spectrum[index].getMZ();
            const Size exp_z = exp_charges[index];
            if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
            {
              if (!peak_matched[index])
              {
                peak_matched[index] = true;
              }
            }
          }
        } 
      }
    } 

  // match y-ions
  // y-H2O and y-NH3
  for (double diff2b : { -18.010565, -17.026549 } ) 
    for (Size z = 1; z <= max_z; ++z)
    {
      for (const RNPxlFragmentAdductDefinition  & fa : partial_loss_modification)
      {
        for (Size i = 1; i < partial_loss_template_z1_y_ions.size(); ++i)  // Note that we start at (i=1 -> y2) as trypsin would otherwise not cut at cross-linking site
        {
          const double theo_mz = (partial_loss_template_z1_y_ions[i] + fa.mass + diff2b
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              peak_matched[index] = true;
            }
          }
        }
      }  
    }
#endif

    UInt y_ion_count(0), b_ion_count(0);
    double b_sum(0.0);
    for (Size i = 0; i != b_ions.size(); ++i) 
    {
      if (b_ions[i] > 0) 
      {
        intensity_sum[i] += b_ions[i];
        b_sum += b_ions[i];
        ++b_ion_count;
      }       
    } 

    double y_sum(0.0);
    for (Size i = 0; i != y_ions.size(); ++i) 
    {
      if (y_ions[i] > 0) 
      {
        intensity_sum[i] += y_ions[i];
        y_sum += y_ions[i];
        ++y_ion_count;
      }       
    }

    const double& TIC = exp_spectrum.getFloatDataArrays()[0][0];

    if (y_ion_count == 0 && b_ion_count == 0) 
    {
      plss_hyperScore = 0;
      plss_MIC = 0;
      plss_Morph = 0;
      plss_err = fragment_mass_tolerance_unit_ppm ? 2000.0 * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
    }
    else
    {
      const double yFact = logfactorial_(y_ion_count);
      const double bFact = logfactorial_(b_ion_count);
      plss_hyperScore = log1p(dot_product) + yFact + bFact;
      plss_MIC = std::accumulate(intensity_sum.begin(), intensity_sum.end(), 0.0);
      for (auto& i : intensity_sum) { i /= TIC; } // scale intensity sum
      plss_MIC /= TIC;
      plss_Morph = b_ion_count + y_ion_count + plss_MIC;
      plss_err = (y_mean_err + b_mean_err)/(b_sum + y_sum);
   }
    
    // match (partially) shifted precusor ions z = 1..pc_charge
    double pc_match_count(0);
    for (double pc_loss : {0.0, -18.010565, -17.026548} ) // normal, loss of water, loss of ammonia
    { 
      const double peptide_mass = peptide_mass_without_NA + pc_loss;
      for (Size z = 1; z <= pc_charge; ++z)
      {
        for (const RNPxlFragmentAdductDefinition & fa : partial_loss_modification)
        {
          const double theo_mz = (peptide_mass + fa.mass + z * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              const double intensity = exp_spectrum[index].getIntensity();
              plss_pc_MIC += intensity;
              pc_match_count += 1.0;
              peak_matched[index] = true;
              matches++;
            }
          }
          ++n_theoretical_peaks;
        }      
      }
    }
    plss_pc_MIC /= TIC;
    plss_pc_MIC += pc_match_count; // Morpheus score

    ////////////////////////////////////////////////////////////////////////////////
    // match shifted immonium ions

    // lambda to match one peak and sum up in score
    auto match_one_peak_z1 = [&](const double& theo_mz, float& score)
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;      
        auto index = exp_spectrum.findNearest(theo_mz);
        if (exp_charges[index] == 1 && 
          std::abs(theo_mz - exp_spectrum[index].getMZ()) < max_dist_dalton) // found peak match
        {
          if (!peak_matched[index])
          {
            score += exp_spectrum[index].getIntensity();      
            peak_matched[index] = true;
            matches++;
          }
        } 
        ++n_theoretical_peaks;
      };

    static const double imY = EmpiricalFormula("C8H10NO").getMonoWeight();
    static const double imW = EmpiricalFormula("C10H11N2").getMonoWeight();
    static const double imF = EmpiricalFormula("C8H10N").getMonoWeight();
    static const double imH = EmpiricalFormula("C5H8N3").getMonoWeight();
    static const double imC = EmpiricalFormula("C2H6NS").getMonoWeight();
    static const double imP = EmpiricalFormula("C4H8N").getMonoWeight();
    static const double imL = EmpiricalFormula("C5H12N").getMonoWeight();
    static const double imK1 = EmpiricalFormula("C5H13N2").getMonoWeight();
    static const double imK2 = EmpiricalFormula("C5H10N1").getMonoWeight();
    static const double imK3 = EmpiricalFormula("C6H13N2O").getMonoWeight();
    static const double imQ = 101.0715;
    static const double imE = 102.0555;
    static const double imM = 104.0534;

    for (const RNPxlFragmentAdductDefinition & fa : partial_loss_modification)
    {
      if (iip.Y) 
      {
        match_one_peak_z1(imY + fa.mass, plss_im_MIC);
      }
      if (iip.W) 
      {
        match_one_peak_z1(imW + fa.mass, plss_im_MIC);
      }
      if (iip.F) 
      {
        match_one_peak_z1(imF + fa.mass, plss_im_MIC);
      }
      if (iip.H) 
      {
        match_one_peak_z1(imH + fa.mass, plss_im_MIC);
      }
      if (iip.C) 
      {
        match_one_peak_z1(imC + fa.mass, plss_im_MIC);
      }
      if (iip.P) 
      {
        match_one_peak_z1(imP + fa.mass, plss_im_MIC);
      }
      if (iip.L) 
      {
        match_one_peak_z1(imL + fa.mass, plss_im_MIC);
      }
      if (iip.K) 
      {
        match_one_peak_z1(imK1 + fa.mass, plss_im_MIC);
        // according to A. Stuetzer mainly observed with C‘-NH3 (94.0167 Da)
        match_one_peak_z1(imK2 + fa.mass, plss_im_MIC);
        // usually only observed without shift (A. Stuetzer)
        // TODO: only enable for DNA? get's sometimes matched by chance  
        // match_one_peak_z1(imK3 + fa.mass, plss_im_MIC); 
      }
      if (iip.M) 
      {
        match_one_peak_z1(imM + fa.mass, plss_im_MIC);
      }
      if (iip.Q) 
      {
        match_one_peak_z1(imQ + fa.mass, plss_im_MIC);
      }
      if (iip.E) 
      {
        match_one_peak_z1(imE + fa.mass, plss_im_MIC);
      }
    }
    plss_im_MIC /= TIC;

    // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
    plss_err = plss_Morph > 2 ? plss_err : 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;

    const float fragment_mass_tolerance_Da = 2.0 * fragment_mass_tolerance * 1e-6 * 1000;

#ifdef MODDS_ON_ABY_IONS_ONLY
    plss_modds = matchOddsScore_(
     partial_loss_template_z1_b_ions.size() + partial_loss_template_z1_y_ions.size(), 
     fragment_mass_tolerance_Da,
     exp_spectrum.size(),
     exp_spectrum.back().getMZ(),
     (int)plss_Morph);
#else
    plss_modds = matchOddsScore_(
     n_theoretical_peaks, 
     fragment_mass_tolerance_Da,
     exp_spectrum.size(),
     exp_spectrum.back().getMZ(),
     matches);
#endif

  } 

/*
*  Combine subscores of all-ion scoring.
*/
  static float calculateCombinedScore(const AnnotatedHit& ah, 
    const bool isXL, 
    const double nucleotide_mass_tags
    //, const double fraction_of_top50annotated
    )
  {
// Tie-braker score
//    return + 1.0 * ah.total_loss_score + ah.total_MIC + 0.1 * ah.mass_error_p 
//           - 0.01 * ah.isotope_error - 10.0 * ah.err - 10.0 * ah.pl_err;

/*
     double score = 2.52872532
                   +0.38318303 * nucleotide_mass_tags;
           
    if (!isXL)
    {
       score += ah.mass_error_p * 0.45919677
+ ah.err          * 0.01016288
+ ah.modds        * -0.02450589
+ ah.immonium_score  * 0.26555840
+ ah.precursor_score * 0.06148951
+ ah.MIC             * 0.91845925
+ ah.sequence_score  * 0.23213255
+ ah.total_loss_score * 8.69880275;
    }
    else
    {
score += ah.mass_error_p     *   1.15386068
+ah.err         *       -0.75849696
+ah.pl_err       *       0.01731052
+ah.marker_ions_score *  0.40870416
+ah.total_loss_score *   4.92806210
+ah.modds           *    0.96869679
+ah.immonium_score   *   0.14292426
+ah.precursor_score *   -0.05822564
+ah.MIC            *     0.73432514
+ah.pl_MIC         *     0.27953670
+ah.pl_modds       *     0.03810840
+ah.pl_pc_MIC      *     0.15083043
+ah.pl_im_MIC      *    -0.12681649
+ah.sequence_score  *    0.46228471;
    }

    return score;
*/

//TODO: check Alex    if (ah.Morph + ah.pl_Morph < 5.03) return 0;
//    return ah.total_MIC + ah.marker_ions_score + fraction_of_top50annotated;


/*
    return 
             10.0 * ah.total_loss_score + ah.partial_loss_score 
           + 0.01 * ah.mass_error_p 
           - 10.0 * ah.err 
           - 10.0 * ah.pl_err
           + 3.0 * ah.pl_MIC
           + isXL * 3.0 * ah.marker_ions_score
           + 3.0 * ah.total_MIC + ah.ladder_score;
*/


/*
    return 
              47.29 * ah.total_loss_score
           +  16.93 * ah.partial_loss_score
           +  69.35 * ah.mass_error_p 
           -  72.84 * ah.err 
           -  75.65 * ah.pl_err
           +  26.34 * ah.total_MIC
           +  21.85 * (ah.pl_pc_MIC > 0 || ah.precursor_score > 0)
           +  39.64 * ah.ladder_score
           +  60.15 * (ah.tag_shifted >= 1);
*/
/*
    return 
              64.4 * ah.total_loss_score
           +  21.3 * ah.partial_loss_score
           +  98.97 * ah.mass_error_p 
           -  57.57 * ah.err 
           -  19.23 * ah.pl_err
           +  0.23 * ah.total_MIC
           +  16.47 * (ah.pl_pc_MIC > 0 || ah.precursor_score > 0)
           +  20.26 * ah.ladder_score
           +  4.25 * (ah.tag_shifted >= 1);
*/
/*
    return 
              0.90 * ah.total_loss_score
           +  0.12 * ah.partial_loss_score
           +  2.19 * ah.mass_error_p 
           -  54.60 * ah.err 
           -  71.06 * ah.pl_err
           +  67.47 * ah.total_MIC
           +  66.27 * (ah.pl_pc_MIC > 0 || ah.precursor_score > 0)
           +  11.45 * ah.ladder_score
           +  28.09 * (ah.tag_shifted >= 1);
*/
           const double ret = 1.0 * ah.total_loss_score
           +  0.5 * ah.partial_loss_score
           +  1.0 * ah.mass_error_p 
           -  0.1 * ah.err 
           -  0.1 * ah.pl_err
           + isXL * 0.1 * ah.marker_ions_score
           +  1.0 * ah.total_MIC
           + isXL * (ah.pl_pc_MIC + ah.precursor_score)
//           +  1.0 * (ah.pl_pc_MIC > 0 || ah.precursor_score > 0)
           +  0.1 * ah.ladder_score
           +  1.0 * (ah.tag_shifted >= 1);
           return ret;


/*
 * [96.45142397 49.37352197  2.61098131 -1.21298196 -4.09934265 28.64381592
 *  53.84125899 41.37187898  7.66281392]
 */


/*
 * good on IN VITRO:
              1.0 * ah.total_loss_score
           +  0.5 * ah.partial_loss_score
           +  1.0 * ah.mass_error_p 
           -  0.1 * ah.err 
           -  0.1 * ah.pl_err
           + isXL * 0.1 * ah.marker_ions_score
           +  1.0 * ah.total_MIC
           +  1.0 * (ah.pl_pc_MIC > 0 || ah.precursor_score > 0)
           +  0.1 * ah.ladder_score
           +  1.0 * (ah.tag_shifted >= 1);
*/
/*
            -10.9457 + 1.1836 * isXL + 1.6076 * ah.mass_error_p - 579.912 * ah.err
           + 52.2888 * ah.pl_err - 0.0105 * ah.modds
           + 88.7997 * ah.immonium_score- 0.88823 * ah.partial_loss_score + 14.2052 * ah.pl_MIC
           + 0.61144 * ah.pl_modds + 10.07574543 * ah.pl_pc_MIC -28.05701 * ah.pl_im_MIC
           + 2.59655 * ah.total_MIC + 2.38320 * ah.ladder_score + 0.65422535 * (ah.total_loss_score + ah.partial_loss_score)
           4267 without perc / 5306 with perc auf 1-
           5010 with perc auf 2-
*/

/*
    return 
       (ah.MIC - 0.00693372) * 0.00317232 
    + (ah.Morph - 2.00693) * -5.67908 
    + (ah.err - 1.8356e-05) * 0.0074405
    + (ah.immonium_score) *   3.87159
    + (isXL) *  -2.8901
 +   (ah.ladder_score - 0.0381864) * 2.09632
 +   (ah.marker_ions_score) *  2.53656
 +   (ah.mass_error_p - 0.6439) * 0.189849
 +   (ah.modds - 3.37613e-07) * 9.42353
 +   (ah.partial_loss_score) * 4.22367
 +   (ah.pl_MIC) * 0.0770935
 +   (ah.pl_Morph) * -4.25703
 +   (ah.pl_err - 0.000214374) * -0.059863
 +   (ah.pl_im_MIC) *  -0.367491
 +   (ah.pl_modds) *  2.29632
 +   (ah.pl_pc_MIC) *  1.13888
 +   (ah.precursor_score) * -0.597915
 +   (ah.rank_product - 1) * 0.460718
 +   (ah.sequence_score - 0.0218818) * -4.68762
 +   (ah.tag_XLed) * 2.61103
 +   (ah.tag_shifted) * -0.10144
 +   (ah.tag_unshifted - 1) * 10.8324
 +   (ah.total_MIC - 0.0100453) * 0.836277
 +   (ah.total_loss_score - 0.106861) * 34.8899
 +   (ah.wTop50 - 2) * -0.248061
 +   (ah.sequence.size() - 6) * 0.00656285; // length
*/
  }


  static float calculateFastScore(const AnnotatedHit& ah)
  {
    return + 1.0 * ah.total_loss_score
/*               + 1.0 * ah.total_MIC         
               + 0.333 * ah.mass_error_p*/;
  } 
/*
*  Score fragments carrying NA adducts 
*/
static void scoreXLIons_(
                         const vector<RNPxlFragmentAdductDefinition> &partial_loss_modification,
                         const ImmoniumIonsInPeptide& iip,
                         const PeakSpectrum &exp_spectrum,
                         const double peptide_mass_without_NA,
                         double fragment_mass_tolerance,
                         bool fragment_mass_tolerance_unit_ppm,
                         const vector<double> &partial_loss_template_z1_b_ions,
                         const vector<double> &partial_loss_template_z1_y_ions,
                         const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                         vector<double>& intensity_sum,
                         vector<double>& b_ions,
                         vector<double>& y_ions,
                         vector<bool>& matched_peaks,
                         float &partial_loss_sub_score,
                         float &marker_ions_sub_score,
                         float &plss_MIC, 
                         float &plss_err, 
                         float &plss_Morph,
                         float &plss_modds,
                         float &plss_pc_MIC,
                         float &plss_im_MIC,
                         size_t &n_theoretical_peaks,
                         bool is_decoy)
  {
    OPENMS_PRECONDITION(!partial_loss_template_z1_b_ions.empty(), "Empty partial loss spectrum provided.");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_y_ions.size(), "Sum array needs to be of same size as y-ion array");

    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
    //const double exp_pc_mz = exp_spectrum.getPrecursors()[0].getMZ();

    if (!marker_ions_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance * 2.0,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                                             marker_ions_sub_score_spectrum_z1,
                                             marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[IA_CHARGE_INDEX]);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;

      // count marker ions
      n_theoretical_peaks += marker_ions_sub_score_spectrum_z1.size();
    }

    scoreShiftedLadderIons_(
                      partial_loss_modification,                        
                      partial_loss_template_z1_b_ions,
                      partial_loss_template_z1_y_ions,
                      peptide_mass_without_NA,
                      exp_pc_charge,
                      iip,
                      fragment_mass_tolerance, 
                      fragment_mass_tolerance_unit_ppm, 
                      exp_spectrum, 
                      exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                      intensity_sum,
                      b_ions,
                      y_ions,
                      matched_peaks,
                      partial_loss_sub_score,
                      plss_MIC,
                      plss_Morph,
                      plss_err,
                      plss_modds,
                      plss_pc_MIC,
                      plss_im_MIC,
                      n_theoretical_peaks,
                      is_decoy);
#ifdef DEBUG_OpenNuXL
    LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
    // cap plss_err to something larger than the mean_mz * max_ppm_error
    float ft_da = fragment_mass_tolerance_unit_ppm ? fragment_mass_tolerance * 1e-6 * 1000.0 : fragment_mass_tolerance;
    if (plss_err > ft_da) plss_err = ft_da;
  }

  // De novo tagger
  class OpenNuXLTagger
  {
    public:

    // initalize tagger with minimum/maximum tag length and +- tolerance ppm
  explicit OpenNuXLTagger(float tol = 0.05, size_t min_tag_length = 0, size_t max_tag_length = 65535)
  {
    tol_ = tol;
    min_tag_length_ = min_tag_length;
    max_tag_length_ = max_tag_length;
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural19WithoutI");
    for (const auto& r : aas)
    {
      const char letter = r->getOneLetterCode()[0]; 
      const float mass = r->getMonoWeight(Residue::Internal);
      mass2aa[mass] = letter;
#ifdef DEBUG_OPENNUXL_TAGGER
      cout << "Mass: " << mass << "\t" << letter << endl; 
#endif
    }

    min_gap_ = mass2aa.begin()->first - tol;
    max_gap_ = mass2aa.rbegin()->first + tol;

#ifdef DEBUG_OPENNUXL_TAGGER
    TheoreticalSpectrumGenerator test;
    auto param = test.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "false");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    test.setParameters(param);
    MSSpectrum test_s;
    test.getSpectrum(test_s, AASequence::fromString("TESTPEPTIDE"), 1, 1); 
    cout << "should be ESTPEPTIDE:" << getLongestTag(test_s) << endl; 
#endif
  }

    void getTag(const std::vector<float>& mzs, std::set<std::string>& tags) const 
    {
      // start peak
      if (min_tag_length_ > mzs.size()) return; // avoid segfault
    
      std::string tag;
      for (size_t i = 0; i < mzs.size() - min_tag_length_; ++i)
      {
        getTag_(tag, mzs, i, tags);
        tag.clear();
      }
    }

    // generate tags from mass vector @mzs using the standard residues in ResidueDB
    void getTag(const MSSpectrum& spec, std::set<std::string>& tags) const
    {
      const size_t N = spec.size();
      if (N < min_tag_length_) { return; }
      // copy to float vector (speed)
      std::vector<float> mzs;
      mzs.reserve(N);
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      getTag(mzs, tags); 
    }

    // generate tags from mass vector @mzs using the standard residues in ResidueDB
    std::string getLongestTag(const MSSpectrum& spec) const
    {
      std::set<std::string> tags;
      getTag(spec, tags);
      if (tags.empty()) return "";
      //cout << "Ntags:" << tags.size() << endl;
      //for (const auto & s: tags) { cout << s << endl; }
      const auto longest = std::max_element(tags.cbegin(), tags.cend(),
        [](const std::string& lhs, const std::string& rhs) { return lhs.size() < rhs.size(); });
      //cout << "longest:" << *longest << endl;
      return *longest;
    }

    // note: this is much more efficient
    size_t getLongestTagLength(const MSSpectrum& spec) const
    {
      // simple DP to detect longest tag
      const size_t N = spec.size();
      if (N < 2) return 0;
      std::vector<float> mzs;
      mzs.reserve(N);
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      std::vector<size_t> max_tag(N, 0); // maximum tag length up to this peak
      size_t longest_tag = 0;
      for (size_t i = 0; i < N - 1; ++i)
      {
        for (size_t k = i + 1; k < N; ++k)
        {
          const double gap = mzs[k] - mzs[i];
          if (gap > max_gap_) { break; }
          const char aa = getAAByMass_(gap);
          if (aa == ' ') { continue; } // can't extend tag to k-th peak
          if (max_tag[k] < max_tag[i] + 1) // check if we found a longer tag to k-th peak
          {
            ++max_tag[k]; // update longest tag to this peak
            if (longest_tag < max_tag[k]) { longest_tag = max_tag[k]; }
          }
        }
      }
      return longest_tag;
    }

    private:
      float min_gap_; // will be set to smallest residue mass in ResidueDB
      float max_gap_; // will be set to highest residue mass in ResidueDB
      float tol_; // < tolerance
      size_t min_tag_length_; // < minimum tag length
      size_t max_tag_length_; // < maximum tag length
      std::map<float, char> mass2aa;

      char getAAByMass_(float m) const
      {
        // fast check for border cases
        if (m < min_gap_ || m > max_gap_) return ' ';
        auto left = mass2aa.lower_bound(m - tol_);
        //if (left == mass2aa.end()) return ' '; // cannot happen, since we checked boundaries above
        if (fabs(left->first - m) < tol_) return left->second; 
        return ' ';
      }        

    void getTag_(std::string & tag, const std::vector<float>& mzs, const size_t i, std::set<std::string>& tags) const
    {
      const size_t N = mzs.size();
      size_t j = i + 1;
      // recurse for all peaks in distance < max_gap
      while (j < N) 
      {
        if (tag.size() == max_tag_length_) { return; } // maximum tag size reached? - continue with next parent

        const float gap = mzs[j] - mzs[i];
        if (gap > max_gap_) { return; } // already too far away - continue with next parent
        const char aa = getAAByMass_(gap);
#ifdef DEBUG_OPENNUXL_TAGGER
        cout << i << "\t" << j << "\t" << mzs[i] << "\t" << mzs[j] << "\t" << gap << "\t'" << aa << "'" << endl;
#endif
      
        if (aa == ' ') { ++j; continue; } // can't extend tag
        tag += aa;
        getTag_(tag, mzs, j, tags);
        if (tag.size() >= min_tag_length_) tags.insert(tag);
        tag.pop_back();  // remove last char
        ++j;
      }         
    }
  };

  struct RankScores
  {
    double explained_peak_fraction = 0;
    double wTop50 = 0;
  };

  class SmallestElements
  {
  private:
    int max_size_;
  public:
    priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> pq;
    explicit SmallestElements(size_t size): 
      max_size_(size)
    {
    }

    void tryAdd(size_t v)
    {
       if (pq.size() < max_size_)
       {
         pq.push(v);
         return;
       }
       if (v < pq.top())
       {
         pq.pop(); //get rid of the root
         pq.push(v); //priority queue will automatically restructure
       }
    }
};

  RankScores rankScores_(const MSSpectrum& spectrum, vector<bool> peak_matched)
  {
    RankScores r;
    if (spectrum.empty()) return r;
    const double matched = std::accumulate(peak_matched.begin(), peak_matched.end(), 0);
    vector<double> matched_ranks;
    for (size_t i = 0; i != peak_matched.size(); ++i)
    {
      if (!peak_matched[i]) { continue; }
      matched_ranks.push_back(spectrum.getIntegerDataArrays()[IA_RANK_INDEX][i]);
    }
    std::sort(matched_ranks.begin(), matched_ranks.end());

    // optimal ranking would be 0,1, ... , (number_of_matched_peaks - 1)
    // calculate number of "insertions" of higher-intensity peaks compared to optimal ranking
    size_t sum_rank_diff(matched_ranks[0]);
    for (size_t i = 1; i != matched_ranks.size(); ++i)
    {
      sum_rank_diff += matched_ranks[i] - matched_ranks[0] - 1;
    }

    r.wTop50 = log10(1.0 + sum_rank_diff);
    r.explained_peak_fraction = matched / (double)spectrum.size();
    return r;
  }
/*
  RankScores rankScores_(const MSSpectrum& spectrum, vector<bool> peak_matched)
  {
    double matched = std::accumulate(peak_matched.begin(), peak_matched.end(), 0);
    SmallestElements top7ranks(7);
    RankScores r;
    for (size_t i = 0; i != peak_matched.size(); ++i)
    {
      if (!peak_matched[i]) 
      {
        continue;
      }
      else
      { 
        const double rank = 1 + spectrum.getIntegerDataArrays()[IA_RANK_INDEX][i]; // ranks start at 0 -> add 1
        r.rp += 1.0/matched * log((double)rank);
        top7ranks.tryAdd(rank);         
      }
    }

    size_t median = peak_matched.size() / 2; // init to number of peaks / 2
    for (size_t i = 1; i <= 4;  ++i)
    {
      if (top7ranks.pq.empty()) break;
      median = top7ranks.pq.top();
      top7ranks.pq.pop();
    }
    r.wTop50 = median;
    r.rp = exp(r.rp - 1.0 / matched * lgamma(matched+1)); // = rp / lowest possible rp given number of matches
    return r;
  }
*/

#ifdef FILTER_AMBIGIOUS_PEAKS
  static map<double, double> mass2high_frequency_;
#endif

  void calculateNucleotideTags_(PeakMap& exp, 
    const double fragment_mass_tolerance, 
    const bool fragment_mass_tolerance_unit_ppm,
    const RNPxlParameterParsing::NucleotideToFragmentAdductMap &  nucleotide_to_fragment_adducts)
  {
    // set of all possibly observable fragment adduct masses
    set<double> adduct_mass;
    for (const auto & p : nucleotide_to_fragment_adducts)
    {
      for (const auto & fa : p.second)
      {
        adduct_mass.insert(fa.mass);
      }
    }

    // mass shift to residue + adduct (including no adduct, see below)
    map<double, map<const Residue*, double> > aa_plus_adduct_mass;
    auto residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

    for (const double d : adduct_mass)
    {
      for (const Residue* r : residues)
      {
        double m = d + r->getMonoWeight(Residue::Internal);
        aa_plus_adduct_mass[m][r] = d; // mass, residue, shift mass
      }
    }
    // add mass shits of plain residues
    for (const Residue* r : residues)
    {
      double m = r->getMonoWeight(Residue::Internal);
      aa_plus_adduct_mass[m][r] = 0;
    }

    if (debug_level_ > 0)
    {
      // output ambigious masses
      ofstream of;
      of.open(getStringOption_("in") + ".ambigious_masses.csv");
      of << "Ambigious residues (+adduct) masses that exactly match to other masses." << endl;
      of << "Total\tResidue\tAdduct" << endl;
      for (auto& m : aa_plus_adduct_mass)
      {
        double mass = m.first;
        if (m.second.size() == 1) continue; 
        // more than one residue / adduct registered for that mass
        for (auto& a : m.second)
        {
          of << mass << "\t" << a.first->getOneLetterCode() << "\t" << a.second << "\n";
        }
      }
      of.close(); 
    }

    map<double, size_t> adduct_mass_count;
    map<double, size_t> aa_plus_adduct_mass_count;

    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;
      // faster
      vector<double> mzs;
      vector<double> charges;
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      for (auto const& p : spec.getIntegerDataArrays()[IA_CHARGE_INDEX]) { charges.push_back(p); }
 
      size_t match(0);
      size_t in_mass_range(0);

      // for all peak pairs
      for (Size i = 0; i != mzs.size(); ++i)
      {
        for (Size j = i+1; j < mzs.size(); ++j)
        {
          if (charges[i] != charges[j]) continue;

          double m = mzs[j];
          double dm = m - mzs[i];

          const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, m) : fragment_mass_tolerance;

          if (dm * charges[i] > *adduct_mass.rbegin() + tolerance) break;

          auto left = adduct_mass.lower_bound((dm * charges[i]) - tolerance);
          if (left == adduct_mass.end()) continue;
          ++in_mass_range;
          // count if distance matches to adduct mass
          if (fabs(*left - (dm * charges[i])) < tolerance )
          {
            ++match;
            ++adduct_mass_count[*left];
          }
        } 
      } 

      // count how often a shift matches a residue + adduct mass (including mass 0 for unmodified residue)
      size_t aa_plus_adduct_in_mass_range(0);
      for (Size i = 0; i != mzs.size(); ++i)
      {
        for (Size j = i+1; j < mzs.size(); ++j)
        {
          double m =  mzs[j];
          double dm = m - mzs[i];

          if (charges[i] != charges[j]) continue;

          const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, m) : fragment_mass_tolerance;

          // find all (shifted/normal) residues that match to an observed shift
          auto left = aa_plus_adduct_mass.lower_bound((dm * charges[i]) - tolerance);
          auto right = aa_plus_adduct_mass.upper_bound((dm * charges[i]) + tolerance);
          for (; left != right; ++left)
          {
            ++aa_plus_adduct_in_mass_range;
            if (fabs(left->first - (dm * charges[i])) < tolerance )
            {
              ++aa_plus_adduct_mass_count[left->first];
            }
          }
        } 
      } 

      spec.getFloatDataArrays().resize(2);

      spec.getFloatDataArrays()[1].resize(1);
      spec.getFloatDataArrays()[1][0] = (double)match / (double)in_mass_range;
      spec.getFloatDataArrays()[1].setName("nucleotide_mass_tags");


    }

    // calculate ranks
    cout << "Calculating ranks..." << endl;
    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;

      // initialize original index locations
      vector<size_t> idx(spec.size());
      std::iota(idx.begin(), idx.end(), 0);
        
      // sort indexes based on comparing intensity values (0 = highest intensity)
      sort(idx.begin(), idx.end(),
        [&spec](size_t i1, size_t i2) { return spec[i1].getIntensity() > spec[i2].getIntensity(); });

      spec.getIntegerDataArrays().resize(IA_RANK_INDEX + 1);
      spec.getIntegerDataArrays()[IA_RANK_INDEX].clear();
      for (int rank : idx) { spec.getIntegerDataArrays()[IA_RANK_INDEX].push_back(rank); }
      spec.getIntegerDataArrays()[IA_RANK_INDEX].setName("intensity_rank");
    }
    cout << " done!" << endl;

    cout << "Calculating longest mass tags..." << endl;
    OpenNuXLTagger tagger(0.03, 3);
    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;
      spec.getIntegerDataArrays().resize(IA_DENOVO_TAG_INDEX + 1);
      spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX].resize(1);
      spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] = 0;
      #ifdef CALCULATE_LONGEST_TAG
      size_t longest_tag = tagger.getLongestTagLength(spec);
      spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] = longest_tag; 
      //spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] = tagger.getLongestTag(spec).size(); // slow
      #endif
      spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX].setName("longest_tag");
    }
    cout << " done!" << endl;


    if (debug_level_ > 0) 
    {
      cout << "Distinct residue + adduct masses (including residues without shift): " << aa_plus_adduct_mass_count.size() << endl; 
      // Calculate background statistics on shifts
      cout << "mass\tresidue\tshift:" << endl;
      for (const auto& mra : aa_plus_adduct_mass)
      {
        double m = mra.first;
        const map<const Residue*, double>& residue2adduct = mra.second;
        for (auto& r2a : residue2adduct)
        {
          cout << m << "\t" << r2a.first->getOneLetterCode() << "\t" << r2a.second << endl;
        }
      }
    }

    // reformat to get: amino acid, mass, count statistics
    map<const Residue*, map<double, size_t> > aa2mass2count;
    for (const auto& mc : aa_plus_adduct_mass_count)
    {
      double mass = mc.first;
      size_t count = mc.second;

      auto it = aa_plus_adduct_mass.lower_bound(mass - 1e-6); // "exact" match
      if (it == aa_plus_adduct_mass.end()) continue;

      const map<const Residue*, double>& residue2adduct = it->second;
      for (auto& r2a : residue2adduct)
      {
        const Residue* residue = r2a.first;
        String name = residue->getName();
        aa2mass2count[residue][mass] = count;
      }
    }

    if (debug_level_ > 0) { cout << "Total counts per residue:" << endl; }

#ifdef FILTER_AMBIGIOUS_PEAKS
    mass2high_frequency_.clear();
#endif
    for (const auto& aa2 : aa2mass2count)
    {
      auto& mass2count = aa2.second;
      for (const auto& m2c : mass2count)
      {
        double current_mass = m2c.first;
        size_t current_residue_count = m2c.second;
#ifdef FILTER_AMBIGIOUS_PEAKS
        size_t unmodified_residue_count = mass2count.begin()->second;
#endif
        if (debug_level_ > 0)
        {
          cout << aa2.first->getName() << "\t" << current_mass << "\t" << current_residue_count << endl; // aa, mass, count  
        }

#ifdef FILTER_AMBIGIOUS_PEAKS
        double frequency_normalized = (double)current_residue_count / unmodified_residue_count;  // frequency relative to unmodified residue
        if (frequency_normalized > 0.5 && frequency_normalized != 1.0) // residue with shift as frequent as unmodified residue
        {
          mass2high_frequency_[m2c.first] = frequency_normalized;
        }
#endif
      }
    }

    if (debug_level_ > 0)
    {
      cout << "Normalized counts per residue:" << endl;
      for (const auto& aa2 : aa2mass2count)
      {
        auto& mass2count = aa2.second;
        for (const auto& m2c : mass2count)
        {
          // normalize by counts
          double current_mass = m2c.first;
          size_t current_residue_count = m2c.second;
          size_t unmodified_residue_count = mass2count.begin()->second;
          double frequency_normalized = (double)current_residue_count / unmodified_residue_count;
          cout << aa2.first->getName() << "\t" << current_mass << "\t" << frequency_normalized << endl; // aa mass count
        }
      }
    }

#ifdef FILTER_AMBIGIOUS_PEAKS
    cout << "Frequent background mass shifts (mass vs. freq):" << endl;
    for (auto & hf : mass2high_frequency_)
    {
      cout << hf.first << "\t" << hf.second << endl;
    }
#endif

  }

  /* @brief Filter spectra to remove noise.

     - Remove zero intensities
     - Scale by root to reduce impact of high-intensity peaks
     - Normalize max intensity to 1.0
     - Remove isotopic peaks and determine charge
     - Set Unknown charge to z=1. Otherwise we get a lot of spurious matches 
     - Keep 20 highest-intensity peaks in 100 m/z windows
     - Keep max. 400 peaks per spectrum
       to highly charged fragments in the low m/z region
     - Calculate TIC of filtered spectrum
   */
  void preprocessSpectra_(PeakMap& exp, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    bool single_charge_spectra, 
    bool annotate_charge,
    double window_size,
    size_t peakcount)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      MSSpectrum & spec = exp[exp_index];

      // sort by mz
      spec.sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(spec, 
                                         0.01,
                                         false,
                                         1, 3, 
                                         false, 
                                         2, 10, 
                                         single_charge_spectra, 
                                         annotate_charge,
                                         false, // no iso peak count annotation
                                         true, // decreasing isotope model
                                         2, // enforce only starting from second peak
                                         true); // add up intensities
    }


    SqrtMower sqrt_mower_filter;
    sqrt_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", peakcount, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

    #ifdef DEBUG_OpenNuXL
    BinnedSpectrum peak_density(MSSpectrum(), 0.03, false, 0, 0);
    #endif

#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      MSSpectrum & spec = exp[exp_index];
      // sort by mz
      spec.sortByPosition();

      if (annotate_charge)
      { 
        // set Unknown charge to z=1. Otherwise we get a lot of spurious matches 
        // to highly charged fragments in the low m/z region
        DataArrays::IntegerDataArray& ia = spec.getIntegerDataArrays()[IA_CHARGE_INDEX]; // charge array
        for (int & z : ia) { if (z == 0) { z = 1; } }
      } 
    #ifdef DEBUG_OpenNuXL
      cout << "after deisotoping..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          cout  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
      cout << endl;
    #endif

      // remove noise
      window_mower_filter.filterPeakSpectrum(spec);

    #ifdef DEBUG_OpenNuXL
      cout << "after mower..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) cout << spec[i].getMZ() << "\t" << spec[i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          cout  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif
    
      nlargest_filter.filterPeakSpectrum(spec);

    #ifdef DEBUG_OpenNuXL
      cout << "after nlargest..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) cout << spec[i].getMZ() << "\t" << spec[i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          cout  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif
 
      // sort (nlargest changes order)
      spec.sortByPosition();
  
    #ifdef DEBUG_OpenNuXL
      cout << "after sort..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) cout << spec[i].getMZ() << "\t" << spec.getIntensity() << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          cout  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif

#ifdef DEBUG_OpenNuXL
      BinnedSpectrum bs(spec, 0.03, false, 0, 0);
      bs.getBins().coeffs().cwiseMax(1);

#pragma omp critical (peak_density_access)
      peak_density.getBins() += bs.getBins();
#endif

      // calculate TIC and store in float data array
      double TIC = std::accumulate(spec.begin(), spec.end(), 0.0, 
        [&](double a, const Peak1D& b) { return a + b.getIntensity(); });
      spec.getFloatDataArrays().clear();
      spec.getFloatDataArrays().resize(1);
      spec.getFloatDataArrays()[0].push_back(TIC);
      spec.getFloatDataArrays()[0].setName("TIC");
    }

#ifdef DEBUG_OpenNuXL
    ofstream dist_file;
    dist_file.open(getStringOption_("in") + ".fragment_dist.csv");
    dist_file << "m/z\tfragments" << "\n";
    for (double mz = 0; mz < 2500.0; mz+=0.03)
    {
      dist_file << mz << "\t" << peak_density.getBinIntensity(mz) << "\n";
    }
    dist_file.close();
#endif
    if (debug_level_ > 10) MzMLFile().store("debug_filtering.mzML", exp); 
  }

  void filterTopNAnnotations_(vector<vector<AnnotatedHit>>& ahs, const Size top_hits)
  {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)ahs.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      const Size topn = top_hits > ahs[scan_index].size() ? ahs[scan_index].size() : top_hits;
      std::partial_sort(ahs[scan_index].begin(), ahs[scan_index].begin() + topn, ahs[scan_index].end(), AnnotatedHit::hasBetterScore);
      ahs[scan_index].resize(topn);
      ahs[scan_index].shrink_to_fit();
    }
  }

  void rescoreFastHits_(
    const PeakMap& exp, 
    vector<vector<AnnotatedHit>>& annotated_hits,
    const RNPxlModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const RNPxlParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    auto param = partial_loss_spectrum_generator.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    partial_loss_spectrum_generator.setParameters(param);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      vector<AnnotatedHit> new_hits;

      // for each PSM of this spectrum
      for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
      {
        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, annotated_hits[scan_index][i].NA_mod_index); // advance to sum formula at index NA_mod_index
        const String precursor_rna_adduct = *mod_combinations_it->second.begin(); // get first nucleotide formula e.g. UU-H2O  TODO:is this correct? only handle the first?????????????????????????????????????????????????????????????????????
        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;
        
        if (precursor_rna_adduct == "none") 
        {
          new_hits.push_back(annotated_hits[scan_index][i]);
        }
        else
        {
          // if we have a cross-link, copy PSM information for each cross-linkable nucleotides
          for (auto const & c : feasible_MS2_adducts)
          {
            AnnotatedHit a(annotated_hits[scan_index][i]);
            a.cross_linked_nucleotide = c.first; // nucleotide
            new_hits.push_back(a);
          }
        }
      }
      annotated_hits[scan_index].swap(new_hits);
    }

    // fill in values of slow scoring so they can be used in percolator
    for (Size scan_index = 0; scan_index != annotated_hits.size(); ++scan_index)
    {
      // for each PSM of this spectrum
      for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
      {
        AnnotatedHit& ah = annotated_hits[scan_index][i];

        // reconstruct fixed and variable modified peptide sequence (without NA)
        const String& unmodified_sequence = ah.sequence.getString();
        AASequence aas = AASequence::fromString(unmodified_sequence);
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);
        AASequence fixed_and_variable_modified_peptide = all_modified_peptides[ah.peptide_mod_index]; 
        double current_peptide_mass_without_NA = fixed_and_variable_modified_peptide.getMonoWeight();

        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, ah.NA_mod_index);
        const String precursor_rna_adduct = *mod_combinations_it->second.begin();
        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;
        const vector<RNPxlFragmentAdductDefinition>& marker_ions = all_feasible_adducts.at(precursor_rna_adduct).marker_ions;
        const double precursor_rna_mass = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

        if (precursor_rna_adduct == "none") 
        {
          const double tags = exp[scan_index].getFloatDataArrays()[1][0];
          ah.score = OpenNuXL::calculateCombinedScore(ah, false, tags);
          continue;
        }

        // determine current nucleotide and associated partial losses
        vector<RNPxlFragmentAdductDefinition> partial_loss_modification;
        for (auto const & nuc_2_adducts : feasible_MS2_adducts)
        {
          if (nuc_2_adducts.first == ah.cross_linked_nucleotide)
          {
            partial_loss_modification = nuc_2_adducts.second;
          } 
        } 

        // TODO: not needed to generate all templates (but will make not much of a difference
        // as this is only done a few thousand times during post-scoring). That way,
        // the code is basically the same as in the main scoring loop.
        PeakSpectrum  partial_loss_template_z1, 
          partial_loss_template_z2, 
          partial_loss_template_z3;
     
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3); 

        PeakSpectrum marker_ions_sub_score_spectrum_z1, 
          partial_loss_spectrum_z1, 
          partial_loss_spectrum_z2;

        // nucleotide is associated with certain NA-related fragment losses?
        if (!partial_loss_modification.empty())
        {
          // shifted b- / y- / a-ions
          // generate shifted_immonium_ions_sub_score_spectrum.empty
          RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                      current_peptide_mass_without_NA,
                                      precursor_rna_adduct,
                                      precursor_rna_mass,
                                      1,
                                      partial_loss_modification,
                                      partial_loss_template_z1,
                                      partial_loss_template_z2,
                                      partial_loss_template_z3,
                                      partial_loss_spectrum_z1);

          RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                      current_peptide_mass_without_NA,
                                      precursor_rna_adduct,
                                      precursor_rna_mass,
                                      2, // don't know the charge of the precursor at that point
                                      partial_loss_modification,
                                      partial_loss_template_z1,
                                      partial_loss_template_z2,
                                      partial_loss_template_z3,
                                      partial_loss_spectrum_z2);
        }

        // add shifted marker ions
        marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
        marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation
        RNPxlFragmentIonGenerator::addMS2MarkerIons(
          marker_ions,
          marker_ions_sub_score_spectrum_z1,
          marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[IA_CHARGE_INDEX],
          marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

        const PeakSpectrum& exp_spectrum = exp[scan_index];
        float  partial_loss_sub_score(0), 
          marker_ions_sub_score(0),
          plss_MIC(0), 
          plss_err(1.0), 
          plss_Morph(0), 
          plss_modds(0);

//TODO: dieser Teil ist anders
        postScorePartialLossFragments_( unmodified_sequence.size(),
                                    exp_spectrum,
                                    fragment_mass_tolerance, 
                                    fragment_mass_tolerance_unit_ppm,
                                    partial_loss_spectrum_z1, 
                                    partial_loss_spectrum_z2,
                                    marker_ions_sub_score_spectrum_z1,
                                    partial_loss_sub_score,
                                    marker_ions_sub_score,
                                    plss_MIC, 
                                    plss_err, 
                                    plss_Morph,
                                    plss_modds);


        // fill in missing scores not considered in fast scoring
        ah.pl_MIC = plss_MIC;
        ah.pl_err = plss_err;
        ah.pl_Morph = plss_Morph;
        ah.pl_modds = plss_modds;
        // add extra matched ion current
// TODO: dieser Teil ist anders
        ah.total_MIC += plss_MIC + marker_ions_sub_score; 
        // scores from shifted peaks
        ah.marker_ions_score = marker_ions_sub_score;
        ah.partial_loss_score = partial_loss_sub_score;
        // combined score
        const double tags = exp[scan_index].getFloatDataArrays()[1][0];
        ah.score = OpenNuXL::calculateCombinedScore(ah, true, tags);
      } 
    } 
  }

  void annotateAndLocate_(
    const PeakMap& exp, 
    vector<vector<AnnotatedHit>>& annotated_hits,
    const RNPxlModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const RNPxlParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)   
  {
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    auto param = partial_loss_spectrum_generator.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    param.setValue("add_internal_fragments", "true"); 
    partial_loss_spectrum_generator.setParameters(param);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) { continue; }

      const PeakSpectrum & exp_spectrum = exp[scan_index];
      const Size & precursor_charge = exp_spectrum.getPrecursors()[0].getCharge();

      for (auto & a : annotated_hits[scan_index])
      {
        // get unmodified string
        const String unmodified_sequence = a.sequence.getString();

        // initialize result fields
        a.best_localization = unmodified_sequence;
        a.best_localization_score = 0;

        AASequence aas(AASequence::fromString(unmodified_sequence));

        // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

        // sequence with modifications - note: reannotated version requires much more memory heavy AASequence object
        const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[a.peptide_mod_index];
        const double fixed_and_variable_modified_peptide_weight = fixed_and_variable_modified_peptide.getMonoWeight();

        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.cbegin();
        std::advance(mod_combinations_it, a.NA_mod_index);
        const String precursor_rna_adduct = *mod_combinations_it->second.begin(); // TODO: check if it is enough to consider only first precursor adduct ????????????????????????????????????????????????????????
        const double precursor_rna_mass = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

        // we don't localize on non-cross-links
        if (precursor_rna_adduct == "none") { continue; }

        // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
        // 1. get all possible NA fragment shifts in the MS2 (based on the precursor RNA/DNA)
        OPENMS_LOG_DEBUG << "Precursor NA adduct: "  << precursor_rna_adduct << endl;

        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;

        if (feasible_MS2_adducts.empty()) { continue; } // should not be the case - check case of no nucleotide but base fragment ?

        // 2. retrieve the (nucleotide specific) fragment adducts for the cross-linked nucleotide (annotated in main search)
        auto nt_to_adducts = std::find_if(feasible_MS2_adducts.begin(), feasible_MS2_adducts.end(),
          [&a](NucleotideToFeasibleFragmentAdducts const & item)
          {
            return (item.first == a.cross_linked_nucleotide);
          });

        OPENMS_POSTCONDITION(nt_to_adducts != feasible_MS2_adducts.end(), "Nucleotide not found in mapping to feasible adducts.")

        const vector<RNPxlFragmentAdductDefinition>& partial_loss_modification = nt_to_adducts->second;

        // get marker ions (these are not specific to the cross-linked nucleotide but also depend on the whole oligo bound to the precursor)
        const vector<RNPxlFragmentAdductDefinition>& marker_ions = all_feasible_adducts.at(precursor_rna_adduct).marker_ions;
        OPENMS_LOG_DEBUG << "Marker ions used for this Precursor NA adduct: "  << endl;
        for (auto & fa : marker_ions)
        {
          OPENMS_LOG_DEBUG << fa.name << " " << fa.mass << endl;
        }

        // generate total loss spectrum for the fixed and variable modified peptide (without NAs) (using the settings for partial loss generation)
        // but as we also add the abundant immonium ions for charge 1 and precursor ions for all charges to get a more complete annotation
        // (these have previously not been used in the scoring of the total loss spectrum)
        PeakSpectrum total_loss_spectrum;
        {
          TheoreticalSpectrumGenerator tmp_generator;
          Param new_param(partial_loss_spectrum_generator.getParameters());
          new_param.setValue("add_all_precursor_charges", "true");
          new_param.setValue("add_abundant_immonium_ions", "true");
          new_param.setValue("add_losses", "true");
          new_param.setValue("add_term_losses", "true");
          new_param.setValue("add_a_ions", "true");
          new_param.setValue("add_internal_fragments", "true");
          tmp_generator.setParameters(new_param);
          tmp_generator.getSpectrum(total_loss_spectrum, fixed_and_variable_modified_peptide, 1, precursor_charge);
        }

        // add special immonium ions
        RNPxlFragmentIonGenerator::addSpecialLysImmonumIons(
          unmodified_sequence,
          total_loss_spectrum, 
          total_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
          total_loss_spectrum.getStringDataArrays()[0]);
        total_loss_spectrum.sortByPosition(); // need to resort after adding special immonium ions

        PeakSpectrum partial_loss_spectrum, 
          partial_loss_template_z1, 
          partial_loss_template_z2, 
          partial_loss_template_z3;
       
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3); 
        RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                    fixed_and_variable_modified_peptide_weight,
                                    precursor_rna_adduct,
                                    precursor_rna_mass,
                                    precursor_charge,
                                    partial_loss_modification,
                                    partial_loss_template_z1,
                                    partial_loss_template_z2,
                                    partial_loss_template_z3,
                                    partial_loss_spectrum);

         // add shifted marker ions
         RNPxlFragmentIonGenerator::addMS2MarkerIons(
           marker_ions,
           partial_loss_spectrum,
           partial_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
           partial_loss_spectrum.getStringDataArrays()[0]);

        partial_loss_spectrum.sortByPosition(); // need to resort after adding marker ions

        // fill annotated spectrum information
        set<Size> peak_is_annotated;  // experimental peak index

        // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)
        using MapIonIndexToFragmentAnnotation = map<Size, vector<RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_> >;
        MapIonIndexToFragmentAnnotation unshifted_b_ions, unshifted_y_ions, unshifted_a_ions, shifted_b_ions, shifted_y_ions, shifted_a_ions;
        vector<PeptideHit::PeakAnnotation> shifted_immonium_ions,
          unshifted_loss_ions,
          annotated_marker_ions,
          annotated_precursor_ions,
          annotated_immonium_ions;

        // first annotate total loss peaks (these give no information where the actual shift occured)
        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << "Annotating ion (total loss spectrum): " << fixed_and_variable_modified_peptide.toString()  << endl;
        #endif
        vector<pair<Size, Size>> alignment;

        // align spectra (only allow matching charges)
        DataArrays::FloatDataArray ppm_error_array; // not needed here but filled by alignment
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, 
          fragment_mass_tolerance, 
          fragment_mass_tolerance_unit_ppm, 
          total_loss_spectrum, 
          exp_spectrum, 
          total_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX], 
          exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX], 
          ppm_error_array);

        const PeakSpectrum::StringDataArray& total_loss_annotations = total_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& total_loss_charges = total_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX];

        for (auto const & aligned : alignment)
        {
          // information on the experimental fragment in the alignment
          const Size& fragment_index = aligned.second;
          const Peak1D& fragment = exp_spectrum[fragment_index];
          const double fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double fragment_mz = fragment.getMZ();
         

          const String& ion_name = total_loss_annotations[aligned.first];
          const int charge = total_loss_charges[aligned.first];

          OPENMS_PRECONDITION(exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX][fragment_index] == charge, "Charges in alignment must match.");

          // define which ion names are annotated
          if (ion_name[0] == 'y')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
            const bool ion_has_neutral_loss = (loss_first != string::npos);

            if (ion_has_neutral_loss) // ion with neutral loss e.g. water
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = charge;
              fa.annotation = ion_name;
              unshifted_loss_ions.push_back(fa);
              peak_is_annotated.insert(aligned.second);
            }
            else // no neutral loss
            {
              String ion_nr_string = ion_name.substr(1, charge_pos - 1);
              Size ion_number = (Size)ion_nr_string.toInt();
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_y_ions[ion_number].push_back(d);
              #ifdef DEBUG_OpenNuXL
                const AASequence& peptide_sequence = fixed_and_variable_modified_peptide.getSuffix(ion_number);
                OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              peak_is_annotated.insert(aligned.second);
            }
          }
          else if (ion_name[0] == 'b')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
            const bool ion_has_neutral_loss = (loss_first != string::npos);

            if (ion_has_neutral_loss)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = charge;
              fa.annotation = ion_name;
              unshifted_loss_ions.push_back(fa);
              peak_is_annotated.insert(aligned.second);
            }
            else
            {
              String ion_nr_string = ion_name.substr(1, charge_pos - 1);
              Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_OpenNuXL
                const AASequence& peptide_sequence = aas.getPrefix(ion_number);
                OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_b_ions[ion_number].push_back(d);
              peak_is_annotated.insert(aligned.second);
            }
          }
          else if (ion_name[0] == 'a')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
            const bool ion_has_neutral_loss = (loss_first != string::npos);

            if (ion_has_neutral_loss)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = charge;
              fa.annotation = ion_name;
              unshifted_loss_ions.push_back(fa);
              peak_is_annotated.insert(aligned.second);
            }
            else
            {
              String ion_nr_string = ion_name.substr(1, charge_pos - 1);
              auto ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_OpenNuXL
                const AASequence& peptide_sequence = aas.getPrefix(ion_number);
                OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_a_ions[ion_number].push_back(d);
              peak_is_annotated.insert(aligned.second);
            }
          }
          else if (ion_name.hasPrefix("[M+")) // precursor ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            fa.annotation = ion_name;
            annotated_precursor_ions.push_back(fa);
            peak_is_annotated.insert(aligned.second);
          }
          else if (ion_name.hasPrefix("i")) // immonium ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            fa.annotation = ion_name;
            annotated_immonium_ions.push_back(fa);
            peak_is_annotated.insert(aligned.second);
          }
          else if (isupper(ion_name[0])) // internal ions
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            fa.annotation = ion_name;
            annotated_immonium_ions.push_back(fa);  //TODO: add to annotated_internal_fragment_ions or rename vector
            peak_is_annotated.insert(aligned.second);
          }
        }

        // generate fragment annotation strings for unshifted ions
        vector<PeptideHit::PeakAnnotation> fas;
        if (!unshifted_b_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("b", unshifted_b_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!unshifted_y_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("y", unshifted_y_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!unshifted_a_ions.empty())
        {
          const vector<PeptideHit::PeakAnnotation>& fas_tmp = RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("a", unshifted_a_ions);
          fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
        }
        if (!annotated_immonium_ions.empty())
        {
          fas.insert(fas.end(), annotated_immonium_ions.begin(), annotated_immonium_ions.end());          
        }
        if (!unshifted_loss_ions.empty())
        {
          fas.insert(fas.end(), unshifted_loss_ions.begin(), unshifted_loss_ions.end());          
        }

        vector<double> sites_sum_score(aas.size(), 0);

        /////////////////
        // Align partial-loss-spectrum to the experimental measured one
        alignment.clear();
        ppm_error_array.clear();

        // align spectra (only allow matching charges)
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, 
          fragment_mass_tolerance, 
          fragment_mass_tolerance_unit_ppm, 
          partial_loss_spectrum, 
          exp_spectrum, 
          partial_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX], 
          exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX], 
          ppm_error_array);

        #ifdef OPENUXL_DEBUG
        for (size_t i = 0; i != exp_spectrum.size(); ++i)
        {
          OPENMS_LOG_DEBUG << "exp: " << exp_spectrum[i].getMZ() << "\t" << exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX][i] << endl;
        }
        #endif

        const PeakSpectrum::StringDataArray& partial_loss_annotations = partial_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& partial_loss_charges = partial_loss_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX];

        if (alignment.empty())
        {
          a.fragment_annotations = fas;
          continue;
        }

        /* uncomment to write all annotations to a file(only makes sense if a single spectrum is searched)
        MSExperiment tmp_exp;
        tmp_exp.addSpectrum(total_loss_spectrum);
        tmp_exp.addSpectrum(partial_loss_spectrum);
        MzMLFile().store("theoretical_loss_spectrum.mzML", tmp_exp);
        */

        for (auto pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
        {
          // only annotate experimental peaks with shift - i.e. do not annotated complete loss peaks again
          if (peak_is_annotated.find(pair_it->second) != peak_is_annotated.end()) { continue; }

          // information on the experimental fragment in the alignment
          const Size & fragment_index = pair_it->second;
          const Peak1D & fragment = exp_spectrum[fragment_index];
          const double & fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double & fragment_mz = fragment.getMZ();
          const int & fragment_charge = exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX][fragment_index];
          
          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "fragment_mz:" << fragment_mz << " fragment_charge:" << fragment_charge << endl; 
          #endif

          String ion_name = partial_loss_annotations[pair_it->first];
          const int charge = partial_loss_charges[pair_it->first];

          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "theo_name:" << ion_name  << " theo_charge:" << charge << endl; 
          #endif
          vector<String> f;

          ion_name.split(' ', f);  // e.g. "y3 C3O" or just "y2"
          String fragment_shift_name;
          if (f.size() == 2) { fragment_shift_name = f[1]; }

          String fragment_ion_name = f[0]; // e.g. y3

          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
          #endif

          // define which ion names are annotated
          if (fragment_ion_name.hasPrefix("y"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("y", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            shifted_y_ions[ion_number].push_back(d);
          }
          else if (fragment_ion_name.hasPrefix("b"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("b", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            shifted_b_ions[ion_number].push_back(d);
          }
          else if (fragment_ion_name.hasPrefix("a"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("a", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            shifted_a_ions[ion_number].push_back(d);
          }
          else if (ion_name.hasPrefix(RNPxlFragmentIonGenerator::ANNOTATIONS_MARKER_ION_PREFIX))
          {
            OPENMS_LOG_DEBUG << "Marker ion aligned: " << ion_name << " fragment_mz: " << fragment_mz << " fragment_charge: " << fragment_charge << endl;
            if (fragment_charge == 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              annotated_marker_ions.push_back(fa);
            }
          }
          else if (ion_name.hasPrefix("i"))
          {
            if (fragment_charge == 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              shifted_immonium_ions.push_back(fa);
            }
          }
          else if (ion_name.hasPrefix("[M+"))
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            fa.annotation = ion_name;  
            annotated_precursor_ions.push_back(fa);
          }
          else if (isupper(ion_name[0])) // shifted internal ions
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            String with_plus = ion_name;
            with_plus.substitute(' ', '+');
            fa.annotation = with_plus; // turn "PEPT U-H2O" into "PEPT+U-H20"
            shifted_immonium_ions.push_back(fa);  //TODO: add to shfited_internal_fragment_ions or rename vector
          }
        }

        // track shifts in n- and c-term ladders (in AAs coordinates)
        // n_shifts and c_shifts will contain the summed intensities over all observed shifts at that position
        // the distinction allows to easily detect prefix and suffix ladders in the next step
        vector<double> n_shifts(sites_sum_score.size(), 0);
        vector<double> c_shifts(sites_sum_score.size(), 0);

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_b_ions.find(i + 1) == shifted_b_ions.end()) { continue; }
          for (auto& k : shifted_b_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_a_ions.find(i + 1) == shifted_a_ions.end()) { continue; }
          for (auto& k : shifted_a_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_shifts.size(); ++i)
        {
          const Size ion_index = c_shifts.size() - i;
          if (shifted_y_ions.find(ion_index) == shifted_y_ions.end()) { continue; }
          for (auto& k : shifted_y_ions[ion_index]) { c_shifts[i] += k.intensity; }
        }

        vector<double> n_noshifts(sites_sum_score.size(), 0);
        vector<double> c_noshifts(sites_sum_score.size(), 0);
        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_b_ions.find(i + 1) == unshifted_b_ions.end()) { continue; }
          for (auto& k : unshifted_b_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_a_ions.find(i + 1) == unshifted_a_ions.end()) { continue; }
          for (auto& k : unshifted_a_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_noshifts.size(); ++i)
        {
          const Size ion_index = c_noshifts.size() - i;
          if (unshifted_y_ions.find(ion_index) == unshifted_y_ions.end()) { continue; }
          for (auto& k : unshifted_y_ions[ion_index]) { c_noshifts[i] += k.intensity; }
        }

#ifdef DEBUG_OpenNuXL
        cout << "n:";
        for (auto& k : n_shifts) cout << k << " ";
        cout << endl;
        cout << "c:";
        for (auto& k : c_shifts) cout << k << " ";
        cout << endl;
        cout << "n0:";
        for (auto& k : n_noshifts) cout << k << " ";
        cout << endl;
        cout << "c0:";
        for (auto& k : c_noshifts) cout << k << " ";
        cout << endl;
#endif

        // Rules implemented:
        // 1. if cross-link on AA, then the prefix or suffix ending at this AA must be shifted
        // 2. if the previous AA in the prefix / suffix had a stronger shifted signal, then the current on is not the correct one
        // 3. if the current AA is cross-linked, then the previous AA is not cross-linked and we should observe an unshifted prefix / suffix ion
        // 4. if the current AA is cross-linked, we should observe a shifted prefix / suffix ion for the next AA, too
        // 5. Sum up all intensities of shifted prefix / suffix ions
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          sites_sum_score[i] = 0.0;
          if (n_shifts[i] == 0 && c_shifts[i] == 0) { continue; } // no shifts? no cross-link at this AA

          if (n_shifts[i] > 0)
          {
            if (i >= 1 && n_shifts[i - 1] > n_shifts[i]) continue; // Strong signal from shifted AA before the current one? Then skip it.
            if (i >= 1 && n_noshifts[i - 1] == 0) continue; // continue if unshifted AA is missing before (left of) the shifted one.
            if (i < n_shifts.size()-1 && n_shifts[i + 1] == 0) continue; // Need a shifted ladder after (maybe too conservative?)
            // sum up all intensities from this position and all longer prefixes that also carry the NA
            for (Size j = i; j != sites_sum_score.size(); ++j) { sites_sum_score[i] += n_shifts[j]; }
          }

          if (c_shifts[i] > 0)
          {
            if (i < c_shifts.size()-1 && c_shifts[i + 1] > c_shifts[i]) continue; // AA after already shifted. So ignore this one.
            if (i < c_noshifts.size()-1 && c_noshifts[i + 1] == 0) continue; // continue if unshifted AA is missing before (right of) the shifted one.
            if (i >=1 && c_shifts[i - 1] == 0) continue; // Need a shifted ladder after (maybe too conservative?)
            sites_sum_score[i] += c_shifts[i];
            // sum up all intensities from this position and all longer suffixes that also carry the NA
            for (int j = i; j >= 0; --j) { sites_sum_score[i] += c_shifts[j]; }
          }
        }
#ifdef DEBUG_OpenNuXL
        cout << "site sum score (shifted a/b/y-ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << "Localisation based on immonium ions: ";
        #endif
        String aas_unmodified = aas.toUnmodifiedString();
        for (Size i = 0; i != aas_unmodified.size(); ++i)
        {
          String origin = String(aas_unmodified[i]);

          for (auto& a : shifted_immonium_ions)
          {
            // compare origin (the AA) of immonium ion to current AA
            if (a.annotation[1] == aas_unmodified[i])
            {
              sites_sum_score[i] += a.intensity;
            }
          }
        }
#ifdef DEBUG_OpenNuXL
        cout << "site sum score (shifted a/b/y-ions & immonium ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        String best_localization = unmodified_sequence;
        int best_localization_position = -1; // UNKNOWN
        double best_localization_score = 0;
        String localization_scores;
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          if (sites_sum_score[i] > best_localization_score) { best_localization_score = sites_sum_score[i]; }
        }

        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << String::number(100.0 * sites_sum_score[i], 2);
          #endif

          if (i != 0) localization_scores += ',';
          if (sites_sum_score[i] > 0 )
          {
            localization_scores += String::number(100.0 * sites_sum_score[i], 2);
          }
          else
          {
            localization_scores += "0";
          }

          if (best_localization_score > 0.0 && sites_sum_score[i] >= best_localization_score - 1e-6)
          {
            best_localization[i] = tolower(best_localization[i]);
            best_localization_position = i; // Note: check if there are situations where multiple have the same score
          }
        }
        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << endl;
        #endif

        // create annotation strings for shifted fragment ions
        RNPxlFragmentAnnotationHelper::addShiftedPeakFragmentAnnotation_(shifted_b_ions,
                                          shifted_y_ions,
                                          shifted_a_ions,
                                          shifted_immonium_ions,
                                          annotated_marker_ions,
                                          annotated_precursor_ions,
                                          fas);

        // store score of best localization(s)
        a.localization_scores = localization_scores;
        a.best_localization = best_localization;
        a.best_localization_score = best_localization_score;
        a.best_localization_position = best_localization_position;
        a.fragment_annotations = fas;

        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << "Ion centric annotation: " << endl;
          OPENMS_LOG_DEBUG << "unshifted b ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", unshifted_b_ions) << endl;
          OPENMS_LOG_DEBUG << "unshifted y ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", unshifted_y_ions) << endl;
          OPENMS_LOG_DEBUG << "unshifted a ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", unshifted_a_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted b ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", shifted_b_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted y ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", shifted_y_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted a ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", shifted_a_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted immonium ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(shifted_immonium_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted marker ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(annotated_marker_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted precursor ions: " << endl;
          OPENMS_LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(annotated_precursor_ions) << endl;
          OPENMS_LOG_DEBUG << "Localization scores: ";
          OPENMS_LOG_DEBUG << localization_scores << endl;
          OPENMS_LOG_DEBUG << "Localisation based on ion series and immonium ions of all observed fragments: ";
          OPENMS_LOG_DEBUG << best_localization << endl;
        #endif
      }
    }
  }

  /* @brief Localization step of the cross-link identification engine.
    1. Generates all fragment adducts based on the attached precursor adduct
    2. Calculates an additive score that considers the presence or absence of evidence for a cross-linking site
    3. Add additional meta information for PSM.
   */
  void postScoreHits_(const PeakMap& exp, 
                      vector<vector<AnnotatedHit> >& annotated_XL_hits, 
                      vector<vector<AnnotatedHit> >& annotated_peptide_hits, 
                      const RNPxlModificationMassesResult& mm, 
                      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
                      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
                      Size max_variable_mods_per_peptide, 
                      double fragment_mass_tolerance, 
                      bool fragment_mass_tolerance_unit_ppm, 
                      const RNPxlParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    assert(exp.size() == annotated_XL_hits.size());
    assert(exp.size() == annotated_peptide_hits.size());

    // If we did a (total-loss) only fast scoring, PSMs were not associated with a nucleotide.
    // To make the localization code work for both fast and slow (all-shifts) scoring,
    // we copy PSMs for every cross-linkable nucleotide present in the precursor.
    // Then we recalculate the XL specific scores not conisdered in fast scoring
    if (fast_scoring_)
    {
      rescoreFastHits_(exp, annotated_XL_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,  all_feasible_adducts);
      rescoreFastHits_(exp, annotated_peptide_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, all_feasible_adducts);
    }

    annotateAndLocate_(exp, annotated_XL_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,  all_feasible_adducts);
    annotateAndLocate_(exp, annotated_peptide_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, all_feasible_adducts);
  }

  void fillSpectrumID_(
    const vector<AnnotatedHit>& ahs, 
    PeptideIdentification& pi, 
    const RNPxlModificationMassesResult& mm, 
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    const Size max_variable_mods_per_peptide,
    const Size scan_index, 
    const MSSpectrum& spec,
    const map<String, PrecursorPurity::PurityScores>& purities,
    const vector<size_t>& nr_candidates)
  {
    pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
    pi.setMetaValue("spectrum_reference", spec.getNativeID());
    pi.setScoreType("NuXLScore");
    pi.setHigherScoreBetter(true);
    pi.setRT(spec.getRT());
    pi.setMZ(spec.getPrecursors()[0].getMZ());
    double precursor_intensity_log10 = log10(1.0 + spec.getPrecursors()[0].getIntensity());
    pi.setMetaValue("precursor_intensity_log10", precursor_intensity_log10);
    Size charge = spec.getPrecursors()[0].getCharge();

    // create full peptide hit structure from annotated hits
    vector<PeptideHit> phs = pi.getHits();
    for (auto const & ah : ahs)
    {
      PeptideHit ph;
      ph.setCharge(charge);

      // get unmodified string
      const String & s = ah.sequence.getString();

      OPENMS_POSTCONDITION(!s.empty(), "Error: empty sequence in annotated hits.");
      AASequence aas = AASequence::fromString(s);

      // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
      vector<AASequence> all_modified_peptides;
      ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
      ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

      // reannotate much more memory heavy AASequence object
      AASequence fixed_and_variable_modified_peptide = all_modified_peptides[ah.peptide_mod_index];           
      ph.setScore(ah.score);
      ph.setMetaValue(String("NuXL:score"), ah.score); // important for Percolator feature set because the PeptideHit score might be overwritten by a q-value

      // - # of variable mods 
      // - Phosphopeptide
      int is_phospho(0);
      int n_var_mods = 0;
      for (Size i = 0; i != fixed_and_variable_modified_peptide.size(); ++i)
      { 
        const Residue& r = fixed_and_variable_modified_peptide[i];
        if (!r.isModified()) continue;
        if (variable_modifications.val.find(r.getModification()) != variable_modifications.val.end())
        {
          ++n_var_mods;
        }

        if (r.getModification()->getId() == "Phospho") { is_phospho = 1; }
      }
      auto n_term_mod = fixed_and_variable_modified_peptide.getNTerminalModification();
      auto c_term_mod = fixed_and_variable_modified_peptide.getCTerminalModification();

      if (n_term_mod != nullptr &&
        variable_modifications.val.find(n_term_mod) != variable_modifications.val.end()) ++n_var_mods;
      if (c_term_mod != nullptr &&
        variable_modifications.val.find(c_term_mod) != variable_modifications.val.end()) ++n_var_mods;

      ph.setMetaValue(String("variable_modifications"), n_var_mods);
      ph.setMetaValue(String("n_theoretical_peaks"), ah.n_theoretical_peaks);

      // determine NA modification from index in map
      auto mod_combinations_it = mm.mod_combinations.cbegin();
      std::advance(mod_combinations_it, ah.NA_mod_index);
      ph.setMetaValue(String("NuXL:mass_error_p"), ah.mass_error_p);
      ph.setMetaValue(String("NuXL:total_loss_score"), ah.total_loss_score);
      ph.setMetaValue(String("NuXL:immonium_score"), ah.immonium_score);
      ph.setMetaValue(String("NuXL:precursor_score"), ah.precursor_score);
      ph.setMetaValue(String("NuXL:marker_ions_score"), ah.marker_ions_score);
      ph.setMetaValue(String("NuXL:partial_loss_score"), ah.partial_loss_score);

      // total loss and partial loss (pl) related subscores (matched ion current, avg. fragment error, morpheus score)
      ph.setMetaValue(String("NuXL:MIC"), ah.MIC);
      ph.setMetaValue(String("NuXL:err"), ah.err);
      ph.setMetaValue(String("NuXL:Morph"), ah.Morph);
      ph.setMetaValue(String("NuXL:modds"), ah.modds);
      ph.setMetaValue(String("NuXL:pl_MIC"), ah.pl_MIC);
      ph.setMetaValue(String("NuXL:pl_err"), ah.pl_err);
      ph.setMetaValue(String("NuXL:pl_Morph"), ah.pl_Morph);
      ph.setMetaValue(String("NuXL:pl_modds"), ah.pl_modds);
      ph.setMetaValue(String("NuXL:pl_pc_MIC"), ah.pl_pc_MIC);
      ph.setMetaValue(String("NuXL:pl_im_MIC"), ah.pl_im_MIC);
      ph.setMetaValue(String("NuXL:total_Morph"), ah.Morph + ah.pl_Morph);
      ph.setMetaValue(String("NuXL:total_HS"), ah.total_loss_score + ah.partial_loss_score);

      ph.setMetaValue(String("NuXL:tag_XLed"), ah.tag_XLed);
      ph.setMetaValue(String("NuXL:tag_unshifted"), ah.tag_unshifted);
      ph.setMetaValue(String("NuXL:tag_shifted"), ah.tag_shifted);
      
      ph.setMetaValue(String("NuXL:total_MIC"), ah.total_MIC);  // fraction of matched ion current from total + partial losses

      String NA = *mod_combinations_it->second.begin(); // TODO: consider only first (potentially ambigious one) ??????????????????????????????????????
      ph.setMetaValue(String("NuXL:NA"), NA); // return first nucleotide formula matching the index of the empirical formula

      double na_mass_z0 = EmpiricalFormula(mod_combinations_it->first).getMonoWeight(); // NA uncharged mass via empirical formula
      // length of oligo
      size_t NA_length = NA.find_first_of("+-");
      if (NA_length == std::string::npos)
      {
        if (na_mass_z0 > 0)
        {
          ph.setMetaValue(String("NuXL:NA_length"), NA.size());
        }
        else
        {
          ph.setMetaValue(String("NuXL:NA_length"), 0);
        }
      }
      else
      {
        ph.setMetaValue(String("NuXL:NA_length"), NA_length);
      }

      ph.setMetaValue("NuXL:NT", String(ah.cross_linked_nucleotide));  // the cross-linked nucleotide
      ph.setMetaValue("NuXL:NA_MASS_z0", na_mass_z0); // NA uncharged mass via empirical formula
      ph.setMetaValue("NuXL:isXL", na_mass_z0 > 0); 
      ph.setMetaValue("NuXL:isPhospho", is_phospho); 

      ph.setMetaValue("NuXL:best_localization_score", ah.best_localization_score);
      if (!ah.localization_scores.empty())
      {
        ph.setMetaValue("NuXL:localization_scores", ah.localization_scores);
      }
      else
      {
        ph.setMetaValue("NuXL:localization_scores", "NA");
      }
      ph.setMetaValue("NuXL:best_localization", ah.best_localization);
      ph.setMetaValue("NuXL:best_localization_position", ah.best_localization_position);

      // one-hot encoding of cross-linked nucleotide
      const String can_cross_link = getStringOption_("RNPxl:can_cross_link");
      for (const auto& c : can_cross_link) 
      {
        if (c == ah.cross_linked_nucleotide)
        { 
          ph.setMetaValue(String("NuXL:XL_" + String(c)), 1);
        }
        else
        {
          ph.setMetaValue(String("NuXL:XL_" + String(c)), 0);
        }
      }

      // also annotate PI to hit so it is available to percolator
      ph.setMetaValue("precursor_intensity_log10", precursor_intensity_log10);

      if (!purities.empty())
      {
        ph.setMetaValue("precursor_purity", purities.at(spec.getNativeID()).signal_proportion);
      }

      ph.setMetaValue("nucleotide_mass_tags", (double)spec.getFloatDataArrays()[1][0]);
      int maxtag = spec.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0];
      ph.setMetaValue("NuXL:aminoacid_max_tag", maxtag);
      const double id2maxtag = maxtag == 0 ? 0 : (ah.ladder_score * s.size()) / (double)maxtag; 
      ph.setMetaValue("NuXL:aminoacid_id_to_max_tag_ratio", id2maxtag);
      ph.setMetaValue("nr_candidates", nr_candidates[scan_index]);
      ph.setMetaValue("NuXL:explained_peak_fraction", ah.explained_peak_fraction);
      ph.setMetaValue("NuXL:wTop50", ah.wTop50);

      ph.setPeakAnnotations(ah.fragment_annotations);
      ph.setMetaValue("isotope_error", static_cast<int>(ah.isotope_error));
      ph.setMetaValue(String("NuXL:ladder_score"), ah.ladder_score);
      ph.setMetaValue(String("NuXL:sequence_score"), ah.sequence_score);
      ph.setMetaValue(String("CalcMass"), + (fixed_and_variable_modified_peptide.getMonoWeight(Residue::Full, charge) + na_mass_z0)/charge); // overwrites CalcMass in PercolatorAdapter
      // set the amino acid sequence (for complete loss spectra this is just the variable and modified peptide. For partial loss spectra it additionally contains the loss induced modification)
      ph.setSequence(fixed_and_variable_modified_peptide);
      phs.push_back(ph);  // add new hit
    }

    pi.setHits(phs);
    pi.assignRanks();    

    // assign (unique) ranks
    phs = pi.getHits();
    for (Size r = 0; r != phs.size(); ++r) { phs[r].setMetaValue("rank", static_cast<int>(r)); }
    pi.setHits(phs);
  }


  /**
    1. Reconstruct original peptide from memory efficient structure
    2. Add additional meta information for PSM.
  */
  void postProcessHits_(const PeakMap& exp, 
    vector<vector<AnnotatedHit> >& annotated_XL_hits, 
    vector<vector<AnnotatedHit> >& annotated_peptide_hits, 
    vector<ProteinIdentification>& protein_ids, 
    vector<PeptideIdentification>& peptide_ids, 
    const RNPxlModificationMassesResult& mm, 
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide,
    const map<String, PrecursorPurity::PurityScores>& purities,
    const vector<size_t>& nr_candidates)
  {
    assert(annotated_XL_hits.size() == annotated_peptide_hits.size());
    SignedSize hit_count = static_cast<SignedSize>(annotated_XL_hits.size());

    for (SignedSize scan_index = 0; scan_index < hit_count; ++scan_index)
    {
      const MSSpectrum& spec = exp[scan_index];
      vector<AnnotatedHit>& ahs_XL = annotated_XL_hits[scan_index];
      vector<AnnotatedHit>& ahs_peptide = annotated_peptide_hits[scan_index];

      if (ahs_XL.empty() && ahs_peptide.empty()) continue;

      // create empty PeptideIdentification object and fill meta data
      peptide_ids.push_back(PeptideIdentification());

      if (!ahs_XL.empty())
      {
        fillSpectrumID_(
          ahs_XL, 
          peptide_ids.back(), // append hits
          mm,
          fixed_modifications, 
          variable_modifications, 
          max_variable_mods_per_peptide,
          scan_index, 
          spec,
          purities,
          nr_candidates);
      }

      if (!ahs_peptide.empty())
      {
        fillSpectrumID_(
          ahs_peptide, 
          peptide_ids.back(), // append hits
          mm, 
          fixed_modifications, 
          variable_modifications, 
          max_variable_mods_per_peptide,
          scan_index, 
          spec,
          purities,
          nr_candidates);
      }
    }
    // hits have rank and are sorted by score

    map<String, Size> sequence_is_topPSM;
    map<String, set<int>> sequence_charges; // of top PSM
    map<String, Size> sequence_is_XL;
    map<String, Size> sequence_is_peptide;
    for (const auto & pid : peptide_ids)
    {
      if (pid.getHits().empty()) continue;
      const auto & top_hit = pid.getHits()[0];
      const String& unmodified_sequence = top_hit.getSequence().toUnmodifiedString();
      ++sequence_is_topPSM[unmodified_sequence];
      sequence_charges[unmodified_sequence].insert(top_hit.getCharge());
      if (static_cast<int>(top_hit.getMetaValue("NuXL:isXL")) == 1)
      {
        ++sequence_is_XL[unmodified_sequence];
      }
      else
      {
        ++sequence_is_peptide[unmodified_sequence];
      }
    }
    for (auto & pid : peptide_ids)
    {
      for (auto & ph : pid.getHits())
      {
        const String& unmodified_sequence = ph.getSequence().toUnmodifiedString();
        if (sequence_is_topPSM.find(unmodified_sequence) != sequence_is_topPSM.end())
        {  
          ph.setMetaValue("CountSequenceIsTop", sequence_is_topPSM[unmodified_sequence]);
          ph.setMetaValue("CountSequenceCharges", sequence_charges[unmodified_sequence].size());
          ph.setMetaValue("CountSequenceIsXL", sequence_is_XL[unmodified_sequence]);
          ph.setMetaValue("CountSequenceIsPeptide", sequence_is_peptide[unmodified_sequence]);
        }
      }
    }
    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenNuXL");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    search_parameters.charges = String(getIntOption_("precursor:min_charge")) + ":" + String(getIntOption_("precursor:max_charge"));
    search_parameters.fixed_modifications = getStringList_("modifications:fixed");
    search_parameters.variable_modifications = getStringList_("modifications:variable");
    search_parameters.missed_cleavages = getIntOption_("peptide:missed_cleavages");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("peptide:enzyme"));

#ifdef OPENNUXL_SEPARATE_FEATURES
    for (auto & pi : peptide_ids)
    {
      for (auto & ph : pi.getHits())
      {
        if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
        {
          // fill XL related feature columns with zeros
          for (const auto & s : feature_set_)
          {
            if (s.hasPrefix("NuXL")) ph.setMetaValue("XL_" + s, 0.0);
          }
        }
        else  // XL
        {
          // fill linear peptide releated feature columns with zero (after copying them over)
          for (const auto& s : feature_set)
          {
            if (s.hasPrefix("NuXL"))
            { 
              ph.setMetaValue("XL_" + s, static_cast<double>(ph.getMetaValue(s)));
              ph.setMetaValue(s, 0.0);
            }
          }
        }
      }
    }

    // we duplicated the feature set
    auto XL_columns = feature_set;
    for (const auto& s : feature_set) { if (s.hasPrefix("NuXL")) XL_columns.push_back("XL_" + s); }
    feature_set_ = XL_columns;
#endif

    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set_, ","));

    protein_ids[0].setSearchParameters(search_parameters);
  }

  void mapPrecursorMassesToScans(const Int min_precursor_charge,
                                 const Int max_precursor_charge,
                                 const IntList &precursor_isotopes,
                                 const double small_peptide_mass_filter_threshold,
                                 const Size peptide_min_size,
                                 const PeakMap & spectra,
                                 multimap<double, pair<Size, int>> & multimap_mass_2_scan_index) const
  {
    Size fractional_mass_filtered(0), small_peptide_mass_filtered(0);

    for (MSExperiment::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();

        if (precursor_charge < min_precursor_charge
         || precursor_charge > max_precursor_charge)
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();

        // map (corrected) precursor mass to spectra
        for (int i : precursor_isotopes)
        {
          double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

          // corrected for monoisotopic misassignments of the precursor annotation
          if (i != 0) { precursor_mass -= i * Constants::C13C12_MASSDIFF_U; }

          if (getFlag_("RNPxl:filter_fractional_mass"))
          {
            if (precursor_mass < 1750.0 && precursor_mass - floor(precursor_mass) < 0.2)
            {
              fractional_mass_filtered++;
              continue;
            }
          }

          if (precursor_mass < small_peptide_mass_filter_threshold)
          {
            small_peptide_mass_filtered++;
            continue;
          }

          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, make_pair(scan_index, i)));
        }
      }
    }
  }

  // calculate PSMs using total loss scoring (no NA-shifted fragments) - used in fast scoring
  static void addPSMsTotalLossScoring_(
    const PeakSpectrum& exp_spectrum,
    const StringView sequence,
    const Size & mod_pep_idx,
    const Size & rna_mod_idx,
    const double & current_peptide_mass,
    const double & current_peptide_mass_without_NA,
    const double & exp_pc_mass,
    const ImmoniumIonsInPeptide & iip, 
    const int & isotope_error,
    const vector<double> & total_loss_template_z1_b_ions, 
    const vector<double> & total_loss_template_z1_y_ions, 
    const boost::math::normal & gaussian_mass_error,
    const double & fragment_mass_tolerance,
    const bool & fragment_mass_tolerance_unit_ppm,
    vector<AnnotatedHit> & annotated_hits,
#ifdef _OPENMP
    omp_lock_t & annotated_hits_lock,
#endif
    const Size& report_top_hits)
  {  
    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    float total_loss_score(0),
      tlss_MIC(0),
      tlss_err(1.0),
      tlss_Morph(0),
      tlss_modds(0),
      pc_MIC(0),
      im_MIC(0);

    size_t n_theoretical_peaks(0);

    vector<double> intensity_sum(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<bool> peak_matched(exp_spectrum.size(), false);

    scorePeptideIons_(
      exp_spectrum, 
      exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
      total_loss_template_z1_b_ions,
      total_loss_template_z1_y_ions,
      current_peptide_mass_without_NA,
      exp_pc_charge,
      iip, 
      fragment_mass_tolerance, 
      fragment_mass_tolerance_unit_ppm,
      intensity_sum,
      b_ions,
      y_ions,
      peak_matched,
      total_loss_score,
      tlss_MIC,
      tlss_Morph,
      tlss_modds,
      tlss_err,
      pc_MIC,
      im_MIC,
      n_theoretical_peaks
    );

    const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));

    //total_loss_score = total_loss_score - 0.22 * (double)sequence.size();

    // early-out if super bad score
    if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_modds, tlss_total_MIC)) { return; }

    const double mass_error_ppm = (current_peptide_mass - exp_pc_mass) / exp_pc_mass * 1e6;
    const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);

    // add peptide hit
    AnnotatedHit ah;
    ah.mass_error_p = mass_error_score;

    ah.sequence = sequence; // copy StringView
    ah.peptide_mod_index = mod_pep_idx;
    ah.total_loss_score = total_loss_score;

    ah.MIC = tlss_MIC;
    ah.err = tlss_err;
    ah.Morph = tlss_Morph;
    ah.modds = tlss_modds;
    ah.immonium_score = im_MIC;
    ah.precursor_score = pc_MIC;

    ah.total_MIC = tlss_total_MIC;

    ah.NA_mod_index = rna_mod_idx;
    ah.isotope_error = isotope_error;
    ah.n_theoretical_peaks = n_theoretical_peaks;
    auto range = make_pair(intensity_sum.begin(), intensity_sum.end());
    ah.ladder_score = ladderScore_(range) / (double)intensity_sum.size(); 
    range = longestCompleteLadder_(intensity_sum.begin(), intensity_sum.end());
    if (range.second != range.first) // see Mascot Percolator paper
    {
      ah.sequence_score = ladderScore_(range) / (double)intensity_sum.size();
    }

    // simple combined score in fast scoring:
    ah.score = calculateFastScore(ah); 

  #ifdef DEBUG_OpenNuXL
    LOG_DEBUG << "best score in pre-score: " << score << endl;
  #endif

  #ifdef _OPENMP
    omp_set_lock(&(annotated_hits_lock));
  #endif
    {
      annotated_hits.emplace_back(move(ah));

      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
      if (annotated_hits.size() >= 2 * report_top_hits)
      {
        std::partial_sort(annotated_hits.begin(), annotated_hits.begin() + report_top_hits, annotated_hits.end(), AnnotatedHit::hasBetterScore);
        annotated_hits.resize(report_top_hits); 
      }
    }
  #ifdef _OPENMP
    omp_unset_lock(&(annotated_hits_lock));
  #endif
  }

  // check for misannotation (absolute m/z instead of offset) and correct
  void checkAndCorrectIsolationWindows_(MSExperiment& e)
  {
    int isolation_windows_reannotated(0);
    int isolation_windows_reannotation_error(0);

    for (MSSpectrum & s : e)
    {
      if (s.getMSLevel() == 2 && s.getPrecursors().size() == 1)
      {
        Precursor& p = s.getPrecursors()[0];
        if (p.getIsolationWindowLowerOffset() > 100.0 && p.getIsolationWindowUpperOffset() > 100.0)
        {
          // in most cases lower and upper offset contain the absolute values.
          // if that is the case we use those
          double left = -(p.getIsolationWindowLowerOffset() - p.getMZ());
          double right = p.getIsolationWindowUpperOffset() - p.getMZ();
          if (left > 0.0 && right > 0.0)
          {
            p.setIsolationWindowLowerOffset(left);
            p.setIsolationWindowUpperOffset(right);
          }
          else // in some files from PD the target m/z is sometimes outside the isolation window (bug?)
          {
            double half_w = (right - left) / 2.0;
            left = p.getMZ() - half_w;
            right = p.getMZ() + half_w;
            p.setIsolationWindowLowerOffset(left);
            p.setIsolationWindowUpperOffset(right);
            isolation_windows_reannotation_error++;
          }
          isolation_windows_reannotated++;          
        }
      }
    }

    if (isolation_windows_reannotated > 0)
    {
      OPENMS_LOG_WARN << "Isolation windows format was incorrect. Reannotated " << isolation_windows_reannotated << " precursors windows. " << endl;
      if (isolation_windows_reannotation_error > 0)
      {
        OPENMS_LOG_WARN << "Reannotation failed for " << isolation_windows_reannotation_error 
          << " precursors windows because the target m/z was outside of boundaries." << endl;
      }
    }
  }

  // returns iterator on start of longest non-zero sequence and end on one-after non-zero sequence (complete ladder)
  template<class Iterator>
  static pair<Iterator, Iterator> longestCompleteLadder_(Iterator b, Iterator e)
  {
    int max_l = 0;
    Iterator best_start(b);
    for (auto i = b; i != e;) // iterate once over vector
    {
      for (; i != e && *i <= 0.0; ++i) {}; // skip zeros
      if (i == e) // end?
      {
        return make_pair(best_start, best_start + max_l);
      }
      unsigned int l = 0;
      Iterator start(i);
      for (; i != e && *i > 0.0; ++i) { ++l; } // count sequence of non-zeros
      if (l > max_l) // longer sequence found?  
      {
        best_start = start;
        max_l = l;
      }
      if (i == e) // end?
      {
        return make_pair(best_start, best_start + max_l);
      }
    }
    return make_pair(best_start, best_start + max_l);
  }

  template<class Iterator>
  static float ladderScore_(pair<Iterator, Iterator> p)
  {
    float MIC(0);
    int count(0);
    for (; p.first != p.second; ++p.first)
    {
      if (*p.first > 0.0)
      {
        MIC += *p.first;
        ++count;
      }
    }
    return count + MIC; // Morph score of matched (complete / partial) ladder
  }

  String convertRawFile_(const String& in, bool no_peak_picking = false)
  {
    writeLog_("RawFileReader reading tool. Copyright 2016 by Thermo Fisher Scientific, Inc. All rights reserved");
    String net_executable = getStringOption_("NET_executable");
    TOPPBase::ExitCodes exit_code;
    QStringList arguments;
    String out = in + ".mzML";
    // check if this file exists and not empty so we can skip further conversions
    if (!File::empty(out)) { return out; }
#ifdef OPENMS_WINDOWSPLATFORM      
    if (net_executable.empty())
    { // default on Windows: if no mono executable is set use the "native" .NET one
      arguments << String("-i=" + in).toQString()
                << String("--output_file=" + out).toQString()
                << String("-f=2").toQString() // indexedMzML
                << String("-e").toQString(); // ignore instrument errors
      if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
      exit_code = runExternalProcess_(getStringOption_("ThermoRaw_executable").toQString(), arguments);
    }
    else
    { // use e.g., mono
      arguments << getStringOption_("ThermoRaw_executable").toQString()
                << String("-i=" + in).toQString()
                << String("--output_file=" + out).toQString()
                << String("-f=2").toQString()
                << String("-e").toQString();
      if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
      exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
    }      
#else
    // default on Mac, Linux: use mono
    net_executable = net_executable.empty() ? "mono" : net_executable;
    arguments << getStringOption_("ThermoRaw_executable").toQString()
              << String("-i=" + in).toQString()
              << String("--output_file=" + out).toQString()
              << String("-f=2").toQString()
              << String("-e").toQString();
    if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
    exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
#endif
    if (exit_code != ExitCodes::EXECUTION_OK)
    {
      OPENMS_LOG_ERROR << "File conversion from RAW file to mzML failed." << endl;
    }
    else
    {
      OPENMS_LOG_INFO << "Raw File successfuly converted to mzML." << endl;
      OPENMS_LOG_INFO << "Please delete it if not needed anymore." << endl;
    }
    return out;
  }

  // datastructure to store longest tag in unshifte/shifted sequence and tag spanning the XL position
  struct XLTags
  {
    size_t tag_unshifted = 0;
    size_t tag_shifted = 0;
    size_t tag_XLed = 0;  // tag that contains the transition from unshifted to shifted
  };


  XLTags getLongestABYLadderWithShift(
       const vector<double>& ab, 
       const vector<double>& y, 
       const vector<double>& ab_xl, 
       const vector<double>& y_xl)
  {
    OPENMS_PRECONDITION(ab.size() == y.size(), "b and y ion arrays must have same size");
    OPENMS_PRECONDITION(ab_xl.size() == y_xl.size(), "cross-linked b and y ion arrays must have same size");

    XLTags tags;

    const size_t n = ab.size();

    // calculate longest consecutive unshifted / shifted sequence and longest sequence spanning unshifted + shifted residues
    vector<int> runAB(n, 0);
    size_t run(0);
    size_t max_ab_run(0);
    for (int l = 0; l != n; ++l)
    {
      if (ab[l] == 0) { run = 0; continue; }
      ++run;
      runAB[l] = run;
      if (run > max_ab_run) max_ab_run = run;
    }
    // runAB[i] now contains current run length e.g.: 000123400100 for prefix ions

    vector<int> runY(n, 0);
    run = 0;
    size_t max_y_run(0);
    for (int l = (int)n - 1; l >= 0; --l)
    {
      if (y[l] == 0) { run = 0; continue; }
      ++run;
      runY[l] = run;
      if (run > max_y_run) max_y_run = run;
    }
    // runY[i] now contains current run length e.g.: 000432100100 for suffix ions

    tags.tag_unshifted = std::max(max_ab_run, max_y_run);

    const size_t n_xl = ab_xl.size();
    if (n_xl != 0)
    {
      OPENMS_PRECONDITION(n_xl == n, "xl and non-xl arrays need to have same size");

      // for XL we calculate the runs in reverse order so we can later quickly calculate the maximum run
      // through non-cross-linked and cross-linked ions

      vector<int> runAB_XL(n_xl, 0);
      run = 0;
      size_t max_ab_shifted(0);
      for (int x = (int)n_xl - 1; x >= 0; --x) // note the reverse order
      {
        if (ab_xl[x] == 0) { run = 0; continue; }
        ++run;
        runAB_XL[x] = run;
        if (run > max_ab_shifted) max_ab_shifted = run;
      }
      // max_ab_shifted contains longest run of shifted ab ions
      // runAB_XL[i] now contains the longest run in X starting at position i e.g.: 00003210000 for prefix ions
    
      vector<int> runY_XL(n_xl, 0);
      run = 0;
      size_t max_y_shifted(0);
      for (int x = 0; x != n_xl; ++x)
      {
        if (y_xl[x] == 0) { run = 0; continue; }
        ++run;
        runY_XL[x] = run;
        if (run > max_y_shifted) max_y_shifted = run;
      }
      // runY_XL[i] now contains the longest run in X starting at position i e.g.: 00001230000 for suffix ions

      tags.tag_shifted = std::max(max_ab_shifted, max_y_shifted);

      size_t maximum_ab_tag_length(0);

      // calculate maximum tag that spans linear intensities and at least one XLed amino acid for prefix ions
      for (Size i = 0; i < n_xl - 1; ++i)
      {
        if (runAB[i] == 0 || runAB_XL[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = runAB[i] + runAB_XL[i + 1]; // tag length if cross-link is introduced at amino acid i+1
        if (tag_length > maximum_ab_tag_length) maximum_ab_tag_length = tag_length; 
      }

      size_t maximum_y_tag_length(0);

      // same for suffix ions
      for (Size i = 0; i < n_xl - 1; ++i)
      {
        if (runY_XL[i] == 0 || runY[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = runY_XL[i] + runY[i + 1]; // tag length with cross-linked part and non-cross-linked
        if (tag_length > maximum_y_tag_length) maximum_y_tag_length = tag_length; 
      }
      tags.tag_XLed = std::max(maximum_ab_tag_length, maximum_y_tag_length);
    }
     
    return tags;
  }

  XLTags getLongestLadderWithShift(const vector<double>& intL, const vector<double>& intXL)
  {
    // calculate longest consecutive unshifted / shifted sequence and longest sequence spanning unshifted + shifted residues
    XLTags tags;             
    vector<int> prefixRunL(intL.size(), 0);
    size_t run(0);
    for (int l = 0; l != intL.size(); ++l)
    {
      if (intL[l] == 0) { run = 0; continue; }
      ++run;
      prefixRunL[l] = run;
      if (run > tags.tag_unshifted) tags.tag_unshifted = run;
    }
    // tags.tag_unshifted contains longest run
    // prefixRunL[i] now contains current run length e.g.: 000123400100 for prefix ions

    vector<int> suffixRunL(intL.size(), 0);
    run = 0;
    for (int l = (int)intL.size() - 1; l >= 0; --l)
    {
      if (intL[l] == 0) { run = 0; continue; }
      ++run;
      suffixRunL[l] = run;
    }
    // suffixRunL[i] now contains current run length e.g.: 000432100100 for suffix ions

    if (!intXL.empty())
    {
      // for XL we calculate the runs in reverse order so we can later quickly calculate the maximum run
      // through non-cross-linked and cross-linked ions

      vector<int> prefixRunX(intXL.size(), 0);
      run = 0;
      for (int x = (int)intXL.size() - 1; x >= 0; --x) // note the reverse order
      {
        if (intXL[x] == 0) { run = 0; continue; }
        ++run;
        prefixRunX[x] = run;
        if (run > tags.tag_shifted) tags.tag_shifted = run;
      }
      // tags.tag_shifted contains longest run
      // prefixRunX[i] now contains the longest run in X starting at position i e.g.: 00003210000 for prefix ions
    
      vector<int> suffixRunX(intXL.size(), 0);
      run = 0;
      for (int x = 0; x != (int)intXL.size(); ++x)
      {
        if (intXL[x] == 0) { run = 0; continue; }
        ++run;
        suffixRunX[x] = run;
      }
      // suffixRunX[i] now contains the longest run in X starting at position i e.g.: 00001230000 for suffix ions

      size_t maximum_tag_length(0);

      // calculate maximum tag that spans linear intensities and at least one XLed amino acid for prefix ions
      for (Size i = 0; i < intXL.size() - 1; ++i)
      {
        if ( prefixRunL[i] == 0 || prefixRunX[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = prefixRunL[i] + prefixRunX[i + 1]; // tag length if cross-link is introduced at amino acid i+1
        if (tag_length > maximum_tag_length) maximum_tag_length = tag_length; 
      }

      // same for suffix ions
      for (Size i = 0; i < intXL.size() - 1; ++i)
      {
        if (suffixRunX[i] == 0 || suffixRunL[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = suffixRunX[i] + suffixRunL[i + 1]; // tag length with cross-linked part and non-cross-linked
        if (tag_length > maximum_tag_length) maximum_tag_length = tag_length; 
      }
      tags.tag_XLed = maximum_tag_length;
    }
    return tags;
  }

  ExitCodes correctPrecursors(MSExperiment& ms_centroided)
  {
    //-------------------------------------------------------------
    // HighRes Precursor Mass Correction
    //-------------------------------------------------------------
    std::vector<double> deltaMZs, mzs, rts;
    std::set<Size> corrected_to_highest_intensity_peak = PrecursorCorrection::correctToHighestIntensityMS1Peak(
      ms_centroided, 
      0.01, // check if we can estimate this from data (here it is given in m/z not ppm)
      false, // is ppm = false
      deltaMZs, 
      mzs, 
      rts
      );      
    writeLog_("Info: Corrected " + String(corrected_to_highest_intensity_peak.size()) + " precursors.");
    if (!deltaMZs.empty())
    {
      vector<double> deltaMZs_ppm, deltaMZs_ppmabs;
      for (Size i = 0; i != deltaMZs.size(); ++i)
      {
        deltaMZs_ppm.push_back(Math::getPPM(mzs[i], mzs[i] + deltaMZs[i]));
        deltaMZs_ppmabs.push_back(Math::getPPMAbs(mzs[i], mzs[i] + deltaMZs[i]));
      }

      double median = Math::median(deltaMZs_ppm.begin(), deltaMZs_ppm.end());
      double MAD =  Math::MAD(deltaMZs_ppm.begin(), deltaMZs_ppm.end(), median);
      double median_abs = Math::median(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end());
      double MAD_abs = Math::MAD(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end(), median_abs);
      writeLog_("Precursor correction to highest intensity peak:\n  median delta m/z  = " 
        + String(median) + " ppm  MAD = " + String(MAD)
        + "\n  median delta m/z (abs.) = " + String(median_abs) 
        + " ppm  MAD = " + String(MAD_abs));
    }

      FeatureMap features;    
    {
      MSExperiment e(ms_centroided); // FFM seems to delete passed spectra
      FeatureFinderMultiplexAlgorithm algorithm;
      Param p = algorithm.getParameters();
      p.setValue("algorithm:labels", ""); // label-free
      p.setValue("algorithm:charge", "2:5");
      p.setValue("algorithm:rt_typical", 30.0);
      p.setValue("algorithm:rt_band", 3.0); // max 3 seconds shifts between isotopic traces
      p.setValue("algorithm:rt_min", 4.0);
      p.setValue("algorithm:spectrum_type", "centroid");
      algorithm.setParameters(p);
      algorithm.run(e, true);
      features = algorithm.getFeatureMap(); 
      writeLog_("Detected peptides: " + String(features.size()));
    }

    set<Size> correct_to_nearest_feature = PrecursorCorrection::correctToNearestFeature(
      features, 
      ms_centroided, 
      20.0, 
      0.01, 
      false, 
      true, 
      false, 
      false, 
      3, 
      10);
    writeLog_("Precursor correction to feature:\n  succesful in = " 
      + String(correct_to_nearest_feature.size()) + " cases.");

    return EXECUTION_OK;
  }


  void mapAccessionToTDProteins(ProteinIdentification& prot_id, map<String, ProteinHit*>& acc2protein_targets, map<String, ProteinHit*>& acc2protein_decoys)
  {
    vector<ProteinHit>& proteins = prot_id.getHits();
    for (ProteinHit& protein : proteins)
    {
      if (protein.getMetaValue("target_decoy").toString().hasPrefix("target"))
      {
        acc2protein_targets[protein.getAccession()] = &protein;
      }
      else
      {
        acc2protein_decoys[protein.getAccession()] = &protein;
      }
    }
  }

  void annotateProteinModificationForTopHits_(vector<ProteinIdentification>& prot_ids, const vector<PeptideIdentification>& peps, TextFile& tsv_file, bool report_decoys = false)
  {
    assert(prots.size() == 1); // support for one run only

    // protein identification run
    ProteinIdentification& prot_id = prot_ids[0];

    // create lookup accession -> protein
    map<String, ProteinHit*> acc2protein_targets;
    map<String, ProteinHit*> acc2protein_decoys;
    mapAccessionToTDProteins(prot_id, acc2protein_targets, acc2protein_decoys);

    struct ModifiedRegion
    {
      ResidueModification* xl;
      int first; 
      int last;
      bool operator<(const ModifiedRegion& rhs) const
      {
        return std::tie(first, last, *xl) < std::tie(rhs.first, rhs.last, *rhs.xl);
      }
    };

    // protein -> regions with XLs without AA level localization and in how many PSMs it was detected
    map<String, map<ModifiedRegion, size_t>> modified_region_xls_targets;    
    map<String, map<ModifiedRegion, size_t>> modified_region_xls_decoys;

    map<String, ResidueModification*> name2mod; // used to free temporary residues

    // store modification statistic for every protein    
    for (const PeptideIdentification& pep : peps)
    {
      auto& hits = pep.getHits();
      if (hits.empty()) continue;
      const PeptideHit& ph = hits[0]; // only consider top hit
      const std::vector<PeptideEvidence>& ph_evidences = ph.getPeptideEvidences();
      const int best_localization = ph.getMetaValue("NuXL:best_localization_position");

      // create a user defined modification (at most once)
      ResidueModification* xl;
      const String NA = ph.getMetaValue("NuXL:NA");
      auto it = name2mod.find(NA);
      if (it != name2mod.end())
      {
        xl = it->second;
      }
      else
      {
        xl = new ResidueModification();
        xl->setFullId(NA);
        xl->setOrigin('X'); // any AA
        name2mod[NA] = xl;
      }

      for (auto& ph_evidence : ph_evidences)
      {
        const String& acc = ph_evidence.getProteinAccession();
        const int peptide_start_in_protein = ph_evidence.getStart();

        if (peptide_start_in_protein < 0) continue;
 
        const int xl_pos_in_protein = peptide_start_in_protein + best_localization;

        bool is_target = acc2protein_targets.find(acc) != acc2protein_targets.end();

        // not localized? annotate region
        if (best_localization < 0)
        {
          ModifiedRegion mr;
          mr.xl = xl;
          mr.first = ph_evidence.getStart();
          mr.last = ph_evidence.getEnd();
          if (is_target)
          {
            modified_region_xls_targets[acc][mr]++;
          }
          else
          {
            modified_region_xls_decoys[acc][mr]++;
          }
        }

        // retrieve protein the evidence points to
        ProteinHit* protein = is_target ? acc2protein_targets[acc] : acc2protein_decoys[acc];
        if (xl_pos_in_protein < protein->getSequence().size())
        {
          auto mods = protein->getModifications();  // TODO: add mutable reference access
          mods.AALevelSummary[xl_pos_in_protein][xl].count++;
          protein->setModifications(mods);
        }
      }
    }
   
    map<char, size_t> aa2protein_count_targets; // how many proteins have that AA cross-linked
    map<char, size_t> aa2psm_count_targets; // how many PSMs have that AA cross-linked
    map<char, size_t> aa2protein_count_decoys; // how many proteins have that AA cross-linked
    map<char, size_t> aa2psm_count_decoys; // how many PSMs have that AA cross-linked

    map<char, double> aa2background_freq; // AA background distribution for normalization
 
    for (const ProteinIdentification& prot_id : prot_ids)
    { 
      const vector<ProteinHit>& phs = prot_id.getHits();
      for (const ProteinHit& protein : phs)
      { // for all identified proteins
        const String& acc = protein.getAccession();        
        bool is_target = acc2protein_targets.find(acc) != acc2protein_targets.end();

        if (!is_target && !report_decoys) continue; // skip decoys if option is set

        const String& seq = protein.getSequence();

        for (const char c : seq) { aa2background_freq[c] += 1.0; }
        tsv_file.addLine(acc);

        auto mods = protein.getModifications();

        // count how many PSMs support XL-ed position in the current protein
        map<size_t, size_t> position2psm_count;
        for (const auto& p2ms : mods.AALevelSummary)
        {
          size_t position = p2ms.first;
          for (const auto& m2s : p2ms.second)
          { 
            const auto stat = m2s.second;
            position2psm_count[position] += stat.count;            
          }
        }

        // put string with AA, position, PSM count in sequence: e.g.: [Y123,2]
        String annotated_sequence;
        size_t p{0};

        if (!position2psm_count.empty())
        {       
          for (const auto p2psm : position2psm_count)
          {
            while (p < p2psm.first) { annotated_sequence += seq[p]; ++p; }
            // p now points to the modified AA
            annotated_sequence += String("[") + seq[p] + String(p2psm.first + 1) + "," + String(p2psm.second) + "]";
            ++p;
          }
          while (p < seq.size()) { annotated_sequence += seq[p]; ++p; } // output AAs after last modification

          tsv_file.addLine(annotated_sequence);
        }
        else
        {
          tsv_file.addLine(seq);
        }

        // output modification AA and count for this protein
        tsv_file.addLine("Cross-link localizations:");
        set<char> already_counted_at_protein_level;
        for (const auto& p2ms : mods.AALevelSummary)
        {
          size_t position = p2ms.first;

          // retrieve AA a xl position
          char AA_at_position = '?';
          if (seq.size() > position)
          {
            AA_at_position = seq[position];
          }

          // count a modified AA only once per protein for protein count
          if (already_counted_at_protein_level.count(AA_at_position) == 0)
          {        
            if (is_target)
            {
              aa2protein_count_targets[AA_at_position]++;
            }
            else
            {
              aa2protein_count_decoys[AA_at_position]++;
            }
            already_counted_at_protein_level.insert(AA_at_position);
          }

          for (const auto& m2s : p2ms.second)
          {
            auto stat = m2s.second;
            size_t count = stat.count;
            if (is_target)
            {
              aa2psm_count_targets[AA_at_position] += count;
            }
            else
            {
              aa2psm_count_decoys[AA_at_position] += count;
            }
            tsv_file.addLine(m2s.first->getFullId() + ":" + AA_at_position + String(position + 1) + "(" + String(count) + ")");
          }
        }
        
        // output list of modified regions (without AA level localization)
        tsv_file.addLine("Cross-linked peptides without localization at single AA-level:");
        for (const auto& p2xlregions : modified_region_xls_targets)
        {
          for (const auto& region : p2xlregions.second)
          {
            const ModifiedRegion& r = region.first;
            const size_t count = region.second;
            tsv_file.addLine(r.xl->getFullId() + ":" + String(r.first + 1) + "-" + String(r.last + 1) + "(" + String(count) + ")");
          }
        }
      }  
    }
 
    double sum{};
    for (const auto& m : aa2background_freq) { sum += m.second; }
    for (auto& m : aa2background_freq) { m.second /= sum; }

    if (report_decoys)
    {
      tsv_file.addLine("=============================================================");
      tsv_file.addLine("Protein summary for decoy proteins:");
      tsv_file.addLine("AA:proteins");
      for (const auto& a2p : aa2protein_count_decoys)
      {
        tsv_file.addLine(String(a2p.first)  + " : " + String(a2p.second));
      }

      tsv_file.addLine("=============================================================");
      tsv_file.addLine("PSM summary for decoy PSMs:");
      tsv_file.addLine("AA:PSMs");
      for (const auto& a2p : aa2psm_count_decoys)
      {
        tsv_file.addLine(String(a2p.first)  + " : " + String(a2p.second));
      }
    }

    tsv_file.addLine("=============================================================");
    tsv_file.addLine("Protein summary for target proteins:");
    tsv_file.addLine("AA:proteins");
    for (const auto& a2p : aa2protein_count_targets)
    {
      tsv_file.addLine(String(a2p.first)  + " : " + String(a2p.second) + "\tfreq. normalized: " + String(a2p.second * aa2background_freq[a2p.first]));
    }

    tsv_file.addLine("=============================================================");
    tsv_file.addLine("PSM summary for target PSMs:");
    tsv_file.addLine("AA:PSMs");
    for (const auto& a2p : aa2psm_count_targets)
    {
      tsv_file.addLine(String(a2p.first)  + " : " + String(a2p.second) + "\tfreq. normalized: " + String(a2p.second * aa2background_freq[a2p.first]));
    }
    
    for (auto m : name2mod) { delete(m.second); } // free memory
  }

  ExitCodes main_(int, const char**) override
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    // Parameter: Input
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(getStringOption_("in"));

    String in_mzml; 
    if (in_type == FileTypes::MZML)
    {
      in_mzml = getStringOption_("in");
    }
    else if (in_type == FileTypes::RAW)
    {
      in_mzml = convertRawFile_(getStringOption_("in"));
    }

    String out_idxml = getStringOption_("out");
    String in_db = getStringOption_("database");

    //
    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    bool generate_decoys = getStringOption_("RNPxl:decoys") == "true";

    StringList filter = getStringList_("filter");
    bool filter_pc_mass_error = find(filter.begin(), filter.end(), "filter_pc_mass_error") != filter.end();
    bool impute_decoy_medians = find(filter.begin(), filter.end(), "impute_decoy_medians") != filter.end();
    bool filter_bad_partial_loss_scores = find(filter.begin(), filter.end(), "filter_bad_partial_loss_scores") != filter.end();
    bool autotune = find(filter.begin(), filter.end(), "autotune") != filter.end();
    bool idfilter = find(filter.begin(), filter.end(), "idfilter") != filter.end();
    bool spectrumclusterfilter = find(filter.begin(), filter.end(), "spectrumclusterfilter") != filter.end();
    bool pcrecalibration = find(filter.begin(), filter.end(), "pcrecalibration") != filter.end();


    if (pcrecalibration) 
    {
      MSExperiment e;
      MzMLFile().load(in_mzml, e);
      correctPrecursors(e);
      in_mzml = FileHandler::stripExtension(in_mzml) + "_pc.mzML";
      OPENMS_LOG_INFO << "Writing calibrated file to: " << in_mzml << endl;
      MzMLFile().store(in_mzml, e);
    }

    InternalCalibration ic; // only filled if pcrecalibration is set and there are enough calibrants


    // autotune (only works if non-XL peptides present)
    set<String> skip_peptide_spectrum;
    double global_fragment_error(0);

    if (autotune || idfilter)
    {
      SimpleSearchEngineAlgorithm sse;
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      Param p = sse.getParameters();
      p.setValue("precursor:mass_tolerance", precursor_mass_tolerance);
      p.setValue("precursor:mass_tolerance_unit", getStringOption_("precursor:mass_tolerance_unit"));
      p.setValue("fragment:mass_tolerance", fragment_mass_tolerance);
      p.setValue("fragment:mass_tolerance_unit", getStringOption_("fragment:mass_tolerance_unit"));
      StringList var_mods = getStringList_("modifications:variable");
      if (find(var_mods.begin(), var_mods.end(), "Phospho (S)") == var_mods.end()) { var_mods.push_back("Phospho (S)"); }
      if (find(var_mods.begin(), var_mods.end(), "Phospho (T)") == var_mods.end()) { var_mods.push_back("Phospho (T)"); }
      if (find(var_mods.begin(), var_mods.end(), "Phospho (Y)") == var_mods.end()) { var_mods.push_back("Phospho (Y)"); }
      if (find(var_mods.begin(), var_mods.end(), "Oxidation (M)") == var_mods.end()) { var_mods.push_back("Oxidation (M)"); }
      StringList fixed_mods = getStringList_("modifications:fixed");
      p.setValue("modifications:fixed", fixed_mods);
      p.setValue("modifications:variable", var_mods);
      p.setValue("modifications:variable_max_per_peptide", 2);
      p.setValue("peptide:missed_cleavages", 2);
      p.setValue("precursor:isotopes", IntList{0, 1});
      p.setValue("decoys", generate_decoys ? "true" : "false");
      p.setValue("enzyme", getStringOption_("peptide:enzyme"));
      p.setValue("annotate:PSM", 
        StringList{
          Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM, 
          Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
          Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
          Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION
        });
      sse.setParameters(p);
      cout << "Running autotune..." << endl;
      sse.search(in_mzml, in_db, prot_ids, pep_ids);      

/// try to run percolator
      {    
        vector<ProteinIdentification> perc_prot_ids;
        vector<PeptideIdentification> perc_pep_ids;

        const String percolator_executable = getStringOption_("percolator_executable");
        bool sufficient_PSMs_for_score_recalibration = pep_ids.size() > 1000;
        if (!percolator_executable.empty() && sufficient_PSMs_for_score_recalibration) // only try to call percolator if we have some PSMs
        {
          IdXMLFile().store(out_idxml, prot_ids, pep_ids);

          // run percolator on idXML
          String perc_out = out_idxml;
          perc_out.substitute(".idXML", "_sse_perc.idXML");
           
          String weights_out = out_idxml;
          weights_out.substitute(".idXML", "_sse.weights");

          QStringList process_params;
          process_params << "-in" << out_idxml.toQString()
                       << "-out" << perc_out.toQString()
                       << "-percolator_executable" << percolator_executable.toQString()
                       << "-train_best_positive" 
                       << "-score_type" << "q-value"
                       << "-post_processing_tdc"
                       << "-weights" << weights_out.toQString()
                       << "-nested_xval_bins" << "3";

          if (getStringOption_("peptide:enzyme") == "Lys-C")
          {
            process_params << "-enzyme" << "lys-c";
          }
                       
          TOPPBase::ExitCodes exit_code = runExternalProcess_(QString("PercolatorAdapter"), process_params);

          if (exit_code != EXECUTION_OK) 
          { 
            OPENMS_LOG_WARN << "Score recalibration failed in IDFilter. Using original results." << endl; 
          }
          else
          { 
            // load back idXML
            IdXMLFile().load(perc_out, perc_prot_ids, perc_pep_ids);
 
            // generate filtered results
            IDFilter::keepNBestHits(perc_pep_ids, 1);
            IDFilter::removeUnreferencedProteins(perc_prot_ids, perc_pep_ids);
          }
        }

        cout << "Filtering ..." << endl;
        IDFilter::filterHitsByScore(perc_pep_ids, 0.01); // 1% PSM-FDR
        IDFilter::removeEmptyIdentifications(perc_pep_ids);
        cout << "Peptide PSMs at 1% FDR: " << perc_pep_ids.size() << endl;

        // ID-filter part for linear peptides
        if (idfilter)
        {
          for (const auto& pi : perc_pep_ids)
          {
            skip_peptide_spectrum.insert((String)pi.getMetaValue("spectrum_reference")); // get native id
          }
        }

        if (spectrumclusterfilter)
        {
          Size skipped_similar_spectra(0);
          // load MS2 map
          PeakMap spectra;
          MzMLFile f;
          f.setLogType(log_type_);
          PeakFileOptions options;
          options.clearMSLevels();
          options.addMSLevel(2);
          f.getOptions() = options;
          f.load(in_mzml, spectra);
          spectra.sortSpectra(true);
          SpectrumLookup lookup;
          lookup.readSpectra(spectra);
          // build kdtree
          Param p;
          p.setValue("rt_tol", 60.0);
          p.setValue("mz_tol", precursor_mass_tolerance);
          p.setValue("mz_unit", "ppm");
          FeatureMap fmap;
          for (Size i = 0; i != spectra.size(); ++i)
          {
            const MSSpectrum& s = spectra[i];
            
            Feature feat;
            feat.setMZ(s.getPrecursors()[0].getMZ());
            feat.setRT(s.getRT());
            feat.setMetaValue("native_id", s.getNativeID());
            fmap.push_back(feat);
          }
          vector<FeatureMap> fmaps;
          fmaps.push_back(std::move(fmap));
          KDTreeFeatureMaps kdtree(fmaps, p);

          // filter all coeluting MS2 with high spectral similarity to identified one
          for (const auto& pi : perc_pep_ids)
          {
            String this_native_id = (String)pi.getMetaValue("spectrum_reference");

            std::vector<Size> result_indices;

            // find neighbors
            double m = Math::ppmToMass(precursor_mass_tolerance, pi.getMZ());
            kdtree.queryRegion(pi.getRT() - 60.0, pi.getRT() + 60.0, pi.getMZ() - m, pi.getMZ() + m, result_indices);
            
            if (result_indices.size() > 1)
            {
              for (Size ix : result_indices)
              {
                auto f = kdtree.feature(ix);
                const String other_native_id = f->getMetaValue("native_id");

                // skip self-comparison and already identified spectra
                if (this_native_id == other_native_id || skip_peptide_spectrum.count(other_native_id) > 0) continue;

                const MSSpectrum& this_spec = spectra[lookup.findByNativeID(this_native_id)];
                const MSSpectrum& other_spec = spectra[lookup.findByNativeID(other_native_id)];
                BinnedSpectrum bs1 (this_spec, BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
                BinnedSpectrum bs2 (other_spec, BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES); 
                const float contrast_angle = BinnedSpectralContrastAngle()(bs1, bs2);
                if (contrast_angle > 0.9) 
                {
                  skip_peptide_spectrum.insert(other_native_id);
                  skipped_similar_spectra++;
                }
              }
            }
          }
          cout << "Excluded coelution precursors with high spectral similarity: " << skipped_similar_spectra << endl;           
        }
      }
////////// end percolator part

      cout << "Calculating FDR..." << endl;
      FalseDiscoveryRate fdr;
      fdr.apply(pep_ids);
      cout << "Filtering ..." << endl;
      IDFilter::filterHitsByScore(pep_ids, 0.01); // 1% PSM-FDR
      IDFilter::removeEmptyIdentifications(pep_ids);
      cout << "Peptide PSMs at 1% FDR (no percolator): " << pep_ids.size() << endl;
 
      if (pep_ids.size() > 100)
      {
        vector<double> median_fragment_error_ppm_abs;
        vector<double> median_fragment_error_ppm;
        vector<double> precursor_error_ppm;
        for (const auto& pi : pep_ids)
        {
          const PeptideHit& ph = pi.getHits()[0];
          if (ph.metaValueExists(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM))
          {
            //cout << ph.getMetaValue("median_fragment_error_ppm") << endl;
            double fragment_error = (double)ph.getMetaValue(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
            median_fragment_error_ppm_abs.push_back(fabs(fragment_error));
            median_fragment_error_ppm.push_back(fragment_error);
          }
          if (ph.metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
          {
            precursor_error_ppm.push_back((double)ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM));
          }
        }
        sort(median_fragment_error_ppm_abs.begin(), median_fragment_error_ppm_abs.end());
        sort(median_fragment_error_ppm.begin(), median_fragment_error_ppm.end());
        sort(precursor_error_ppm.begin(), precursor_error_ppm.end());

        // use 68-percentile as in identipy
        double new_fragment_mass_tolerance = 4.0 * median_fragment_error_ppm_abs[median_fragment_error_ppm_abs.size() * 0.68]; 
        global_fragment_error = median_fragment_error_ppm[median_fragment_error_ppm_abs.size() * 0.5]; // median of all fragment errors
        double left_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.005];
        double median_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.5];
        double right_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.995];

        if (autotune)
        {
          fragment_mass_tolerance = new_fragment_mass_tolerance; // set new fragment mass tolerance
        }
        cout << "New fragment mass tolerance (ppm): " << new_fragment_mass_tolerance << endl;
        cout << "Global fragment mass shift (ppm): " << global_fragment_error << endl;
        cout << "Estimated precursor mass tolerance (ppm): " << left_precursor_mass_tolerance << "\t" << median_precursor_mass_tolerance << "\t" << right_precursor_mass_tolerance << endl;
      }
      else
      {
         OPENMS_LOG_INFO << "autotune: too few non-cross-linked peptides found. Will keep parameters as-is." << endl;
      }

      if (pcrecalibration)
      {
        ic.setLogType(log_type_);
        ic.fillCalibrants(pep_ids, precursor_mass_tolerance);
        if (global_fragment_error != 0)
        {
          PeakMap spectra;
          MzMLFile f;
          f.load(in_mzml, spectra);
          spectra.sortSpectra(true);
          for (auto & s : spectra)
          {
            if (s.getMSLevel() != 2) continue;
            for (auto & p : s)
            {
              p.setMZ(p.getMZ() - Math::ppmToMass(global_fragment_error, p.getMZ())); // correct global fragment error
            }
          }
         f.store(in_mzml, spectra);
        }
      }
    }
   
    OPENMS_LOG_INFO << "IDFilter excludes " << skip_peptide_spectrum.size() << " spectra." << endl;
 
    String out_tsv = getStringOption_("out_tsv");

    fast_scoring_ = getStringOption_("RNPxl:scoring") == "fast" ? true : false;

    // true positives: assumed gaussian distribution of mass error
    // with sigma^2 = precursor_mass_tolerance
    boost::math::normal gaussian_mass_error(0.0, sqrt(precursor_mass_tolerance));

    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    IntList precursor_isotopes = getIntList_("precursor:isotopes");

    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    double marker_ions_tolerance = getDoubleOption_("RNPxl:marker_ions_tolerance");

    double small_peptide_mass_filter_threshold = getDoubleOption_("RNPxl:filter_small_peptide_mass");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    if (fixed_unique.size() != fixedModNames.size())
    {
      OPENMS_LOG_WARN << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      OPENMS_LOG_WARN << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(fixedModNames);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    size_t report_top_hits = (size_t)getIntOption_("report:top_hits");
    double peptide_FDR = getDoubleOption_("report:peptideFDR");
    DoubleList XL_FDR = getDoubleList_("report:xlFDR");

    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides = getStringList_("RNPxl:target_nucleotides");

    StringList nt_groups = getStringList_("RNPxl:nt_groups");

    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings = getStringList_("RNPxl:mapping");

    // read list of nucleotides that can directly cross-link
    // these are responsible for shifted fragment ions. Their fragment adducts thus determine which shifts will be observed on b-,a-,y-ions
    const String& can_cross_link = getStringOption_("RNPxl:can_cross_link");
    for (const auto& c : can_cross_link) { can_xl_.insert(c); }

    StringList modifications = getStringList_("RNPxl:modifications");

    String sequence_restriction = getStringOption_("RNPxl:sequence");

    Int max_nucleotide_length = getIntOption_("RNPxl:length");

    bool cysteine_adduct = getFlag_("RNPxl:CysteineAdduct");

    // generate mapping from empirical formula to mass and empirical formula to (one or more) precursor adducts
    RNPxlModificationMassesResult mm;
    if (max_nucleotide_length != 0)
    {
      mm = RNPxlModificationsGenerator::initModificationMassesNA(
            target_nucleotides,
            nt_groups, 
            can_xl_,
            mappings,
            modifications, 
            sequence_restriction, 
            cysteine_adduct, 
            max_nucleotide_length);
    }

    if (!getFlag_("RNPxl:only_xl"))
    {
      mm.formula2mass[""] = 0; // insert "null" modification otherwise peptides without NA will not be searched
      mm.mod_combinations[""].insert("none");
    }

    // parse tool parameter and generate all fragment adducts

    // first, we determine which fragments adducts can be generated from a single nucleotide (that has no losses)
    RNPxlParameterParsing::NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts = RNPxlParameterParsing::getTargetNucleotideToFragmentAdducts(getStringList_("RNPxl:fragment_adducts"));

    // calculate all feasible fragment adducts from all possible precursor adducts
    RNPxlParameterParsing::PrecursorsToMS2Adducts all_feasible_fragment_adducts = RNPxlParameterParsing::getAllFeasibleFragmentAdducts(mm, nucleotide_to_fragment_adducts, can_xl_);

    // calculate FDR
    FalseDiscoveryRate fdr;
    Param p = fdr.getParameters();
    p.setValue("add_decoy_proteins", "true"); // we still want decoys in the result (e.g., to run percolator)
    p.setValue("add_decoy_peptides", "true");
    if (report_top_hits >= 2)
    {
      p.setValue("use_all_hits", "true");
    }
    fdr.setParameters(p);

    // load MS2 map
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    // load both MS1 and MS2 for precursor purity annotation
    map<String, PrecursorPurity::PurityScores> purities;
    {
      PeakMap tmp_spectra;
      f.load(in_mzml, tmp_spectra);
      int nMS1 = std::count_if(tmp_spectra.begin(), tmp_spectra.end(), [](MSSpectrum& s){return s.getMSLevel() == 1;});
      if (nMS1 != 0)
      {
        // if isolation windows are properly annotated and correct if necessary
        checkAndCorrectIsolationWindows_(tmp_spectra);
        purities = PrecursorPurity::computePrecursorPurities(tmp_spectra, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
      }
    } // free spectra  

    /////////////////////////////////////////////////////////////////////////
    // define percolator feature set
    /* default features added in PercolatorAdapter:
     * SpecId, ScanNr, ExpMass, CalcMass, mass, 
     * peplen, charge#min..#max, enzN, enzC, enzInt, dm, absdm
     */     
    feature_set_
       << "NuXL:mass_error_p"
       << "NuXL:err"
       << "NuXL:total_loss_score"
       << "NuXL:modds"
       << "NuXL:immonium_score"
       << "NuXL:precursor_score"
       << "NuXL:MIC"
       << "NuXL:Morph"
       << "NuXL:total_MIC"
       << "NuXL:ladder_score"
       << "NuXL:sequence_score"
       << "NuXL:total_Morph"
       << "NuXL:total_HS"
       << "NuXL:tag_XLed"
       << "NuXL:tag_unshifted"
       << "NuXL:tag_shifted"
       << "NuXL:aminoacid_max_tag"
       << "NuXL:aminoacid_id_to_max_tag_ratio"
       << "nr_candidates"
       << "NuXL:explained_peak_fraction"
       << "NuXL:QQ_TIC"
       << "NuXL:QQ_EXPLAINED_FRACTION"
       << "NuXL:wTop50"

       << "NuXL:marker_ions_score"
       << "NuXL:partial_loss_score"
       << "NuXL:pl_MIC"
       << "NuXL:pl_err"
       << "NuXL:pl_Morph"
       << "NuXL:pl_modds"
       << "NuXL:pl_pc_MIC"
       << "NuXL:pl_im_MIC"

       << "NuXL:isPhospho" 
       << "NuXL:isXL" 
       << "NuXL:score"
       << "isotope_error"
       << "variable_modifications"
       << "precursor_intensity_log10"
       << "NuXL:NA_MASS_z0"
       << "NuXL:NA_length"   
       << "nucleotide_mass_tags"
       << "n_theoretical_peaks";

    if (!purities.empty()) feature_set_ << "precursor_purity";

    // one-hot encoding of cross-linked nucleotide
    for (const auto& c : can_cross_link)  { feature_set_ << String("NuXL:XL_" + String(c)); }

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    // only executed if we have a pre-search with enough calibrants
    if (ic.getCalibrationPoints().size() > 1)
    {
      MZTrafoModel::MODELTYPE md = MZTrafoModel::LINEAR;
      bool use_RANSAC = true;

      Size RANSAC_initial_points = (md == MZTrafoModel::LINEAR) ? 2 : 3;
      Math::RANSACParam p(RANSAC_initial_points, 70, 10, 30, true); // TODO: check defaults (taken from tool)
      MZTrafoModel::setRANSACParams(p);

      // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
      MZTrafoModel::setCoefficientLimits(25.0, 25.0, 0.5); 

      IntList ms_level = {1};
      double rt_chunk = 300.0; // 5 minutes covered by each linear model
      String qc_residual_path, qc_residual_png_path;

      if (!ic.calibrate(spectra, ms_level, md, rt_chunk, use_RANSAC, 
              10.0,
              5.0, 
              "",                      
              "",
              qc_residual_path,
              qc_residual_png_path,
              "Rscript"))
      {
         OPENMS_LOG_WARN << "\nCalibration failed. See error message above!" << std::endl;
      }        
    }

    progresslogger.startProgress(0, 1, "Filtering spectra...");
    const bool convert_to_single_charge = false;  // whether to convert fragment peaks with isotopic patterns to single charge
    const double window_size = getDoubleOption_("window_size");
    const size_t peak_count = getIntOption_("peak_count");
    preprocessSpectra_(spectra, 
                       fragment_mass_tolerance, 
                       fragment_mass_tolerance_unit_ppm, 
                       convert_to_single_charge,
                       true, window_size, peak_count); // annotate charge  
    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Calculate Nucleotide Tags...");
    calculateNucleotideTags_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, nucleotide_to_fragment_adducts);
    progresslogger.endProgress();

    // build multimap of precursor mass to scan index (and perform some mass and length based filtering)
    progresslogger.startProgress(0, 1, "Mapping precursors to scan...");
    using MassToScanMultiMap = multimap<double, pair<Size, int>>;
    MassToScanMultiMap multimap_mass_2_scan_index;  // map precursor mass to scan index and (potential) isotopic missassignment
    mapPrecursorMassesToScans(min_precursor_charge,
                              max_precursor_charge,
                              precursor_isotopes,
                              small_peptide_mass_filter_threshold,
                              peptide_min_size,
                              spectra,
                              multimap_mass_2_scan_index);
    progresslogger.endProgress();

    // preallocate storage for PSMs
    vector<size_t> nr_candidates(spectra.size(), 0);
    vector<vector<AnnotatedHit> > annotated_XLs(spectra.size(), vector<AnnotatedHit>());
    for (auto & a : annotated_XLs) { a.reserve(2 * report_top_hits); }
    vector<vector<AnnotatedHit> > annotated_peptides(spectra.size(), vector<AnnotatedHit>());
    for (auto & a : annotated_peptides) { a.reserve(2 * report_top_hits); }
#ifdef _OPENMP     
    // locking is done at the spectrum level to ensure good parallelisation 
    vector<omp_lock_t> annotated_XLs_lock(annotated_XLs.size());
    for (size_t i = 0; i != annotated_XLs_lock.size(); i++) { omp_init_lock(&(annotated_XLs_lock[i])); }
    vector<omp_lock_t> annotated_peptides_lock(annotated_peptides.size());
    for (size_t i = 0; i != annotated_peptides_lock.size(); i++) { omp_init_lock(&(annotated_peptides_lock[i])); }
#endif

#ifdef ANNOTATED_QUANTILES
    vector<quantile_accu_t> annotated_peptides_quantiles_peptides(spectra.size(), quantile_accu_t(quantile_probability = 0.95));
    vector<quantile_accu_t> annotated_peptides_quantiles_XLs(spectra.size(), quantile_accu_t(quantile_probability = 0.95));

    vector<SpectrumLevelScoreQuantiles> QQ_TIC(spectra.size(), SpectrumLevelScoreQuantiles() );
    vector<SpectrumLevelScoreQuantiles> QQ_EXPLAINED_FRACTION(spectra.size(), SpectrumLevelScoreQuantiles() );
#endif

    // load fasta file
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_db, fasta_db);
    progresslogger.endProgress();

    // generate decoy protein sequences by reversing them
    if (generate_decoys)
    {
      progresslogger.startProgress(0, 1, "Generating decoys...");
      ProteaseDigestion digestor;
      const String enzyme = getStringOption_("peptide:enzyme");
      digestor.setEnzyme(getStringOption_("peptide:enzyme"));
      digestor.setMissedCleavages(0);  // for decoy generation disable missed cleavages

      DecoyGenerator dg;
      dg.setSeed(4711);

      // append decoy proteins
      const size_t old_size = fasta_db.size();
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = fasta_db[i];

        std::vector<AASequence> output;
        digestor.digest(AASequence::fromString(e.sequence), output);

        // pseudo reverse protein digest
        e.sequence = "";
        for (const auto & aas : output)
        {
          if (aas.size() <= 2) { e.sequence += aas.toUnmodifiedString(); continue; }

          //cout << "Seq. before reverse: " << s << endl; 
          //auto last = --s.end();
          //std::reverse(s.begin(), last); // standard peptide reverse
          //cout << "Seq. after reverse : " << s << endl; 

          // switch cleavage sites (like in Andromeda or DecoyPyrate)
          // if (s.size() >= 2) std::swap(s[s.size() - 1], s[s.size() - 2]); 
          e.sequence += dg.shufflePeptides(aas, enzyme).toUnmodifiedString();
        }
        e.identifier = "DECOY_" + e.identifier;
        fasta_db.push_back(e);
      }
      // randomize order of targets and decoys to introduce no global 
      // bias in cases where a target has same score as decoy. (we always take the first best scoring one)
      std::random_shuffle(fasta_db.begin(), fasta_db.end());
      progresslogger.endProgress();
    }


    // set up enzyme
    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    ProteaseDigestion digestor;
    digestor.setEnzyme(getStringOption_("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

    progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    set<StringView> processed_petides;

    // set minimum size of peptide after digestion
    Size min_peptide_length = (Size)getIntOption_("peptide:min_size");
    Size max_peptide_length = (Size)getIntOption_("peptide:max_size");

    Size count_proteins(0), count_peptides(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++count_proteins;

      IF_MASTERTHREAD
      {
        progresslogger.setProgress((SignedSize)count_proteins);
      }

      vector<StringView> current_digest;

      auto const & current_fasta_entry = fasta_db[fasta_index];

      bool is_decoy = current_fasta_entry.identifier.hasPrefix("DECOY_");

      digestor.digestUnmodified(current_fasta_entry.sequence, current_digest, min_peptide_length, max_peptide_length);

      for (auto cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          // skip peptide (and all modified variants) if already processed
          if (processed_petides.find(*cit) != processed_petides.end())
          {
            already_processed = true;
          }
          else
          {
            processed_petides.insert(*cit);
          }
        }

        if (already_processed) { continue; }

#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;

        const String unmodified_sequence = cit->getString();

         // only process peptides without ambiguous amino acids (placeholder / any amino acid)
        if (unmodified_sequence.find_first_of("XBZ") != std::string::npos) continue;

        // determine which residues might give rise to an immonium ion
        ImmoniumIonsInPeptide iip(unmodified_sequence);

        if (is_decoy) iip.swap();     

        AASequence aas = AASequence::fromString(unmodified_sequence);
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);
        
        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass_without_NA = fixed_and_variable_modified_peptide.getMonoWeight();

          //create empty theoretical spectrum.  total_loss_spectrum_z2 contains both charge 1 and charge 2 peaks
          vector<double> total_loss_template_z1_b_ions, total_loss_template_z1_y_ions;

          // spectrum containing additional peaks for sub scoring
          PeakSpectrum marker_ions_sub_score_spectrum;

          // iterate over all NA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          Size NA_mod_index = 0;

          for (std::map<String, double>::const_iterator rna_mod_it = mm.formula2mass.begin(); 
            rna_mod_it != mm.formula2mass.end(); 
            ++rna_mod_it, ++NA_mod_index)
          {            
            const double precursor_rna_mass = rna_mod_it->second;
            const double current_peptide_mass = current_peptide_mass_without_NA + precursor_rna_mass; // add NA mass

            // determine MS2 precursors that match to the current peptide mass
            MassToScanMultiMap::const_iterator low_it, up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + current_peptide_mass * precursor_mass_tolerance * 1e-6);
            }
            else // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + precursor_mass_tolerance);
            }

            if (low_it == up_it) { continue; } // no matching precursor in data

            // add peaks for b- and y- ions with charge 1 (sorted by m/z)
            // total / complete loss spectra are generated for fast and (slow) full scoring
            if (total_loss_template_z1_b_ions.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::ResidueType::BIon, total_loss_template_z1_b_ions);
              generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::ResidueType::YIon, total_loss_template_z1_y_ions);
            }


            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // all ion scoring
            //
            if (!fast_scoring_)
            {
              // retrieve NA adduct name
              auto mod_combinations_it = mm.mod_combinations.begin();
              std::advance(mod_combinations_it, NA_mod_index);
              const String& precursor_rna_adduct = *mod_combinations_it->second.begin(); // TODO: check if sufficient to consider only first NA adduct here ?????????????????????????????????????????????????????????????????????

              if (precursor_rna_adduct == "none")
              {
                // score peptide without NA (same method as fast scoring)
                for (auto & l = low_it; l != up_it; ++l)
                {
                  const Size & scan_index = l->second.first;
                  const PeakSpectrum & exp_spectrum = spectra[scan_index];

#ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                  // require at least one mass tag
                  if (exp_spectrum.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
#endif

                  // count candidate for spectrum
#ifdef _OPENMP
                  omp_set_lock(&(annotated_peptides_lock[scan_index]));
                  omp_set_lock(&(annotated_XLs_lock[scan_index]));
                  ++nr_candidates[scan_index];
                  omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                  omp_unset_lock(&(annotated_peptides_lock[scan_index]));
#endif
                  //const double exp_pc_mass = l->first;
                  const int & isotope_error = l->second.second;
                  const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

                  float total_loss_score(0), 
                        tlss_MIC(0),
                        tlss_err(0), 
                        tlss_Morph(0),
                        tlss_modds(0),
                        pc_MIC(0),
                        im_MIC(0);
                  size_t n_theoretical_peaks(0);

                  vector<double> intensity_linear(total_loss_template_z1_b_ions.size(), 0.0);
                  vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0);
                  vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0);

                  vector<bool> peak_matched(exp_spectrum.size(), false);
                  scorePeptideIons_(
                    exp_spectrum, 
                    exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                    total_loss_template_z1_b_ions,
                    total_loss_template_z1_y_ions,
                    current_peptide_mass_without_NA,
                    exp_pc_charge,
                    iip, 
                    fragment_mass_tolerance, 
                    fragment_mass_tolerance_unit_ppm,
                    intensity_linear,
                    b_ions, 
                    y_ions,
                    peak_matched,
                    total_loss_score,
                    tlss_MIC,
                    tlss_Morph,
                    tlss_modds,
                    tlss_err,
                    pc_MIC,
                    im_MIC,
                    n_theoretical_peaks
                  );                  

                  const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));
    
//                  total_loss_score = total_loss_score - 0.22 * (double)cit->size();

                  if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_modds, tlss_total_MIC)) { continue; }

                  const double mass_error_ppm = (current_peptide_mass - l->first) / l->first * 1e6;
                  const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);
                  
                  // add peptide hit
                  AnnotatedHit ah;
                  ah.mass_error_p = mass_error_score;
                  ah.sequence = *cit; // copy StringView
                  ah.peptide_mod_index = mod_pep_idx;
                  ah.MIC = tlss_MIC;
                  ah.err = tlss_err;
                  ah.Morph = tlss_Morph;
                  ah.modds = tlss_modds;
                  ah.total_loss_score = total_loss_score;
                  ah.immonium_score = im_MIC;
                  ah.precursor_score = pc_MIC;
                  ah.total_MIC = tlss_total_MIC;

                  ah.NA_mod_index = NA_mod_index;
                  ah.isotope_error = isotope_error;

                  auto range = make_pair(intensity_linear.begin(), intensity_linear.end());
                  ah.ladder_score = ladderScore_(range) / (double)intensity_linear.size(); 
                  range = longestCompleteLadder_(intensity_linear.begin(), intensity_linear.end());
                  if (range.second != range.first)
                  {
                    ah.sequence_score = ladderScore_(range) / (double)intensity_linear.size();
                  }

                  RankScores rankscores = rankScores_(exp_spectrum, peak_matched);
                  ah.explained_peak_fraction = rankscores.explained_peak_fraction;
                  ah.wTop50 = rankscores.wTop50;

                  // do we have at least one ladder peak
//                  const XLTags longest_tags = getLongestLadderWithShift(intensity_linear, vector<double>());

                  const XLTags longest_tags = getLongestABYLadderWithShift(b_ions, y_ions, vector<double>(), vector<double>());

#ifdef FILTER_BAD_SCORES_ID_TAGS
                  if (longest_tags.tag_unshifted == 0) continue;
#endif
                  ah.tag_XLed = longest_tags.tag_XLed;
                  ah.tag_unshifted = longest_tags.tag_unshifted;
                  ah.tag_shifted = longest_tags.tag_shifted;

                  // combined score
                  const double tags = exp_spectrum.getFloatDataArrays()[1][0];
                  ah.score = OpenNuXL::calculateCombinedScore(ah, false, tags);
                  ah.n_theoretical_peaks = n_theoretical_peaks;
                  //ah.score = OpenNuXL::calculateFastScore(ah); does this work too

#ifdef DEBUG_OpenNuXL
                  OPENMS_LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP 
                  omp_set_lock(&(annotated_peptides_lock[scan_index]));
#endif
                  {
#ifdef ANNOTATED_QUANTILES
                    annotated_peptides_quantiles_peptides[scan_index](ah.total_loss_score);

                    QQ_TIC[scan_index].insert(ah.total_MIC);
                    QQ_EXPLAINED_FRACTION[scan_index].insert(ah.explained_peak_fraction);
#endif
                    annotated_peptides[scan_index].emplace_back(move(ah));

                    // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                    if (annotated_peptides[scan_index].size() >= 2 * report_top_hits)
                    {
                      std::partial_sort(annotated_peptides[scan_index].begin(), annotated_peptides[scan_index].begin() + report_top_hits, annotated_peptides[scan_index].end(), AnnotatedHit::hasBetterScore);
                      annotated_peptides[scan_index].resize(report_top_hits); 
                    }
                  }
#ifdef _OPENMP 
                  omp_unset_lock(&(annotated_peptides_lock[scan_index]));
#endif
                }
              }
              else  // score peptide with NA MS1 adduct
              {
                // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
                // get NA fragment shifts in the MS2 (based on the precursor RNA/DNA)
                auto const & all_NA_adducts = all_feasible_fragment_adducts.at(precursor_rna_adduct);
                const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_NA_adducts.feasible_adducts;
                // get marker ions
                const vector<RNPxlFragmentAdductDefinition>& marker_ions = all_NA_adducts.marker_ions;

                //cout << "'" << precursor_rna_adduct << "'" << endl;
                //OPENMS_POSTCONDITION(!feasible_MS2_adducts.empty(),
                //                String("FATAL: No feasible adducts for " + precursor_rna_adduct).c_str());


                // Do we have (nucleotide) specific fragmentation adducts? for the current NA adduct on the precursor?
                // If so, generate spectra for shifted ion series


                // score individually for every nucleotide
                for (auto const & nuc_2_adducts : feasible_MS2_adducts)
                {
                  // determine current nucleotide and associated partial losses
                  const char& cross_linked_nucleotide = nuc_2_adducts.first;
                  const vector<RNPxlFragmentAdductDefinition>& partial_loss_modification = nuc_2_adducts.second;

                  PeakSpectrum marker_ions_sub_score_spectrum_z1;
                  // add shifted marker ions of charge 1
                  marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
                  marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation                   
                  RNPxlFragmentIonGenerator::addMS2MarkerIons(
                    marker_ions,
                    marker_ions_sub_score_spectrum_z1,
                    marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[IA_CHARGE_INDEX],
                    marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

                  // nucleotide is associated with certain NA-related fragment losses?
                  vector<double> partial_loss_template_z1_bions, partial_loss_template_z1_yions;
                  if (!partial_loss_modification.empty())
                  {
                    generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::BIon, partial_loss_template_z1_bions);
                    generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::YIon, partial_loss_template_z1_yions);
                  } 

                  for (auto & l = low_it; l != up_it; ++l)
                  {
                    const Size & scan_index = l->second.first;
                    const PeakSpectrum & exp_spectrum = spectra[scan_index];
                    //////////////////////////////////////////
                    //               ID-Filter
                    if (skip_peptide_spectrum.find(exp_spectrum.getNativeID()) != skip_peptide_spectrum.end()) { continue; }
#ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                    // require at least one mass tag
                    if (exp_spectrum.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
#endif

#ifdef _OPENMP
                    omp_set_lock(&(annotated_peptides_lock[scan_index]));
                    omp_set_lock(&(annotated_XLs_lock[scan_index]));
                    // count candidate for spectrum
                    ++nr_candidates[scan_index];
                    omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                    omp_unset_lock(&(annotated_peptides_lock[scan_index]));
#endif
                    const int & isotope_error = l->second.second;
                    float tlss_MIC(0), 
                      tlss_err(1.0), 
                      tlss_Morph(0),
                      tlss_modds(0),
                      partial_loss_sub_score(0), 
                      marker_ions_sub_score(0),
                      total_loss_score(0),
                      pc_MIC(0),
                      im_MIC(0);
                    size_t n_theoretical_peaks(0);

                    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

                    vector<double> intensity_linear(total_loss_template_z1_b_ions.size(), 0.0);                    
                    vector<bool> peak_matched(exp_spectrum.size(), false);
                    vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0);  // b & a ions
                    vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0); 

                    scorePeptideIons_(
                      exp_spectrum, 
                      exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                      total_loss_template_z1_b_ions,
                      total_loss_template_z1_y_ions,
                      current_peptide_mass_without_NA,
                      exp_pc_charge,
                      iip, 
                      fragment_mass_tolerance, 
                      fragment_mass_tolerance_unit_ppm,
                      intensity_linear,
                      b_ions,
                      y_ions,
                      peak_matched,
                      total_loss_score,
                      tlss_MIC,
                      tlss_Morph,
                      tlss_modds,
                      tlss_err,
                      pc_MIC,
                      im_MIC,
                      n_theoretical_peaks
                    );

                    const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));

//                    total_loss_score = total_loss_score - 0.22 * (double)cit->size();

                    if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_modds, tlss_total_MIC)) { continue; }

                    vector<double> intensity_xls(total_loss_template_z1_b_ions.size(), 0.0);

                    vector<double> b_xl_ions(b_ions.size(), 0.0);
                    vector<double> y_xl_ions(b_ions.size(), 0.0);

                    float plss_MIC(0), 
                      plss_err(1.0), 
                      plss_Morph(0), 
                      plss_modds(0),
                      plss_pc_MIC(0),
                      plss_im_MIC(0);

                    scoreXLIons_(partial_loss_modification,
                                 iip,
                                 exp_spectrum,
                                 current_peptide_mass_without_NA,
                                 fragment_mass_tolerance, 
                                 fragment_mass_tolerance_unit_ppm,
                                 partial_loss_template_z1_bions, 
                                 partial_loss_template_z1_yions,
                                 marker_ions_sub_score_spectrum_z1,
                                 intensity_xls,
                                 b_xl_ions,
                                 y_xl_ions,
                                 peak_matched,
                                 partial_loss_sub_score,
                                 marker_ions_sub_score,
                                 plss_MIC, 
                                 plss_err, 
                                 plss_Morph,
                                 plss_modds,
                                 plss_pc_MIC,
                                 plss_im_MIC,
                                 n_theoretical_peaks,
                                 is_decoy);

                    const double total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC)) + plss_MIC + (plss_pc_MIC - floor(plss_pc_MIC)) + plss_im_MIC + marker_ions_sub_score;

                    // decreases number of hits (especially difficult cases - but also number of false discoveries)
                    if (filter_bad_partial_loss_scores && badPartialLossScore(tlss_Morph, plss_Morph, plss_MIC, plss_im_MIC, (plss_pc_MIC - floor(plss_pc_MIC)), marker_ions_sub_score))  
                    { 
                      continue; 
                    }

                    const double mass_error_ppm = (current_peptide_mass - l->first) / l->first * 1e6;
                    const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);
                    
                    // add peptide hit
                    AnnotatedHit ah;
                    ah.mass_error_p = mass_error_score;

                    ah.sequence = *cit; // copy StringView
                    ah.peptide_mod_index = mod_pep_idx;
/*
/////////////////////////////////////////////////////////////////////////////// test recalculate hyperscore on merged XL/non-XL ladders
                    size_t y_ion_count = std::count_if(y_ions.begin(), y_ions.end(), [](double d) { return d > 1e-6; });
                    size_t b_ion_count = std::count_if(b_ions.begin(), b_ions.end(), [](double d) { return d > 1e-6; });
                    size_t dot_product = std::accumulate(intensity_xls.begin(), intensity_xls.end(), 0.0);
                    dot_product = std::accumulate(intensity_linear.begin(), intensity_linear.end(), dot_product);
                    const double yFact = logfactorial_(y_ion_count);
                    const double bFact = logfactorial_(b_ion_count);
                    partial_loss_sub_score = log1p(dot_product) + yFact + bFact;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
                    ah.total_loss_score = total_loss_score;
                    ah.MIC = tlss_MIC;
                    ah.immonium_score = im_MIC;
                    ah.precursor_score = pc_MIC;
                    ah.err = tlss_err;
                    ah.Morph = tlss_Morph;
                    ah.modds = tlss_modds;
                    ah.pl_MIC = plss_MIC;
                    ah.pl_err = plss_err;
                    ah.pl_Morph = plss_Morph;
                    ah.pl_modds = plss_modds;
                    ah.pl_pc_MIC = plss_pc_MIC;
                    ah.pl_im_MIC = plss_im_MIC;
                    ah.cross_linked_nucleotide = cross_linked_nucleotide;
                    ah.total_MIC = total_MIC; 
                    // scores from shifted peaks
                    ah.marker_ions_score = marker_ions_sub_score;
                    ah.partial_loss_score = partial_loss_sub_score;

                    ah.NA_mod_index = NA_mod_index;
                    ah.isotope_error = isotope_error;

                    auto range = make_pair(intensity_linear.begin(), intensity_linear.end());
                    ah.ladder_score = ladderScore_(range) / (double)intensity_linear.size(); 
                    range = longestCompleteLadder_(intensity_linear.begin(), intensity_linear.end());
                    if (range.second != range.first)
                    {
                      ah.sequence_score = ladderScore_(range) / (double)intensity_linear.size();
                    }

                    RankScores rankscores = rankScores_(exp_spectrum, peak_matched);
                    ah.explained_peak_fraction = rankscores.explained_peak_fraction;
                    ah.wTop50 = rankscores.wTop50;


                    // does it have at least one shift from non-cross-linked AA to the neighboring cross-linked one
//                    const XLTags longest_tags = getLongestLadderWithShift(intensity_linear, intensity_xls);
                    const XLTags longest_tags = getLongestABYLadderWithShift(b_ions, y_ions, b_xl_ions, y_xl_ions);

#ifdef FILTER_BAD_SCORES_ID_TAGS
                    if (longest_tags.tag_XLed == 0) { continue; }
#endif
                    ah.tag_XLed = longest_tags.tag_XLed;
                    ah.tag_unshifted = longest_tags.tag_unshifted;
                    ah.tag_shifted = longest_tags.tag_shifted;

                    // combined score
                    const double tags = exp_spectrum.getFloatDataArrays()[1][0];
                    ah.score = OpenNuXL::calculateCombinedScore(ah, true, tags);
                    ah.n_theoretical_peaks = n_theoretical_peaks;

#ifdef DEBUG_OpenNuXL
                    OPENMS_LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
                    omp_set_lock(&(annotated_XLs_lock[scan_index]));
#endif
                    {
#ifdef ANNOTATED_QUANTILES
                      annotated_peptides_quantiles_XLs[scan_index](ah.total_loss_score + ah.partial_loss_score);

                      QQ_TIC[scan_index].insert(ah.total_MIC);
                      QQ_EXPLAINED_FRACTION[scan_index].insert(ah.explained_peak_fraction);
#endif
                      annotated_XLs[scan_index].emplace_back(move(ah));

                      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                      if (annotated_XLs[scan_index].size() >= 2 * report_top_hits)
                      {
                        std::partial_sort(annotated_XLs[scan_index].begin(), annotated_XLs[scan_index].begin() + report_top_hits, annotated_XLs[scan_index].end(), AnnotatedHit::hasBetterScore);
                        annotated_XLs[scan_index].resize(report_top_hits); 
                      }
                    }
#ifdef _OPENMP
                    omp_unset_lock(&(annotated_XLs_lock[scan_index]));
#endif
                  }
                } // for every nucleotide in the precursor
              }
            }
            else // fast scoring
            {
              for (auto & l = low_it; l != up_it; ++l)
              {
                const Size & scan_index = l->second.first;
                auto mod_combinations_it = mm.mod_combinations.begin();
                std::advance(mod_combinations_it, NA_mod_index);
                const String& precursor_rna_adduct = *mod_combinations_it->second.begin(); // For fast scoring it should be sufficient to only consider any of the adducts for this mass and formula (e.g., C-H3N vs U-H2O)
                MSSpectrum& exp_spectrum = spectra[scan_index];

                if (precursor_rna_adduct != "none" && skip_peptide_spectrum.find(exp_spectrum.getNativeID()) != skip_peptide_spectrum.end()) continue;


#ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                 // require at least one mass tag
                 if (exp_spectrum.getIntegerDataArrays()[IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
#endif
#ifdef _OPENMP
                omp_set_lock(&(annotated_peptides_lock[scan_index]));
                omp_set_lock(&(annotated_XLs_lock[scan_index]));
                ++nr_candidates[scan_index];
                omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                omp_unset_lock(&(annotated_peptides_lock[scan_index]));
#endif
                const int & isotope_error = l->second.second;
                const double & exp_pc_mass = l->first;

                // generate PSMs for spectrum[scan_index] and add them to annotated hits
                addPSMsTotalLossScoring_(
                  spectra[scan_index],
                  *cit, // string view on unmodified sequence
                  mod_pep_idx, // index of peptide mod
                  NA_mod_index, // index of NA mod
                  current_peptide_mass,
                  current_peptide_mass_without_NA,
                  exp_pc_mass,
                  iip,
                  isotope_error, 
                  total_loss_template_z1_b_ions, 
                  total_loss_template_z1_y_ions,
                  gaussian_mass_error,
                  fragment_mass_tolerance,
                  fragment_mass_tolerance_unit_ppm,
                  annotated_peptides[scan_index],
#ifdef _OPENMP
                  annotated_peptides_lock[scan_index],
#endif
                  report_top_hits
                );
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    OPENMS_LOG_INFO << "Proteins: " << count_proteins << endl;
    OPENMS_LOG_INFO << "Peptides: " << count_peptides << endl;
    OPENMS_LOG_INFO << "Processed peptides: " << processed_petides.size() << endl;

    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;
    progresslogger.startProgress(0, 1, "Post-processing PSMs... (spectra filtering)");

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Localization
    //

    // reload spectra from disc with same settings as before (important to keep same spectrum indices)
    spectra.clear(true);
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);    

    preprocessSpectra_(spectra, 
                       fragment_mass_tolerance, 
                       fragment_mass_tolerance_unit_ppm, 
                       false, // no single charge (false)
                       true, window_size, peak_count); // annotate charge (true)

    calculateNucleotideTags_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, nucleotide_to_fragment_adducts);
    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Post-processing PSMs... (localization of cross-links)");
    assert(spectra.size() == annotated_XLs.size());
    assert(spectra.size() == annotated_peptides.size());

    // remove all but top n scoring for localization (usually all but the first one)
    filterTopNAnnotations_(annotated_XLs, report_top_hits);
    filterTopNAnnotations_(annotated_peptides, report_top_hits);

    postScoreHits_(spectra, 
                   annotated_XLs, 
                   annotated_peptides, 
                   mm, 
                   fixed_modifications, 
                   variable_modifications, 
                   max_variable_mods_per_peptide, 
                   fragment_mass_tolerance, 
                   fragment_mass_tolerance_unit_ppm, 
                   all_feasible_fragment_adducts);

    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Post-processing PSMs... (annotation)");
    // remove all but top n scoring PSMs again
    // Note: this is currently necessary as postScoreHits_ might reintroduce nucleotide specific hits for fast scoring
    filterTopNAnnotations_(annotated_XLs, report_top_hits);
    filterTopNAnnotations_(annotated_peptides, report_top_hits);

    postProcessHits_(spectra, 
                     annotated_XLs, 
                     annotated_peptides, 
                     protein_ids, 
                     peptide_ids, 
                     mm, 
                     fixed_modifications, 
                     variable_modifications, 
                     max_variable_mods_per_peptide,
                     purities,
                     nr_candidates);

    progresslogger.endProgress();

    protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_mzml)});

    // reindex ids
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", getStringOption_("peptide:enzyme"));
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    param_pi.setValue("write_protein_sequence", "true");
    param_pi.setValue("write_protein_description", "true");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
      {
        return INPUT_FILE_EMPTY;       
      }
      else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
      {
        return UNEXPECTED_RESULT;
      }
      else
      {
        return UNKNOWN_ERROR;
      }
    } 

    StringList meta_values_to_export;
    meta_values_to_export.push_back("NuXL:total_loss_score");  
    meta_values_to_export.push_back("NuXL:partial_loss_score");  
    meta_values_to_export.push_back("CountSequenceIsTop");  
    meta_values_to_export.push_back("CountSequenceCharges");  
    meta_values_to_export.push_back("CountSequenceIsXL");  
    meta_values_to_export.push_back("CountSequenceIsPeptide");  
    meta_values_to_export.push_back("NuXL:MIC");  
    meta_values_to_export.push_back("NuXL:pl_pc_MIC");  
    meta_values_to_export.push_back("NuXL:pl_MIC");  
    meta_values_to_export.push_back("nr_candidates");
    meta_values_to_export.push_back("isotope_error");  

    // annotate RNPxl related information to hits and create report
    vector<RNPxlReportRow> csv_rows = RNPxlReport::annotate(spectra, peptide_ids, meta_values_to_export, marker_ions_tolerance);

#ifdef ANNOTATED_QUANTILES
      for (Size scan_index = 0; scan_index != peptide_ids.size(); ++scan_index)
      {
        PeptideIdentification& pi = peptide_ids[scan_index];
//        cout << "score\tpeptides\tXLs\ttype" << endl;
        for (auto & ph : pi.getHits())
        {
          if (ph.getMetaValue("NuXL:isXL") == "1")
          {
            ph.setMetaValue("NuXL:total_HS", (double)ph.getMetaValue("NuXL:total_HS") / (1.0 + p_square_quantile(annotated_peptides_quantiles_XLs[scan_index])));
          }
          else
          {
            ph.setMetaValue("NuXL:total_HS", (double)ph.getMetaValue("NuXL:total_HS") / (1.0 + p_square_quantile(annotated_peptides_quantiles_peptides[scan_index])));
          }
        }

        for (auto & ph : pi.getHits())
        {
          ph.setMetaValue("NuXL:QQ_TIC", QQ_TIC[scan_index].quantileOfValue(ph.getMetaValue("NuXL:total_MIC")));
          ph.setMetaValue("NuXL:QQ_EXPLAINED_FRACTION", QQ_EXPLAINED_FRACTION[scan_index].quantileOfValue(ph.getMetaValue("NuXL:explained_peak_fraction")));
        }
      }
#endif


/*
    // keep 10 bins with best 100 scores using priority queues
    vector<LargestElements> best_score_per_bin(10, LargestElements(100));
    for (size_t index = 0; index != peptide_ids.size(); ++index)
    {
      size_t bin_index = 10.0 * index / (double)peptide_ids.size();
      if (peptide_ids[index].getHits().empty()) continue;
      if (peptide_ids[index].getHits()[0].getMetaValue("target_decoy") == "target")
      {
        best_score_per_bin[bin_index].tryAdd(peptide_ids[index].getHits()[0].getScore());
      }
    }
*/
    if (generate_decoys) 
    {
      map<double, double, std::greater<double>> map_score2ppm;
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
         if (peptide_ids[index].getHits().empty()) continue;
         if (peptide_ids[index].getHits()[0].getMetaValue("target_decoy") == "target")
         {
           double ppm_error = peptide_ids[index].getHits()[0].getMetaValue(OpenMS::Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
           map_score2ppm[peptide_ids[index].getHits()[0].getScore()] = ppm_error; 
         }
      }

      // calculate mean ppm error from top scoring PSMs (max. 1000 considered)
      double mean(0), mean_negative(0), mean_positive(0);
      size_t c(0), c_negative(0), c_positive(0);
      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         mean += it->second;
         ++c;
         if (c > 1000) break;
      }
      if (c != 0) { mean /= c; }

      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         if (it->second > 0) continue; // skip positive ppm
         mean_negative += it->second;
         ++c_negative;
         if (c_negative > 1000) break;
      }
      if (c_negative != 0) { mean_negative /= c_negative; }

      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         if (it->second < 0) continue; // skip negative ppm
         mean_positive += it->second;
         ++c_positive;
         if (c_positive > 1000) break;
      }
      if (c_positive != 0) { mean_positive /= c_positive; }

      double sd(0), sd_negative(0), sd_positive(0);
      auto it = map_score2ppm.begin();
      for (size_t i = 0; i != c; ++i)
      {
         sd += pow(it->second - mean, 2.0);
         if (it->second < 0) sd_negative += pow(it->second - mean, 2.0);
         if (it->second > 0) sd_positive += pow(it->second - mean, 2.0);
         ++it;
      }

      if (c != 0) 
      { 
        sd = sqrt(1.0/static_cast<double>(c) * sd); 
        if (c_negative != 0) 
        { 
          sd_negative = sqrt(1.0/static_cast<double>(c_negative) * sd_negative);
        }
        if (c_positive != 0) 
        { 
          sd_positive = sqrt(1.0/static_cast<double>(c_positive) * sd_positive);
        }
        cout << "mean ppm error: " << mean << " sd: " << sd << " 5*sd: " << 5*sd << " calculated based on " << c << " best ids." << endl;
        cout << "mean negative ppm error: " << mean_negative << " sd: " << sd_negative << " 5*sd: " << 5*sd_negative << " calculated based on " << c_negative << " best ids." << endl;
        cout << "mean positive ppm error: " << mean_positive << " sd: " << sd_positive << " 5*sd: " << 5*sd_positive << " calculated based on " << c_positive << " best ids." << endl;
      } 
  
      if (filter_pc_mass_error && c != 0)
      {
        // as we are dealing with a very large search space, filter out all PSMs with mass error > 5 *sd
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
          vector<PeptideHit>& phs = peptide_ids[index].getHits();
          if (phs.empty()) continue;
          auto new_end = std::remove_if(phs.begin(), phs.end(),
              [&sd, &mean](const PeptideHit & ph) 
            { 
             return fabs((double)ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM)) - fabs(mean) > 5.0*sd; 
            });
          phs.erase(new_end, phs.end());
        }
        IDFilter::removeEmptyIdentifications(peptide_ids);
      }
      map_score2ppm.clear(); 

      if (impute_decoy_medians)
      {
        cout << "Imputing decoy medians." << endl;
        // calculate median score of decoys for specific meta value
        auto metaMedian = [](const vector<PeptideIdentification> & peptide_ids, const String name)->double
        {
          vector<double> decoy_XL_scores;
          for (const auto & pi : peptide_ids)
          {
            for (const auto & ph : pi.getHits())
            {
              const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
              if (!is_XL || ph.getMetaValue("target_decoy") != "decoy") continue;
              double score = ph.getMetaValue(name);
              decoy_XL_scores.push_back(score); 
            }
          }
          std::sort(decoy_XL_scores.begin(), decoy_XL_scores.end(), greater<double>());
          return Math::median(decoy_XL_scores.begin(), decoy_XL_scores.end());
        };
     
        // all medians 
        auto metaMean = [](const vector<PeptideIdentification> & peptide_ids, const String name)->double
        {
          vector<double> decoy_XL_scores;
          for (const auto & pi : peptide_ids)
          {
            for (const auto & ph : pi.getHits())
            {
              const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
              if (!is_XL) continue;
              double score = ph.getMetaValue(name);
              decoy_XL_scores.push_back(score); 
            }
          }
          std::sort(decoy_XL_scores.begin(), decoy_XL_scores.end(), greater<double>());
          return Math::mean(decoy_XL_scores.begin(), decoy_XL_scores.end());
        };

        map<String, double> medians;
        for (const String mn : { "NuXL:marker_ions_score", "NuXL:partial_loss_score", "NuXL:pl_MIC", "NuXL:pl_err", "NuXL:pl_Morph", "NuXL:pl_modds", "NuXL:pl_pc_MIC", "NuXL:pl_im_MIC" })
        {
//          medians[mn] = metaMedian(peptide_ids, mn);
           medians[mn] = metaMean(peptide_ids, mn);
        }

        size_t imputed(0);
        for (auto & pi : peptide_ids)
        {
          for (auto & ph : pi.getHits())
          {
             const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
             if (!is_XL) 
             { 
               for (const String mn : { "NuXL:marker_ions_score", "NuXL:partial_loss_score", "NuXL:pl_MIC", "NuXL:pl_err", "NuXL:pl_Morph", "NuXL:pl_modds", "NuXL:pl_pc_MIC", "NuXL:pl_im_MIC" })
               {
                 ph.setMetaValue(mn, medians[mn]);   // impute missing with medians
               }
               ++imputed;
             }
          }
          pi.assignRanks();
        }
        cout << "Imputed XL features in " << imputed << " linear peptides." << endl;
      }


#ifdef SVM_RECALIBRATE
      ///////////////////////////////////////// SVM score recalibration
      // find size of minority class
      map<double, size_t> pep_t;
      map<double, size_t> pep_d;
      map<double, size_t> XL_t;
      map<double, size_t> XL_d;

      // map scores to index (sorted)
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
         if (peptide_ids[index].getHits().empty()) continue;
         const PeptideHit& ph = peptide_ids[index].getHits()[0]; // get best match to current spectrum
         bool is_target = ph.getMetaValue("target_decoy") == "target";
         bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
         double score = ph.getScore();
         if (is_target)
         {
           if (is_XL) 
           {
             XL_t[score] = index; 
           }
           else 
           {
             pep_t[score] = index; 
           }
         }
         else
         {
           if (is_XL) 
           {
             XL_d[score] = index; 
           }
           else 
           {
             pep_d[score] = index; 
           }
         }
      }
      size_t minority_class = std::min({pep_t.size(), pep_d.size(), XL_t.size(), XL_d.size()});
      cout << "Peptide (target/decoy)\t XL (target/decoy):" << endl;
      cout << pep_t.size() << "\t" << pep_d.size() << "\t" << XL_t.size() << "\t" << XL_d.size() << endl; 

      if (minority_class > 250) minority_class = 250;

      // keep only top elements (=highest scoring hits) for all classes
      pep_t.erase(pep_t.begin(), next(pep_t.begin(), pep_t.size() - minority_class));
      pep_d.erase(pep_d.begin(), next(pep_d.begin(), pep_d.size() - minority_class));
      XL_t.erase(XL_t.begin(), next(XL_t.begin(), XL_t.size() - minority_class));
      XL_d.erase(XL_d.begin(), next(XL_d.begin(), XL_d.size() - minority_class));

      // training indices = index of peptide id that correspond to top scoring hits (of the 4 classes)
      set<size_t> training_indices;
      for (auto & l : {pep_t, pep_d, XL_t, XL_d})
      {
        for (auto & i : l) training_indices.insert(i.second);
      }

      if (minority_class > 100)
      {
        SimpleSVM::PredictorMap predictors;
        map<Size, Int> labels;

        // copy all scores in predictors ("score" + all from feature_set). 
        // Only add labels for balanced training set
        size_t current_row(0);
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
           const vector<PeptideHit>& phits = peptide_ids[index].getHits();
           for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank)
           {
             const PeptideHit& ph = phits[psm_rank];
             bool is_target = ph.getMetaValue("target_decoy") == "target";
             double score = ph.getScore();
             // predictors["score"].push_back(score);
             predictors["length"].push_back(ph.getSequence().size());
             for (auto & f : feature_set_)
             {
               double value = ph.getMetaValue(f);
               predictors[f].push_back(value);
             }
             // only add label for training data (rank = 0 and previously selected for training)
             if (psm_rank == 0 && training_indices.count(index) > 0)
             {
               labels[current_row] = is_target;
             }
             ++current_row;
           }
        }
  
        SimpleSVM svm;
        Param svm_param = svm.getParameters();        
        svm_param.setValue("kernel", "linear");
//        svm_param.setValue("kernel", "RBF");
//        svm_param.setValue("log2_C", ListUtils::create<double>("0"));
//        svm_param.setValue("log2_gamma", ListUtils::create<double>("1"));
        svm.setParameters(svm_param);
        svm.setup(predictors, labels);
        vector<SimpleSVM::Prediction> predictions;
        cout << "Predicting class probabilities:" << endl;
        svm.predict(predictions);
        std::map<String, double> feature_weights;
        svm.getFeatureWeights(feature_weights);
 /*
        cout << "Feature weights:" << endl;
        for (const auto& m : feature_weights)
        {
          cout << m.first << "\t" << m.second << endl;
        }
        cout << "Feature scaling:" << endl;
        auto feature_scaling = svm.getScaling();
        for (const auto& m : feature_scaling)
        {
          cout << m.first << "\t" << m.second.first << "\t" << m.second.second << endl;
        }
*/
        size_t psm_index(0);
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
           vector<PeptideHit>& phits = peptide_ids[index].getHits();
           for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank, ++psm_index)
           {
             PeptideHit& ph = phits[psm_rank];
             ph.setScore(predictions[psm_index].probabilities[1]); // set probability of being a target as score
           }
          peptide_ids[index].assignRanks();    
        }
        // IdXMLFile().store(out_idxml + "_svm.idXML", protein_ids, peptide_ids);
      }
#endif
      fdr.apply(peptide_ids); 
  
      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

      // generate filtered results
#ifdef FILTER_RANKS
      for (auto & pi : peptide_ids)
      {
        auto & phs = pi.getHits();
        if (phs.empty()) continue;
        if (static_cast<int>(phs[0].getMetaValue("NuXL:isXL")) == 0) continue; // only rerank cross-links

        double max_total_Morph = std::max_element(phs.begin(), phs.end(), 
           [] (PeptideHit const& lhs, PeptideHit const& rhs) 
           {
             return (double)lhs.getMetaValue("NuXL:total_Morph") < (double)rhs.getMetaValue("NuXL:total_Morph"); 
           })->getMetaValue("NuXL:total_Morph");

        auto new_end = std::remove_if(phs.begin(), phs.end(),
            [&max_total_Morph](const PeptideHit & ph) 
          { 
           return fabs((double)ph.getMetaValue("NuXL:total_Morph") - max_total_Morph) > 1e-4; 
          });
        phs.erase(new_end, phs.end());
      } 
#else
      IDFilter::keepNBestHits(peptide_ids, 1);
#endif       
      IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);

      // split PSMs into XLs and non-XLs but keep only best one of both
      vector<PeptideIdentification> pep_pi;
      vector<PeptideIdentification> xl_pi;
      for (const auto & pi : peptide_ids)
      {
        vector<PeptideHit> pep_ph;
        vector<PeptideHit> xl_ph;
        for (const auto & ph : pi.getHits())
        {
           if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
           { // only add best hit
             if (pep_ph.empty() && xl_ph.empty()) pep_ph.push_back(ph); 
           }
           else
           {
             if (pep_ph.empty() && xl_ph.empty()) xl_ph.push_back(ph); 
           }
        }
        if (!pep_ph.empty()) { pep_pi.push_back(pi); pep_pi.back().setHits(pep_ph); }
        if (!xl_ph.empty()) { xl_pi.push_back(pi); xl_pi.back().setHits(xl_ph); }
      }

      // calculate FDRs separately
      fdr.apply(xl_pi); 
      fdr.apply(pep_pi);

      if (debug_level_ > 1)
      {

      }

      IDFilter::removeDecoyHits(xl_pi);
      IDFilter::removeDecoyHits(pep_pi);

      if (peptide_FDR > 0.0 && peptide_FDR < 1.0)
      {
         IDFilter::filterHitsByScore(pep_pi, peptide_FDR); 
      }

      {
        String out2(out_idxml);
        out2.substitute(".idXML", "_");
        vector<ProteinIdentification> tmp_prots = protein_ids;
        IDFilter::removeUnreferencedProteins(tmp_prots, pep_pi);
        IdXMLFile().store(out2 + String::number(peptide_FDR, 4) + "_peptides.idXML", tmp_prots, pep_pi);
      }

      std::replace(XL_FDR.begin(), XL_FDR.end(), 0.0, 1.0);
      sort(XL_FDR.begin(), XL_FDR.end(), greater<double>()); // important: sort by threshold (descending) to generate results by applying increasingly stringent q-value filters
      for (double xlFDR : XL_FDR)
      { 
        OPENMS_LOG_INFO << "Writing XL results at xl-FDR: " << xlFDR << endl;
        if (xlFDR > 0.0 && xlFDR < 1.0)
        {
          IDFilter::filterHitsByScore(xl_pi, xlFDR);
        }
        {
          String out2(out_idxml);
          out2.substitute(".idXML", "_");
          vector<ProteinIdentification> tmp_prots = protein_ids;
          IDFilter::removeUnreferencedProteins(tmp_prots, xl_pi);
          // PSM result
          IdXMLFile().store(out2 + String::number(xlFDR, 4) + "_XLs.idXML", tmp_prots, xl_pi);

          // protein result
          if (!out_tsv.empty())
          {
            out2 = out_tsv;
            out2.substitute(".tsv", "_proteins");
            TextFile tsv_file;
            annotateProteinModificationForTopHits_(tmp_prots, xl_pi, tsv_file);
            tsv_file.store(out2 + String::number(xlFDR, 4) + "_XLs.tsv");
          }
        }
      }
     
      String percolator_executable = getStringOption_("percolator_executable");
      bool sufficient_PSMs_for_score_recalibration = (xl_pi.size() + pep_pi.size()) > 1000;
      if (!percolator_executable.empty() && sufficient_PSMs_for_score_recalibration) // only try to call percolator if we have some PSMs
      {
        // run percolator on idXML
        String perc_out = out_idxml;
        perc_out.substitute(".idXML", "_perc.idXML");
        String weights_out = out_idxml;
        weights_out.substitute(".idXML", ".weights");
        String pin = out_idxml;
        pin.substitute(".idXML", ".tsv");

        QStringList process_params;
        process_params << "-in" << out_idxml.toQString()
                       << "-out" << perc_out.toQString()
                       << "-percolator_executable" << percolator_executable.toQString()
                       << "-train_best_positive" 
                       << "-score_type" << "svm"
                       << "-unitnorm"
                       << "-post_processing_tdc"
                       << "-nested_xval_bins" << "3"
                       << "-weights" << weights_out.toQString()
                       << "-out_pin" << pin.toQString();

        if (getStringOption_("peptide:enzyme") == "Lys-C")
        {
          process_params << "-enzyme" << "lys-c";
        }
//        process_params << "-out_pout_target" << "merged_target.tab" << "-out_pout_decoy" << "merged_decoy.tab";

        TOPPBase::ExitCodes exit_code = runExternalProcess_(QString("PercolatorAdapter"), process_params);

        if (exit_code != EXECUTION_OK) 
        { 
          OPENMS_LOG_WARN << "Score recalibration failed." << endl; 
        }
        else
        { 
          // load back idXML
          IdXMLFile().load(perc_out, protein_ids, peptide_ids);
 
          // generate filtered results
          IDFilter::keepNBestHits(peptide_ids, 1);
          IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);

	  // annotate RNPxl related information to hits and create report
          vector<RNPxlReportRow> csv_rows_percolator = RNPxlReport::annotate(spectra, peptide_ids, meta_values_to_export, marker_ions_tolerance);

          // save report
          if (!out_tsv.empty())
          {
            TextFile csv_file;
            csv_file.addLine(RNPxlReportRowHeader().getString("\t", meta_values_to_export));
            for (const RNPxlReportRow r : csv_rows_percolator)
            {
              csv_file.addLine(r.getString("\t"));
            }
            const String out_percolator_tsv = FileHandler::stripExtension(out_tsv) + "_perc.tsv";
            csv_file.store(out_percolator_tsv);
          }


          // split PSMs into XLs and non-XLs but keep only best one of both
          vector<PeptideIdentification> pep_pi;
          vector<PeptideIdentification> xl_pi;
          for (const auto & pi : peptide_ids)
          {
            vector<PeptideHit> pep_ph, xl_ph;
            for (const auto & ph : pi.getHits())
            {
               if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
               { // only add best hit
                 if (pep_ph.empty() && xl_ph.empty()) pep_ph.push_back(ph); 
               }
               else
               {
                 if (pep_ph.empty() && xl_ph.empty()) xl_ph.push_back(ph); 
               }
            }
            if (!pep_ph.empty()) { pep_pi.push_back(pi); pep_pi.back().setHits(pep_ph); }
            if (!xl_ph.empty()) { xl_pi.push_back(pi); xl_pi.back().setHits(xl_ph); }
          }

          // calculate FDRs separately
          fdr.apply(xl_pi); 
          fdr.apply(pep_pi);

          IDFilter::removeDecoyHits(xl_pi);
          IDFilter::removeDecoyHits(pep_pi);
 
          if (peptide_FDR > 0.0 && peptide_FDR < 1.0)
          {
            IDFilter::filterHitsByScore(pep_pi, peptide_FDR); 
          }
          {
            String out2(out_idxml);
            out2.substitute(".idXML", "_perc_");
            vector<ProteinIdentification> tmp_prots = protein_ids;
            IDFilter::removeUnreferencedProteins(tmp_prots, pep_pi);
            IdXMLFile().store(out2 + String::number(peptide_FDR, 4) + "_peptides.idXML", tmp_prots, pep_pi);
          }
 
          std::replace(XL_FDR.begin(), XL_FDR.end(), 0.0, 1.0); // treat disabled filtering as 100% FDR
          sort(XL_FDR.begin(), XL_FDR.end(), greater<double>()); // important: sort by threshold (descending) to generate results by applying increasingly stringent q-value filters
          for (double xlFDR : XL_FDR)
          {
            OPENMS_LOG_INFO << "Writing XL results at xl-FDR: " << xlFDR << endl;
            if (xlFDR > 0.0 && xlFDR < 1.0)
            {
              IDFilter::filterHitsByScore(xl_pi, xlFDR);
            }
            {
              String out2(out_idxml);
              out2.substitute(".idXML", "_perc_");
              vector<ProteinIdentification> tmp_prots = protein_ids;
              IDFilter::removeUnreferencedProteins(tmp_prots, xl_pi);

              // store PSM results
              IdXMLFile().store(out2 + String::number(xlFDR, 4) + "_XLs.idXML", tmp_prots, xl_pi);
              // store protein results
              if (!out_tsv.empty())
              {
                out2 = out_tsv;
                out2.substitute(".tsv", "_perc_proteins");
                TextFile tsv_file;
                annotateProteinModificationForTopHits_(tmp_prots, xl_pi, tsv_file);
                tsv_file.store(out2 + String::number(xlFDR, 4) + "_XLs.tsv");
              }
            }
          }
        }
      }
      else
      {
        if (sufficient_PSMs_for_score_recalibration == false) OPENMS_LOG_WARN << "Too few PSMs for score recalibration. Skipped." << endl;
      }
    }
    else  // no decoys
    {
      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);
    }

    // save report
    if (!out_tsv.empty())
    {
      TextFile csv_file;
      csv_file.addLine(RNPxlReportRowHeader().getString("\t", meta_values_to_export));
      for (Size i = 0; i != csv_rows.size(); ++i)
      {
        csv_file.addLine(csv_rows[i].getString("\t"));
      }
      csv_file.store(out_tsv);
    }

 
 #ifdef _OPENMP
    // free locks
    for (size_t i = 0; i != annotated_XLs_lock.size(); i++) { omp_destroy_lock(&(annotated_XLs_lock[i])); }
    for (size_t i = 0; i != annotated_peptides_lock.size(); i++) { omp_destroy_lock(&(annotated_peptides_lock[i])); }
 #endif

    return EXECUTION_OK;
  }

  static void postScorePartialLossFragments_(const Size peptide_size,
                                  const PeakSpectrum &exp_spectrum,
                                  double fragment_mass_tolerance,
                                  bool fragment_mass_tolerance_unit_ppm,
                                  const PeakSpectrum &partial_loss_spectrum_z1,
                                  const PeakSpectrum &partial_loss_spectrum_z2,
                                  const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                                  float &partial_loss_sub_score,
                                  float &marker_ions_sub_score,
                                  float &plss_MIC, 
                                  float &plss_err, 
                                  float &plss_Morph,
                                  float &plss_modds) 
  {
    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    if (!marker_ions_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                                             marker_ions_sub_score_spectrum_z1,
                                             marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[IA_CHARGE_INDEX]);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }

    if (!partial_loss_spectrum_z1.empty()) // check if we generated partial loss spectra
    {
      vector<double> intensity_sum(peptide_size, 0.0);
      MSSpectrum const * pl_spec = &partial_loss_spectrum_z1;
      if (exp_pc_charge >= 3)
      {
        pl_spec = &partial_loss_spectrum_z2;
      }
      partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, 
                                                    fragment_mass_tolerance_unit_ppm,
                                                    exp_spectrum, 
                                                    exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                                                    *pl_spec,
                                                    pl_spec->getIntegerDataArrays()[IA_CHARGE_INDEX],
                                                    intensity_sum);
                                                    
      auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                          fragment_mass_tolerance_unit_ppm,
                                                          exp_spectrum,
                                                          exp_spectrum.getIntegerDataArrays()[IA_CHARGE_INDEX],
                                                          *pl_spec,
                                                          pl_spec->getIntegerDataArrays()[IA_CHARGE_INDEX]);      
      plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
      plss_Morph = pl_sub_scores.score;

      // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
      plss_err = plss_Morph > 2 ? pl_sub_scores.err : 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;

      const float fragment_mass_tolerance_Da = 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;
      plss_modds = matchOddsScore_(pl_spec->size(), 
        fragment_mass_tolerance_Da,
        exp_spectrum.size(),
        exp_spectrum.back().getMZ(),
        (int)plss_Morph);
    }
#ifdef DEBUG_OpenNuXL
    OPENMS_LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
    // cap plss_err to something larger than the mean_mz * max_ppm_error
    float ft_da = fragment_mass_tolerance_unit_ppm ? fragment_mass_tolerance * 1e-6 * 1000.0 : fragment_mass_tolerance;
    if (plss_err > ft_da) plss_err = ft_da;
  }
};

#ifdef FILTER_AMBIGIOUS_PEAKS
map<double, double> OpenNuXL::mass2high_frequency_ = {};
#endif

int main(int argc, const char** argv)
{
  OpenNuXL tool;
  return tool.main(argc, argv);
}

