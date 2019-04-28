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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>

#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>
#include <OpenMS/ANALYSIS/RNPXL/MorpheusScore.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlFragmentAnnotationHelper.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

// preprocessing and filtering
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <boost/regex.hpp>

#include <map>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//#define DEBUG_RNPXLSEARCH 

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

class RNPxlSearch :
  public TOPPBase
{
  // fast or all-shifts scoring mode
  bool fast_scoring_ = true;

  // nucleotides can form cross-link
  set<char> can_xl_;


public:
  RNPxlSearch() :
    TOPPBase("RNPxlSearch", "Annotate RNA/DNA-peptide cross-links in MS/MS spectra.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("database", "<file>", "", "input file ");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_tsv", "<file>", "", "tsv output file", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance (+/- around precursor m/z)", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, false);

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)", false, false);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 6, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:max_size", "<num>", 1e6, "Maximum size a peptide may have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);

    // RNPxl specific
    registerTOPPSubsection_("RNPxl", "RNPxl Options");
    registerIntOption_("RNPxl:length", "", 2, "Oligonucleotide maximum length. 0 = disable search for RNA variants.", false);

    registerStringOption_("RNPxl:sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)", false);

    registerStringList_("RNPxl:target_nucleotides", 
                        "", 
                        {"A=C10H14N5O7P", "C=C9H14N3O8P", "G=C10H14N5O8P", "U=C9H13N2O9P"}, 
                        "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG", 
                        false, 
                        false);

    registerStringList_("RNPxl:nt_groups",
        "",
        {},
	"Restrict which nucleotides can cooccur in a precursor adduct to be able to search both RNA and DNA (Formate e.g.: AU CG).",
        false,
        false);

    registerStringList_("RNPxl:mapping", "", {"A->A", "C->C", "G->G", "U->U"}, "format: source->target e.g. A->A, ..., U->U, U->X", false, false);

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
                        "format: [target nucleotide]:[formula] or [precursor adduct]->[fragment adduct formula];[name]: e.g., 'U:C9H10N2O5;U-H3PO4' or 'U:U-H2O->C9H11N2O8P1;U-H2O',", 
                        false, 
                        false);

    registerStringList_("RNPxl:modifications", "", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3", false, false);

    registerStringOption_("RNPxl:scoring", "<method>", "fast", "Scoring algorithm used in prescoring (fast: total-loss, slow: all losses).", false, false);
    setValidStrings_("RNPxl:scoring", {"fast", "slow"});

    registerFlag_("RNPxl:decoys", "Generate decoy sequences and spectra.");

    registerFlag_("RNPxl:CysteineAdduct", "Use this flag if the +152 adduct is expected.", true);
    registerFlag_("RNPxl:filter_fractional_mass", "Use this flag to filter non-crosslinks by fractional mass.", true);
    registerFlag_("RNPxl:carbon_labeled_fragments", "Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).", true);
    registerFlag_("RNPxl:only_xl", "Only search cross-links and ignore non-cross-linked peptides.", true);

    registerDoubleOption_("RNPxl:filter_small_peptide_mass", "<threshold>", 600.0, "Filter precursor that can only correspond to non-crosslinks by mass.", false, true);
    registerDoubleOption_("RNPxl:marker_ions_tolerance", "<tolerance>", 0.05, "Tolerance used to determine marker ions (Da).", false, true);
  }

  struct FragmentAdductDefinition_
  {
    EmpiricalFormula formula; // formula
    String name;  // name used in annotation
    double mass = 0;

    FragmentAdductDefinition_() = default;

    FragmentAdductDefinition_(const FragmentAdductDefinition_&) = default;

    FragmentAdductDefinition_(FragmentAdductDefinition_&&) = default;

    FragmentAdductDefinition_& operator=(const FragmentAdductDefinition_&) = default;

    FragmentAdductDefinition_& operator=(FragmentAdductDefinition_&&) = default;

    bool operator<(const FragmentAdductDefinition_& other) const
    {
      String fa = formula.toString();
      String fb = other.formula.toString();
      return std::tie(mass, fa, name) < std::tie(other.mass, fb, other.name);
    }

    bool operator==(const FragmentAdductDefinition_& other) const
    {
      return std::tie(formula, name) == std::tie(other.formula, other.name);
    }

  };

  // fast (flat) data structure to store feasible x-,y-,a-ion fragment adducts and observable marker ions
  using NucleotideToFeasibleFragmentAdducts = pair<char, vector<FragmentAdductDefinition_> >;

  // stores the fragment adducts and marker ions for a given precursor adduct
  struct MS2AdductsOfSinglePrecursorAdduct
  {
    vector<NucleotideToFeasibleFragmentAdducts> feasible_adducts;
    vector<FragmentAdductDefinition_> marker_ions;
  };

  // helper struct to facilitate parsing of parameters (modifications, nucleotide adducts, ...)
  struct RNPxlParameterParsing
  {
    // Map a nucleotide (e.g. U to all possible fragment adducts)
    using NucleotideToFragmentAdductMap = map<char, set<FragmentAdductDefinition_> >;
    // @brief Parse tool parameter to create map from target nucleotide to all its fragment adducts
    // It maps a single letter nucleotide (e.g., 'T', 'C', ...)
    // to the maximum set of fragment adducts that may arise if the nucleotide is cross-linked.
    // Losses, that might reduce this set, are not considered in this data structure and handled later
    // when specific precursor adducts are considered.
    static NucleotideToFragmentAdductMap getTargetNucleotideToFragmentAdducts(StringList fragment_adducts);


    // @brief Determines the fragment adducts and marker ions for a given precursor.
    // The precursor adduct (the oligo including losses, e.g.: "TC-H3PO4") is mapped to all contained nucleotides
    // and their marker ions. In addition, each cross-linkable nucleotide is mapped to its chemically feasible fragment adducts.
    // Chemical feasible means in this context, that the fragment or marker ion adduct is a subformula of the precursor adduct.
    static MS2AdductsOfSinglePrecursorAdduct getFeasibleFragmentAdducts(
      const String& exp_pc_adduct,
      const String& exp_pc_formula,
      const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts,
      const set<char>& can_xl
    );

    // Maps a precursor adduct (e.g.: "UU-H2O") to all chemically feasible fragment adducts.
    using PrecursorsToMS2Adducts = map<string, MS2AdductsOfSinglePrecursorAdduct>;

    // @brief Calculate all chemically feasible fragment adducts for all possible precursor adducts
    // Same as getFeasibleFragmentAdducts but calculated from all precursor adducts
    static PrecursorsToMS2Adducts getAllFeasibleFragmentAdducts(const RNPxlModificationMassesResult& precursor_adducts,
                                                          const NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts,
                                                          const set<char>& can_xl);

  };

  /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
  /// floats need to be initialized to zero as default
  struct AnnotatedHit
  {
    StringView sequence;
    SignedSize peptide_mod_index = 0; // enumeration index of the non-RNA peptide modification
    Size rna_mod_index = 0; // index of the RNA modification
    int isotope_error = 0; // wheter the hit has been matched with isotopic misassignment

    static constexpr const char UNKNOWN_NUCLEOTIDE = '?';
    char cross_linked_nucleotide = UNKNOWN_NUCLEOTIDE;
    // main score
    float score = 0;

    float total_loss_score = 0;
    // total loss morpheus related subscores
    float MIC = 0;
    float err = 0;
    float Morph = 0;

    // partial loss morpheus related subscores
    float pl_MIC = 0;
    float pl_err = 0;
    float pl_Morph = 0;

    // complete TIC fraction of explained peaks
    float total_MIC = 0;

    // subscores
    float partial_loss_score = 0;
    float immonium_score = 0;
    float precursor_score = 0;
    float a_ion_score = 0;
    float marker_ions_score = 0;

    float best_localization_score = 0;
    String localization_scores;
    String best_localization;  
    std::vector<PeptideHit::PeakAnnotation> fragment_annotations;

    static bool hasBetterScore(const AnnotatedHit& a, const AnnotatedHit& b)
    {
      return a.score > b.score;
    }
  };

  static float calculateCombinedScore(const AnnotatedHit& ah, const bool isXL)
  {
	return
	    + 0.995
		+ 7.142 * (  0.058 * ah.total_loss_score - 0.900)
		+ 0.802 * (  33.35 * ah.immonium_score - 1.148)
		+ 0.327 * (  73.64 * ah.precursor_score - 0.821)
		+ 0.748 * ( 22.014 * ah.marker_ions_score - 0.903)
		+ 0.746 * (  0.043 * ah.partial_loss_score - 0.472)
		- 1.788 * (  301.0 * ah.err - 1.771)
		- 1.292 * (240.825 * ah.pl_err - 1.323)
		+ 2.324 * static_cast<int>(isXL);
	/* old version
	return 
	  2.493
	  + 7.239 * (0.058 * ah.total_loss_score - 0.900)
	  + 1.381 * (26.965 * ah.marker_ions_score - 0.300)
	  + 1.178 * (0.043 * ah.partial_loss_score - 0.472)
	  - 1.934 * (300.828 * ah.err - 1.774)
	  - 0.358 * (240.441 * ah.pl_err - 1.316);
	*/
  }



  /* @brief Filter spectra to remove noise.
     Parameter are passed to spectra filter.
   */
  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra, bool annotate_charge = false)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // sort by mz
      exp[exp_index].sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(exp[exp_index], 
                                         fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
                                         1, 3, 
                                         false, 
                                         2, 10, 
                                         single_charge_spectra, 
                                         annotate_charge);
    #ifdef DEBUG_RNPXLSEARCH
      cout << "after deisotoping..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
//      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
      cout << endl;
    #endif

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);

    #ifdef DEBUG_RNPXLSEARCH
      cout << "after mower..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
    
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

    #ifdef DEBUG_RNPXLSEARCH
      cout << "after nlargest..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      cout << "Fragment charges in spectrum: " << exp_index  << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
 
      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
  
    #ifdef DEBUG_RNPXLSEARCH
      cout << "after sort..." << endl;
      cout << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != exp[exp_index].size(); ++i) cout << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << endl;
      if (exp[exp_index].getIntegerDataArrays().size())
        for (Size i = 0; i != exp[exp_index].size(); ++i) 
          cout  << exp[exp_index][i].getMZ() << "\t" << exp[exp_index][i].getIntensity() << "\t"  << exp[exp_index].getIntegerDataArrays()[0][i] << endl;
    #endif
    }

//    MzMLFile().store(String("RNPxlSearch_a_") + String((int)annotate_charge) + ".mzML", exp);
  }



  // prefix used to denote marker ions in fragment  annotations
  static constexpr const char* ANNOTATIONS_MARKER_ION_PREFIX = "MI:";

  // helper class that adds special ions not covered by TheoreticalSpectrumGenerator
  class RNPxlFragmentIonGenerator
  {
    public:
    // add RNA-marker ions of charge 1
    // this includes the protonated nitrogenous base and all shifts (e.g., U-H2O, U'-H20, ...)
    static void addMS2MarkerIons(
      const vector<FragmentAdductDefinition_>& marker_ions,
      PeakSpectrum& spectrum,
      PeakSpectrum::IntegerDataArray& spectrum_charge,
      PeakSpectrum::StringDataArray& spectrum_annotation);

    static void addShiftedImmoniumIons(
      const String & unmodified_sequence,
      const String & fragment_shift_name,
      const double fragment_shift_mass,
      PeakSpectrum & partial_loss_spectrum,
      PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge,
      PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation);

    static void generatePartialLossSpectrum(const String& unmodified_sequence,
                                     const double& fixed_and_variable_modified_peptide_weight,
                                     const String& precursor_rna_adduct,
                                     const double& precursor_rna_weight,
                                     const int& precursor_charge,
                                     const vector<FragmentAdductDefinition_>& partial_loss_modification,
                                     const PeakSpectrum& patial_loss_template_z1,
                                     const PeakSpectrum& patial_loss_template_z2,
                                     const PeakSpectrum& patial_loss_template_z3,
                                     PeakSpectrum& partial_loss_spectrum);
    static void addPrecursorWithCompleteRNA_(const double fixed_and_variable_modified_peptide_weight,
                                      const String & precursor_rna_adduct,
                                      const double precursor_rna_weight,
                                      const int charge,
                                      PeakSpectrum & partial_loss_spectrum,
                                      MSSpectrum::IntegerDataArray & partial_loss_spectrum_charge,
                                      MSSpectrum::StringDataArray & partial_loss_spectrum_annotation);

    static void addSpecialLysImmonumIons(const String& unmodified_sequence,
                                      PeakSpectrum &spectrum,
                                      PeakSpectrum::IntegerDataArray &spectrum_charge, 
                                      PeakSpectrum::StringDataArray &spectrum_annotation);
  };


  /* @brief Localization step of the cross-link identification engine.
   * Given a top scoring candidate (based on total loss spectrum) it:
   *  - generates all fragment adducts based on the attached precursor adduct
   *  - annotated peaks
   *  - calculates an additive score that considers the presence or absence of evidence for a cross-linking site
   *  - the maximum score is reported
   */
  void postScoreHits_(const PeakMap& exp, 
                      vector<vector<AnnotatedHit> >& annotated_hits, 
                      Size top_hits, 
                      const RNPxlModificationMassesResult& mm, 
                      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
                      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
                      Size max_variable_mods_per_peptide, 
                      const TheoreticalSpectrumGenerator& partial_loss_spectrum_generator, 
                      double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, 
                      const RNPxlParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    assert(exp.size() == annotated_hits.size());

    #ifdef DEBUG_RNPXLSEARCH
      LOG_DEBUG << exp.size() << " : " << annotated_hits.size() << endl;
    #endif

    SpectrumAlignment spectrum_aligner;
    Param pa = spectrum_aligner.getParameters();
    pa.setValue("tolerance", fragment_mass_tolerance, "Defines the absolute (in Da) or relative (in ppm) tolerance in the alignment");
    pa.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm ? "true" : "false");  
    spectrum_aligner.setParameters(pa);

    // remove all but top n scoring for localization (usually all but the first one)
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      const Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
      annotated_hits[scan_index].shrink_to_fit();
    }

    // If we did a (total-loss) only fast scoring, PSMs were not associated with a nucleotide.
    // To make the localization code work for both fast and slow (all-shifts) scoring,
    // we copy PSMs for every cross-linkable nucleotide present in the precursor.
    if (fast_scoring_)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
      {
        vector<AnnotatedHit> new_hits;

        // for each PSM
        for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
        {
          // determine RNA on precursor from index in map
          auto mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, annotated_hits[scan_index][i].rna_mod_index);
          const String precursor_rna_adduct = *mod_combinations_it->second.begin();
          const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;

          // just copy non-cross-linked peptide PSMs
          if (precursor_rna_adduct == "none") 
          {
            new_hits.push_back(annotated_hits[scan_index][i]);
            continue;
          }

          // if we have a cross-link, copy PSM information for each cross-linkable nucleotides
          for (auto const & c : feasible_MS2_adducts)
          {
            AnnotatedHit a(annotated_hits[scan_index][i]);
            a.cross_linked_nucleotide = c.first; // nucleotide
            new_hits.push_back(a);
          }
        }
        annotated_hits[scan_index].swap(new_hits);
      }
    }

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

        // determine RNA on precursor from index in map
        std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, a.rna_mod_index);
        const String precursor_rna_adduct = *mod_combinations_it->second.begin();
        const double precursor_rna_weight = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

        // we don't localize on non-cross-links
        if (precursor_rna_adduct == "none") { continue; }

        // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
        // 1. get all possible RNA fragment shifts in the MS2 (based on the precursor RNA/DNA)
        LOG_DEBUG << "precursor_rna_adduct: "  << precursor_rna_adduct << endl;
        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_rna_adduct).feasible_adducts;

        if (feasible_MS2_adducts.empty()) { continue; } // should not be the case - check case of no nucleotide but base fragment ?

        // 2. retrieve the (nucleotide specific) fragment adducts for the cross-linked nucleotide (annotated in main search)
        auto nt_to_adducts = std::find_if(feasible_MS2_adducts.begin(), feasible_MS2_adducts.end(),
          [&a](NucleotideToFeasibleFragmentAdducts const & item)
          {
            return (item.first == a.cross_linked_nucleotide);
          });

        OPENMS_POSTCONDITION(nt_to_adducts != feasible_MS2_adducts.end(), "Nucleotide not found in mapping to feasible adducts.")

        const vector<FragmentAdductDefinition_>& partial_loss_modification = nt_to_adducts->second;

        // get marker ions (these are not specific to the cross-linked nucleotide but also depend on the whole oligo bound to the precursor)
        const vector<FragmentAdductDefinition_>& marker_ions = all_feasible_adducts.at(precursor_rna_adduct).marker_ions;

        // generate total loss spectrum for the fixed and variable modified peptide (without RNA) (using the settings for partial loss generation)
        // but as we also add the abundant immonium ions for charge 1 and precursor ions for all charges to get a more complete annotation
        // (these have previously not been used in the scoring of the total loss spectrum)
        PeakSpectrum total_loss_spectrum;

        TheoreticalSpectrumGenerator tmp_generator;
        Param new_param(partial_loss_spectrum_generator.getParameters());
        new_param.setValue("add_all_precursor_charges", "true");
        new_param.setValue("add_abundant_immonium_ions", "true");
        new_param.setValue("add_losses", "true");
        tmp_generator.setParameters(new_param);
        tmp_generator.getSpectrum(total_loss_spectrum, fixed_and_variable_modified_peptide, 1, precursor_charge);

        // add special immonium ions
        RNPxlFragmentIonGenerator::addSpecialLysImmonumIons(
          unmodified_sequence,
          total_loss_spectrum, 
          total_loss_spectrum.getIntegerDataArrays()[0],
          total_loss_spectrum.getStringDataArrays()[0]);
        total_loss_spectrum.sortByPosition(); // need to resort after adding special immonium ions

        PeakSpectrum partial_loss_spectrum, partial_loss_template_z1, partial_loss_template_z2, partial_loss_template_z3;
       
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2); 
        partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3); 
        RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                    fixed_and_variable_modified_peptide_weight,
                                    precursor_rna_adduct,
                                    precursor_rna_weight,
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
           partial_loss_spectrum.getIntegerDataArrays()[0],
           partial_loss_spectrum.getStringDataArrays()[0]);

        partial_loss_spectrum.sortByPosition(); // need to resort after adding marker ions

        // fill annotated spectrum information
        set<Size> peak_is_annotated;  // experimental peak index

        // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)
        using MapIonIndexToFragmentAnnotation = map<Size, vector<RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_> >;
        MapIonIndexToFragmentAnnotation unshifted_b_ions, unshifted_y_ions, unshifted_a_ions, shifted_b_ions, shifted_y_ions, shifted_a_ions;
        vector<PeptideHit::PeakAnnotation> shifted_immonium_ions;
        vector<PeptideHit::PeakAnnotation> unshifted_loss_ions;
        vector<PeptideHit::PeakAnnotation> annotated_marker_ions;
        vector<PeptideHit::PeakAnnotation> annotated_precursor_ions;
        vector<PeptideHit::PeakAnnotation> annotated_immonium_ions;

        // first annotate total loss peaks (these give no information where the actual shift occured)
        #ifdef DEBUG_RNPXLSEARCH
          LOG_DEBUG << "Annotating ion (total loss spectrum): " << fixed_and_variable_modified_peptide.toString()  << endl;
        #endif
        vector<pair<Size, Size>> alignment;
        spectrum_aligner.getSpectrumAlignment(alignment, total_loss_spectrum, exp_spectrum);

        const PeakSpectrum::StringDataArray& total_loss_annotations = total_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& total_loss_charges = total_loss_spectrum.getIntegerDataArrays()[0];

        for (auto const & aligned : alignment)
        {
          // information on the experimental fragment in the alignment
          const Size& fragment_index = aligned.second;
          const Peak1D& fragment = exp_spectrum[fragment_index];
          const double fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double fragment_mz = fragment.getMZ();
          const int fragment_charge = exp_spectrum.getIntegerDataArrays().back()[fragment_index];

          const String& ion_name = total_loss_annotations[aligned.first];
          const int charge = total_loss_charges[aligned.first];

          // define which ion names are annotated
          if (ion_name[0] == 'y')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end

            if (loss_first != string::npos) // ion with neutral loss e.g. water
            {
              // only allow matching charges (if a fragment charge was assigned)
              if (fragment_charge == 0 || fragment_charge == charge)
              {
                PeptideHit::PeakAnnotation fa;
                fa.mz = fragment_mz;
                fa.intensity = fragment_intensity;
                fa.charge = charge;
                fa.annotation = ion_name;
                unshifted_loss_ions.push_back(fa);
                peak_is_annotated.insert(aligned.second);
              }
            }
            else
            {
            String ion_nr_string = ion_name.substr(1, charge_pos - 1);
            Size ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_RNPXLSEARCH
              const AASequence& peptide_sequence = fixed_and_variable_modified_peptide.getSuffix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_y_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
            }
          }
          else if (ion_name[0] == 'b')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end

            if (loss_first != string::npos)
            {
              // only allow matching charges (if a fragment charge was assigned)
              if (fragment_charge == 0 || fragment_charge == charge)
              {
                PeptideHit::PeakAnnotation fa;
                fa.mz = fragment_mz;
                fa.intensity = fragment_intensity;
                fa.charge = charge;
                fa.annotation = ion_name;
                unshifted_loss_ions.push_back(fa);
                peak_is_annotated.insert(aligned.second);
              }
            }
            else
            {
              String ion_nr_string = ion_name.substr(1, charge_pos - 1);
              Size ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_RNPXLSEARCH
              const AASequence& peptide_sequence = aas.getPrefix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_b_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
            }
          }
          else if (ion_name[0] == 'a')
          {
            Size loss_first = ion_name.find_first_of('-'); // start of loss
            Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end

            if (loss_first != string::npos)
            {
              // only allow matching charges (if a fragment charge was assigned)
              if (fragment_charge == 0 || fragment_charge == charge)
              {
                PeptideHit::PeakAnnotation fa;
                fa.mz = fragment_mz;
                fa.intensity = fragment_intensity;
                fa.charge = charge;
                fa.annotation = ion_name;
                unshifted_loss_ions.push_back(fa);
                peak_is_annotated.insert(aligned.second);
              }
            }
            else
            {
              String ion_nr_string = ion_name.substr(1, charge_pos - 1);
              auto ion_number = (Size)ion_nr_string.toInt();
            #ifdef DEBUG_RNPXLSEARCH
              const AASequence& peptide_sequence = aas.getPrefix(ion_number);
              LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
            #endif
            peak_is_annotated.insert(aligned.second);

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
              unshifted_a_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
            }
          }
          else if (ion_name.hasPrefix("[M+")) // precursor ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1; // for visualion charge is not really important so we set it to 0
            fa.annotation = ion_name;
            peak_is_annotated.insert(aligned.second);
            annotated_precursor_ions.push_back(fa);
          }
          else if (ion_name.hasPrefix("i")) // immonium ion
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1;
            fa.annotation = ion_name;
            peak_is_annotated.insert(aligned.second);
            annotated_immonium_ions.push_back(fa);
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

        spectrum_aligner.getSpectrumAlignment(alignment, partial_loss_spectrum, exp_spectrum);

        const PeakSpectrum::StringDataArray& partial_loss_annotations = partial_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& partial_loss_charges = partial_loss_spectrum.getIntegerDataArrays()[0];

        if (alignment.empty())
        {
          a.fragment_annotations = fas;
          continue;
        }

        for (vector<std::pair<Size, Size> >::const_iterator pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
        {
          // only annotate experimental peaks with shift - i.e. do not annotated complete loss peaks again
          if (peak_is_annotated.find(pair_it->second) != peak_is_annotated.end()) { continue; }

          // information on the experimental fragment in the alignment
          const Size & fragment_index = pair_it->second;
          const Peak1D & fragment = exp_spectrum[fragment_index];
          const double & fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double & fragment_mz = fragment.getMZ();
          const int & fragment_charge = exp_spectrum.getIntegerDataArrays().back()[fragment_index];

          String ion_name = partial_loss_annotations[pair_it->first];
          const int charge = partial_loss_charges[pair_it->first];

          vector<String> f;

          ion_name.split(' ', f);  // e.g. "y3 C3O" or just "y2"
          String fragment_shift_name;
          if (f.size() == 2) { fragment_shift_name = f[1]; }

          String fragment_ion_name = f[0]; // e.g. y3

          #ifdef DEBUG_RNPXLSEARCH
            LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
          #endif

          // define which ion names are annotated
          if (fragment_ion_name.hasPrefix("y"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("y", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_y_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (fragment_ion_name.hasPrefix("b"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("b", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_b_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (fragment_ion_name.hasPrefix("a"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("a", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            // only allow matching charges (if a fragment charge was assigned)
            if (fragment_charge == 0 || fragment_charge == charge)
            {
              RNPxlFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
              shifted_a_ions[ion_number].push_back(d);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << charge << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix(ANNOTATIONS_MARKER_ION_PREFIX))
          {
            if (fragment_charge <= 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              annotated_marker_ions.push_back(fa);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << 1 << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("i"))
          {
            if (fragment_charge <= 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              shifted_immonium_ions.push_back(fa);
            }
            #ifdef DEBUG_RNPXLSEARCH
            else
            {
              LOG_DEBUG << "Charge missmatch in alignment: " << ion_name << " at position: " << fragment_mz << " charge fragment: " << fragment_charge << " theo. charge: " << 1 << endl;
            }
            #endif
          }
          else if (ion_name.hasPrefix("[M+"))
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = 1; // for visualion charge is not really important so we set it to 1, TODO: read out charge from ion name and set here
            fa.annotation = ion_name;
            annotated_precursor_ions.push_back(fa);
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

#ifdef DEBUG_RNPXLSEARCH
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
#ifdef DEBUG_RNPXLSEARCH
        cout << "site sum score (shifted a/b/y-ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        #ifdef DEBUG_RNPXLSEARCH
          LOG_DEBUG << "Localisation based on immonium ions: ";
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
#ifdef DEBUG_RNPXLSEARCH
        cout << "site sum score (shifted a/b/y-ions & immonium ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        String best_localization = unmodified_sequence;
        double best_localization_score = 0;
        String localization_scores;
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          if (sites_sum_score[i] > best_localization_score) { best_localization_score = sites_sum_score[i]; }
        }

        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          #ifdef DEBUG_RNPXLSEARCH
            LOG_DEBUG << String::number(100.0 * sites_sum_score[i], 2);
          #endif

          if (i != 0) localization_scores += ' ';
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
          }
        }
        #ifdef DEBUG_RNPXLSEARCH
          LOG_DEBUG << endl;
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
        a.fragment_annotations = fas;

        #ifdef DEBUG_RNPXLSEARCH
          LOG_DEBUG << "Ion centric annotation: " << endl;
          LOG_DEBUG << "unshifted b ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", unshifted_b_ions) << endl;
          LOG_DEBUG << "unshifted y ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", unshifted_y_ions) << endl;
          LOG_DEBUG << "unshifted a ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", unshifted_a_ions) << endl;
          LOG_DEBUG << "shifted b ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", shifted_b_ions) << endl;
          LOG_DEBUG << "shifted y ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", shifted_y_ions) << endl;
          LOG_DEBUG << "shifted a ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", shifted_a_ions) << endl;
          LOG_DEBUG << "shifted immonium ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(shifted_immonium_ions) << endl;
          LOG_DEBUG << "shifted marker ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(annotated_marker_ions) << endl;
          LOG_DEBUG << "shifted precursor ions: " << endl;
          LOG_DEBUG << RNPxlFragmentAnnotationHelper::shiftedIonsToString(annotated_precursor_ions) << endl;
          LOG_DEBUG << "Localization scores: ";
          LOG_DEBUG << localization_scores << endl;
          LOG_DEBUG << "Localisation based on ion series and immonium ions of all observed fragments: ";
          LOG_DEBUG << best_localization << endl;
        #endif
      }
    }
  }

  /// Filter by top scoring hits, reconstruct original peptide from memory efficient structure, and add additional meta information.
  void postProcessHits_(const PeakMap& exp, 
    vector<vector<AnnotatedHit> >& annotated_hits, 
    vector<ProteinIdentification>& protein_ids, 
    vector<PeptideIdentification>& peptide_ids, 
    Size top_hits, 
    const RNPxlModificationMassesResult& mm, 
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide,
    const vector<PrecursorPurity::PurityScores>& purities)
  {
      // remove all but top n scoring (Note: this is currently necessary as postScoreHits_ might reintroduce nucleotide specific hits for fast scoring)
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
      {
        // sort and keeps n best elements according to score
        Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
        std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
        annotated_hits[scan_index].resize(topn);
        annotated_hits.shrink_to_fit();
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi;
        pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        pi.setScoreType("hyperscore");
        pi.setHigherScoreBetter(true);
        pi.setRT(exp[scan_index].getRT());
        pi.setMZ(exp[scan_index].getPrecursors()[0].getMZ());
        pi.setMetaValue("precursor_intensity", exp[scan_index].getPrecursors()[0].getIntensity());
        Size charge = exp[scan_index].getPrecursors()[0].getCharge();

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        size_t rank(0);
        for (auto const & ah : annotated_hits[scan_index])
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
          ph.setMetaValue(String("RNPxl:score"), ah.score); // important for Percolator feature set because the PeptideHit score might be overwritten by a q-value

          // determine RNA modification from index in map
          std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, ah.rna_mod_index);
          ph.setMetaValue(String("RNPxl:total_loss_score"), ah.total_loss_score);
          ph.setMetaValue(String("RNPxl:immonium_score"), ah.immonium_score);
          ph.setMetaValue(String("RNPxl:precursor_score"), ah.precursor_score);
          ph.setMetaValue(String("RNPxl:a_ion_score"), ah.a_ion_score);
          ph.setMetaValue(String("RNPxl:marker_ions_score"), ah.marker_ions_score);
          ph.setMetaValue(String("RNPxl:partial_loss_score"), ah.partial_loss_score);

          // total loss and partial loss (pl) related subscores (matched ion current, avg. fragment error, morpheus score)
          ph.setMetaValue(String("RNPxl:MIC"), ah.MIC);
          ph.setMetaValue(String("RNPxl:err"), ah.err);
          ph.setMetaValue(String("RNPxl:Morph"), ah.Morph);
          ph.setMetaValue(String("RNPxl:pl_MIC"), ah.pl_MIC);
          ph.setMetaValue(String("RNPxl:pl_err"), ah.pl_err);
          ph.setMetaValue(String("RNPxl:pl_Morph"), ah.pl_Morph);
          ph.setMetaValue(String("RNPxl:total_MIC"), ah.total_MIC);  // fraction of matched ion current from total + partial losses

          ph.setMetaValue(String("RNPxl:RNA"), *mod_combinations_it->second.begin()); // return first nucleotide formula matching the index of the empirical formula
          ph.setMetaValue(String("RNPxl:NT"), String(ah.cross_linked_nucleotide));  // the cross-linked nucleotide
          ph.setMetaValue(String("RNPxl:RNA_MASS_z0"), EmpiricalFormula(mod_combinations_it->first).getMonoWeight()); // RNA uncharged mass via empirical formula

          ph.setMetaValue(String("RNPxl:best_localization_score"), ah.best_localization_score);
          ph.setMetaValue(String("RNPxl:localization_scores"), ah.localization_scores);
          ph.setMetaValue(String("RNPxl:best_localization"), ah.best_localization);

          if (!purities.empty())
          {
            ph.setMetaValue("precursor_purity", purities[scan_index].signal_proportion);
          }

          ph.setPeakAnnotations(ah.fragment_annotations);
          ph.setMetaValue("isotope_error", static_cast<int>(ah.isotope_error));
          ph.setMetaValue("rank", rank);
          // set the amino acid sequence (for complete loss spectra this is just the variable and modified peptide. For partial loss spectra it additionally contains the loss induced modification)
          ph.setSequence(fixed_and_variable_modified_peptide);
          phs.push_back(ph);
          ++rank;
        }

        pi.setHits(phs);
        pi.assignRanks();

#ifdef _OPENMP
#pragma omp critical (peptide_ids_access)
#endif
        {
          peptide_ids.push_back(pi);
        }
      }
    }

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("RNPxlSearch");
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

    /* default features added in PercolatorAdapter:
     * SpecId, ScanNr, ExpMass, CalcMass, mass, 
     * peplen, charge#min..#max, enzN, enzC, enzInt, dm, absdm
     */     
    StringList feature_set;
    feature_set
       << "isotope_error"
       << "RNPxl:score"
       << "RNPxl:total_loss_score"
       << "RNPxl:immonium_score"
       << "RNPxl:precursor_score"
       << "RNPxl:a_ion_score"
       << "RNPxl:marker_ions_score"
       << "RNPxl:partial_loss_score"
       << "RNPxl:MIC"
       << "RNPxl:err"
       << "RNPxl:Morph"
       << "RNPxl:pl_MIC"
       << "RNPxl:pl_err"
       << "RNPxl:pl_Morph"
       << "RNPxl:total_MIC"
       << "RNPxl:RNA_MASS_z0";
    if (!purities.empty()) feature_set << "precursor_purity";

    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));

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

  void initializeSpectrumGenerators(TheoreticalSpectrumGenerator &total_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &partial_loss_spectrum_generator,
                                  TheoreticalSpectrumGenerator &a_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &immonium_ion_sub_score_spectrum_generator,
                                  TheoreticalSpectrumGenerator &precursor_ion_sub_score_spectrum_generator) const
  {
    Param param(total_loss_spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    total_loss_spectrum_generator.setParameters(param);

    param = partial_loss_spectrum_generator.getParameters();
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

    // generator for sub scores for a-ions, immonium and precursor peaksparam = a_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    a_ion_sub_score_spectrum_generator.setParameters(param);
    param = immonium_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "true");
    param.setValue("add_precursor_peaks", "false");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    immonium_ion_sub_score_spectrum_generator.setParameters(param);
    param = precursor_ion_sub_score_spectrum_generator.getParameters();
    param.setValue("add_abundant_immonium_ions", "false");
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "false");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "false");
    param.setValue("add_z_ions", "false");
    param.setValue("add_metainfo", "true");
    precursor_ion_sub_score_spectrum_generator.setParameters(param);
  }

  ExitCodes main_(int, const char**) override
  {
    // force initialization of residue db
    AASequence::fromString("GPAVLIMCFYWHKRQNEDST");

    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in_mzml = getStringOption_("in");
    String in_db = getStringOption_("database");
    String out_idxml = getStringOption_("out");
    String out_csv = getStringOption_("out_tsv");

    fast_scoring_ = getStringOption_("RNPxl:scoring") == "fast" ? true : false;

    bool generate_decoys = getFlag_("RNPxl:decoys");

    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    IntList precursor_isotopes = getIntList_("precursor:isotopes");

    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    double marker_ions_tolerance = getDoubleOption_("RNPxl:marker_ions_tolerance");

    double small_peptide_mass_filter_threshold = getDoubleOption_("RNPxl:filter_small_peptide_mass");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    if (fixed_unique.size() != fixedModNames.size())
    {
      LOG_WARN << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      LOG_WARN << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(fixedModNames);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    size_t report_top_hits = (size_t)getIntOption_("report:top_hits");

    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides = getStringList_("RNPxl:target_nucleotides");

    StringList nt_groups = getStringList_("RNPxl:nt_groups");

    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings = getStringList_("RNPxl:mapping");

    // read list of nucleotides that can directly cross-link
    // these are responsible for shifted fragment ions. Their fragment adducts thus determine which shifts will be observed on b-,a-,y-ions
    String can_cross_link = getStringOption_("RNPxl:can_cross_link");
    for (const auto& c : can_cross_link) { can_xl_.insert(c); }

    StringList modifications = getStringList_("RNPxl:modifications");

    String sequence_restriction = getStringOption_("RNPxl:sequence");

    Int max_nucleotide_length = getIntOption_("RNPxl:length");

    bool cysteine_adduct = getFlag_("RNPxl:CysteineAdduct");

    // generate all precursor adducts
    RNPxlModificationMassesResult mm;
    if (max_nucleotide_length != 0)
    {
      mm = RNPxlModificationsGenerator::initModificationMassesRNA(
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
      mm.mod_masses[""] = 0; // insert "null" modification otherwise peptides without RNA will not be searched
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
    p.setValue("add_decoy_peptides", "true"); // we still want decoys in the result (e.g., to run percolator)
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
    vector<PrecursorPurity::PurityScores> purities;
    {
      PeakMap tmp_spectra;
      f.load(in_mzml, tmp_spectra);
      int nMS1 = std::count_if(tmp_spectra.begin(), tmp_spectra.end(), [](MSSpectrum& s){return s.getMSLevel() == 1;});
      if (nMS1 != 0)
      {
        purities = PrecursorPurity::computePrecursorPurities(tmp_spectra, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
      }
    } // free spectra  

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);

    progresslogger.startProgress(0, 1, "Filtering spectra...");
    const bool convert_to_single_charge = false;  // whether to convert fragment peaks with isotopic patterns to single charge
    const bool annotate_charge = false;  // whether the charge and type is annotated
    preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, convert_to_single_charge, annotate_charge);
    progresslogger.endProgress();

    // build multimap of precursor mass to scan index (and perform some mass and length based filtering)
    using MassToScanMultiMap = multimap<double, pair<Size, int>>;
    MassToScanMultiMap multimap_mass_2_scan_index;  // map precursor mass to scan index and (potential) isotopic missassignment
    mapPrecursorMassesToScans(min_precursor_charge,
                              max_precursor_charge,
                              precursor_isotopes,
                              small_peptide_mass_filter_threshold,
                              peptide_min_size,
                              spectra,
                              multimap_mass_2_scan_index);

    // initialize spectrum generators (generated ions, etc.)
    TheoreticalSpectrumGenerator total_loss_spectrum_generator;
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    TheoreticalSpectrumGenerator a_ion_sub_score_spectrum_generator;
    TheoreticalSpectrumGenerator immonium_ion_sub_score_spectrum_generator;
    TheoreticalSpectrumGenerator precursor_ion_sub_score_spectrum_generator;
    initializeSpectrumGenerators(total_loss_spectrum_generator,
                                 partial_loss_spectrum_generator,
                                 a_ion_sub_score_spectrum_generator,
                                 immonium_ion_sub_score_spectrum_generator,
                                 precursor_ion_sub_score_spectrum_generator);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit> > annotated_hits(spectra.size(), vector<AnnotatedHit>());
    for (auto & a : annotated_hits) { a.reserve(2 * report_top_hits); }

#ifdef _OPENMP     
    // we want to do locking at the spectrum level so we get good parallelisation 
    vector<omp_lock_t> annotated_hits_lock(annotated_hits.size());
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_init_lock(&(annotated_hits_lock[i])); }
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
      progresslogger.startProgress(0, 1, "Generate decoys...");

      // append decoy proteins
      const size_t old_size = fasta_db.size();
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = fasta_db[i];
        e.sequence.reverse();
        e.identifier = "DECOY_" + e.identifier;
        fasta_db.push_back(e);
      }
      progresslogger.endProgress();
    }

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
        }

        if (already_processed) { continue; }

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          processed_petides.insert(*cit);
        }

#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;
        vector<AASequence> all_modified_peptides;

        const String unmodified_sequence = cit->getString();

#ifdef _OPENMP
#pragma omp critical (residuedb_access)
#endif
        {
           // only process peptides without ambiguous amino acids (placeholder / any amino acid)
          if (unmodified_sequence.find_first_of("XBZ") == std::string::npos)
          {
            AASequence aas = AASequence::fromString(unmodified_sequence);
            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);
          }
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass_without_RNA = fixed_and_variable_modified_peptide.getMonoWeight();

          //create empty theoretical spectrum.  total_loss_spectrum_z2 contains both charge 1 and charge 2 peaks
          PeakSpectrum total_loss_spectrum_z1, total_loss_spectrum_z2;

          // spectrum containing additional peaks for sub scoring
          PeakSpectrum immonium_sub_score_spectrum, 
                       a_ion_sub_score_spectrum, 
                       precursor_sub_score_spectrum,
                       marker_ions_sub_score_spectrum;

          // iterate over all RNA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          Size rna_mod_index = 0;

          // TODO: track the XL-able nt here
          for (std::map<String, double>::const_iterator rna_mod_it = mm.mod_masses.begin(); rna_mod_it != mm.mod_masses.end(); ++rna_mod_it, ++rna_mod_index)
          {            
            const double precursor_rna_weight = rna_mod_it->second;
            const double current_peptide_mass = current_peptide_mass_without_RNA + precursor_rna_weight; // add RNA mass
            // TODO: const char xl_nucleotide; // can be none

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
            if (total_loss_spectrum_z1.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              total_loss_spectrum_generator.getSpectrum(total_loss_spectrum_z1, fixed_and_variable_modified_peptide, 1, 1);
              total_loss_spectrum_generator.getSpectrum(total_loss_spectrum_z2, fixed_and_variable_modified_peptide, 1, 2);
              immonium_ion_sub_score_spectrum_generator.getSpectrum(immonium_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
              RNPxlFragmentIonGenerator::addSpecialLysImmonumIons(
                unmodified_sequence, 
                immonium_sub_score_spectrum, 
                immonium_sub_score_spectrum.getIntegerDataArrays()[0], 
                immonium_sub_score_spectrum.getStringDataArrays()[0]);
              immonium_sub_score_spectrum.sortByPosition();
              precursor_ion_sub_score_spectrum_generator.getSpectrum(precursor_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
              a_ion_sub_score_spectrum_generator.getSpectrum(a_ion_sub_score_spectrum, fixed_and_variable_modified_peptide, 1, 1);
            }

            if (!fast_scoring_)
            {
              PeakSpectrum marker_ions_sub_score_spectrum_z1;
              //shifted_immonium_ions_sub_score_spectrum;
              PeakSpectrum partial_loss_spectrum_z1, partial_loss_spectrum_z2;

              // retrieve RNA adduct name
              auto mod_combinations_it = mm.mod_combinations.begin();
              std::advance(mod_combinations_it, rna_mod_index);
              const String& precursor_rna_adduct = *mod_combinations_it->second.begin();

              if (precursor_rna_adduct == "none")
              {
                // score peptide without RNA (same method as fast scoring)
                for (auto l = low_it; l != up_it; ++l) // OMS_CODING_TEST_EXCLUDE
                {
                  //const double exp_pc_mass = l->first;
                  const Size & scan_index = l->second.first;
                  const int & isotope_error = l->second.second;
                  const PeakSpectrum & exp_spectrum = spectra[scan_index];
                  const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                  PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                  float total_loss_score(0), 
                        immonium_sub_score(0), 
                        precursor_sub_score(0), 
                        a_ion_sub_score(0), 
                        tlss_MIC(0),
                        tlss_err(0), 
                        tlss_Morph(0);

                  scoreTotalLossFragments_(exp_spectrum,
                                         total_loss_spectrum,
                                         fragment_mass_tolerance,
                                         fragment_mass_tolerance_unit_ppm,
                                         a_ion_sub_score_spectrum,
                                         precursor_sub_score_spectrum,
                                         immonium_sub_score_spectrum,
                                         total_loss_score,
                                         tlss_MIC,
                                         tlss_err,
                                         tlss_Morph,
                                         immonium_sub_score,
                                         precursor_sub_score,
                                         a_ion_sub_score);


                  // bad score, likely wihout any single matching peak
                  if (total_loss_score < 0.01) { continue; }

                  // add peptide hit
                  AnnotatedHit ah;
                  ah.sequence = *cit; // copy StringView
                  ah.peptide_mod_index = mod_pep_idx;
                  ah.MIC = tlss_MIC;
                  ah.err = tlss_err;
                  ah.Morph = tlss_Morph;
                  ah.total_loss_score = total_loss_score;
                  ah.immonium_score = immonium_sub_score;
                  ah.precursor_score = precursor_sub_score;
                  ah.a_ion_score = a_ion_sub_score;
                  ah.total_MIC = tlss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                  ah.rna_mod_index = rna_mod_index;
                  ah.isotope_error = isotope_error;

                  // combined score
                  ah.score = RNPxlSearch::calculateCombinedScore(ah, false);

#ifdef DEBUG_RNPXLSEARCH
                  LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP 
                  omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                  {
                    annotated_hits[scan_index].emplace_back(move(ah));

                    // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                    if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                    {
                      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
                      annotated_hits[scan_index].resize(report_top_hits); 
                    }
                  }
#ifdef _OPENMP 
                  omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
                }
              }
              else  // score peptide with RNA adduct
              {
                PeakSpectrum partial_loss_template_z1, partial_loss_template_z2, partial_loss_template_z3;
                partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1); 
                partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2); 
                partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3); 

                // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
                // get RNA fragment shifts in the MS2 (based on the precursor RNA/DNA)
                auto const & all_NA_adducts = all_feasible_fragment_adducts.at(precursor_rna_adduct);
                const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_NA_adducts.feasible_adducts;
                // get marker ions
                const vector<FragmentAdductDefinition_>& marker_ions = all_NA_adducts.marker_ions;

                //cout << "'" << precursor_rna_adduct << "'" << endl;
                //OPENMS_POSTCONDITION(!feasible_MS2_adducts.empty(),
                //                String("FATAL: No feasible adducts for " + precursor_rna_adduct).c_str());


                // Do we have (nucleotide) specific fragmentation adducts? for the current RNA adduct on the precursor?
                // If so, generate spectra for shifted ion series

                // score individually for every nucleotide
                for (auto const & nuc_2_adducts : feasible_MS2_adducts)
                {
                  const char& cross_linked_nucleotide = nuc_2_adducts.first;
                  const vector<FragmentAdductDefinition_>& partial_loss_modification = nuc_2_adducts.second;

                  if (!partial_loss_modification.empty())
                  {
                    // shifted b- / y- / a-ions
                    // generate shifted_immonium_ions_sub_score_spectrum.empty
                    RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                                current_peptide_mass_without_RNA,
                                                precursor_rna_adduct,
                                                precursor_rna_weight,
                                                1,
                                                partial_loss_modification,
					        partial_loss_template_z1,
					        partial_loss_template_z2,
                                                partial_loss_template_z3,
                                                partial_loss_spectrum_z1);
                    for (auto& n : partial_loss_spectrum_z1.getStringDataArrays()[0]) { n[0] = 'y'; } // hyperscore hack

                    RNPxlFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                                current_peptide_mass_without_RNA,
                                                precursor_rna_adduct,
                                                precursor_rna_weight,
                                                2, // don't know the charge of the precursor at that point
                                                partial_loss_modification,
					        partial_loss_template_z1,
					        partial_loss_template_z2,
                                                partial_loss_template_z3,
                                                partial_loss_spectrum_z2);
                    for (auto& n : partial_loss_spectrum_z2.getStringDataArrays()[0]) { n[0] = 'y'; } // hyperscore hack
                  }

                  // add shifted marker ions
                  marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
                  marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation
                  RNPxlFragmentIonGenerator::addMS2MarkerIons(
                    marker_ions,
                    marker_ions_sub_score_spectrum_z1,
                    marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[0],
                    marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

                  for (auto l = low_it; l != up_it; ++l) // OMS_CODING_TEST_EXCLUDE
                  {
                    //const double exp_pc_mass = l->first;
                    const Size& scan_index = l->second.first;
                    const int& isotope_error = l->second.second;
                    const PeakSpectrum& exp_spectrum = spectra[scan_index];
                    float tlss_MIC(0), tlss_err(0), tlss_Morph(0),
                      immonium_sub_score(0), precursor_sub_score(0),
                      a_ion_sub_score(0), partial_loss_sub_score(0), marker_ions_sub_score(0),
                      plss_MIC(0), plss_err(0), plss_Morph(0), score;

                    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                    PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                    scoreTotalLossFragments_(exp_spectrum,
                                             total_loss_spectrum,
                                             fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                             a_ion_sub_score_spectrum,
                                             precursor_sub_score_spectrum,
                                             immonium_sub_score_spectrum,
                                             score,
                                             tlss_MIC,
                                             tlss_err,
                                             tlss_Morph,
                                             immonium_sub_score,
                                             precursor_sub_score,
                                             a_ion_sub_score);

                    // bad score, likely wihout any single matching peak
                    if (score < 0.01) { continue; }

                    scorePartialLossFragments_(exp_spectrum,
                                               fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                               partial_loss_spectrum_z1, partial_loss_spectrum_z2,
                                               marker_ions_sub_score_spectrum_z1,
                                               partial_loss_sub_score,
                                               marker_ions_sub_score,
                                               plss_MIC, plss_err, plss_Morph);

                    // add peptide hit
                    AnnotatedHit ah;
                    ah.sequence = *cit; // copy StringView
                    ah.peptide_mod_index = mod_pep_idx;
                    ah.total_loss_score = score;
                    ah.MIC = tlss_MIC;
                    ah.err = tlss_err;
                    ah.Morph = tlss_Morph;
                    ah.pl_MIC = plss_MIC;
                    ah.pl_err = plss_err;
                    ah.pl_Morph = plss_Morph;
                    ah.immonium_score = immonium_sub_score;
                    ah.precursor_score = precursor_sub_score;
                    ah.a_ion_score = a_ion_sub_score;
                    ah.cross_linked_nucleotide = cross_linked_nucleotide;
                    ah.total_MIC = tlss_MIC + plss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                    // scores from shifted peaks
                    ah.marker_ions_score = marker_ions_sub_score;
                    ah.partial_loss_score = partial_loss_sub_score;

                    ah.rna_mod_index = rna_mod_index;
                    ah.isotope_error = isotope_error;

                    // combined score
                    ah.score = RNPxlSearch::calculateCombinedScore(ah, true);

#ifdef DEBUG_RNPXLSEARCH
                    LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
                    omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                    {
                      annotated_hits[scan_index].emplace_back(move(ah));

                      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                      if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                      {
                        std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
                        annotated_hits[scan_index].resize(report_top_hits); 
                      }
                    }
#ifdef _OPENMP
                    omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
                  }
                } // for every nucleotide in the precursor
              }
            }
            else // fast scoring
            {
              for (auto l = low_it; l != up_it; ++l) // OMS_CODING_TEST_EXCLUDE
              {
                //const double exp_pc_mass = l->first;
                const Size &scan_index = l->second.first;
                const int &isotope_error = l->second.second;
                const PeakSpectrum &exp_spectrum = spectra[scan_index];
                float total_loss_score;
                float immonium_sub_score;
                float precursor_sub_score;
                float a_ion_sub_score;
                float tlss_MIC;
                float tlss_err;
                float tlss_Morph;

                const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
                PeakSpectrum & total_loss_spectrum = (exp_pc_charge < 3) ? total_loss_spectrum_z1 : total_loss_spectrum_z2;

                scoreTotalLossFragments_(exp_spectrum, 
                                         total_loss_spectrum, 
                                         fragment_mass_tolerance,
                                         fragment_mass_tolerance_unit_ppm, 
                                         a_ion_sub_score_spectrum,
                                         precursor_sub_score_spectrum, 
                                         immonium_sub_score_spectrum, 
                                         total_loss_score, 
                                         tlss_MIC,
                                         tlss_err,
                                         tlss_Morph,
                                         immonium_sub_score,
                                         precursor_sub_score,
                                         a_ion_sub_score);

                // no good hit
                if (total_loss_score < 0.01) { continue; }

                // add peptide hit
                AnnotatedHit ah;
                ah.sequence = *cit; // copy StringView
                ah.peptide_mod_index = mod_pep_idx;
                ah.total_loss_score = total_loss_score;
                ah.MIC = tlss_MIC;
                ah.err = tlss_err;
                ah.Morph = tlss_Morph;
                ah.immonium_score = immonium_sub_score;
                ah.precursor_score = precursor_sub_score;
                ah.a_ion_score = a_ion_sub_score;

                ah.total_MIC = tlss_MIC + immonium_sub_score + a_ion_sub_score + precursor_sub_score;

                ah.rna_mod_index = rna_mod_index;
                ah.isotope_error = isotope_error;

                // simple combined score in fast scoring:
                ah.score = total_loss_score + ah.total_MIC; 

#ifdef DEBUG_RNPXLSEARCH
                LOG_DEBUG << "best score in pre-score: " << score << endl;
#endif

#ifdef _OPENMP
                omp_set_lock(&(annotated_hits_lock[scan_index]));
#endif
                {
                  annotated_hits[scan_index].emplace_back(move(ah));

                  // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                  if (annotated_hits[scan_index].size() >= 2 * report_top_hits)
                  {
                    std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + report_top_hits, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
                    annotated_hits[scan_index].resize(report_top_hits); 
                  }
                }
#ifdef _OPENMP
                omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    LOG_INFO << "Proteins: " << count_proteins << endl;
    LOG_INFO << "Peptides: " << count_peptides << endl;
    LOG_INFO << "Processed peptides: " << processed_petides.size() << endl;

    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;
    progresslogger.startProgress(0, 1, "Post-processing PSMs...");

    // Localization
    //

    // reload spectra from disc with same settings as before (important to keep same spectrum indices)
    spectra.clear(true);
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);    

    // for post scoring don't convert fragments to single charge. Annotate charge instead to every peak.
    preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, true); // no single charge (false), annotate charge (true)

    progresslogger.startProgress(0, 1, "localization...");

    // create spectrum generator. For convenience we add more peak types here.
    Param param(total_loss_spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "true");
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    total_loss_spectrum_generator.setParameters(param);

    postScoreHits_(spectra, 
                   annotated_hits, 
                   report_top_hits, 
                   mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, 
                   partial_loss_spectrum_generator, 
                   fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
                   all_feasible_fragment_adducts);

    progresslogger.startProgress(0, 1, "Post-processing and annotation...");
    postProcessHits_(spectra, 
                     annotated_hits, 
                     protein_ids, peptide_ids, 
                     report_top_hits, 
                     mm, 
                     fixed_modifications, variable_modifications, 
                     max_variable_mods_per_peptide,
                     purities);
    progresslogger.endProgress();

    // reindex ids
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", getStringOption_("peptide:enzyme"));
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    param_pi.setValue("log", getStringOption_("log"));
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

    // annotate RNPxl related information to hits and create report
    vector<RNPxlReportRow> csv_rows = RNPxlReport::annotate(spectra, peptide_ids, marker_ions_tolerance);

    if (generate_decoys)	
    {
      fdr.apply(peptide_ids);	
    }

    // write ProteinIdentifications and PeptideIdentifications to IdXML
    IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

    // save report
    if (!out_csv.empty())
    {
      TextFile csv_file;
      csv_file.addLine(RNPxlReportRowHeader().getString("\t"));
      for (Size i = 0; i != csv_rows.size(); ++i)
      {
        csv_file.addLine(csv_rows[i].getString("\t"));
      }
      csv_file.store(out_csv);
    }
 
 #ifdef _OPENMP
    // free locks
    for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_destroy_lock(&(annotated_hits_lock[i])); }
 #endif

    return EXECUTION_OK;
  }

  // determine main score and sub scores of peaks without shifts
  void scoreTotalLossFragments_(const PeakSpectrum &exp_spectrum,
                                const PeakSpectrum &total_loss_spectrum,
                                double fragment_mass_tolerance,
                                bool fragment_mass_tolerance_unit_ppm,
                                const PeakSpectrum &a_ion_sub_score_spectrum,
                                const PeakSpectrum &precursor_sub_score_spectrum,
                                const PeakSpectrum &immonium_sub_score_spectrum,
                                float &total_loss_score,
                                float &tlss_MIC,
                                float &tlss_err,
                                float &tlss_Morph,
                                float &immonium_sub_score,
                                float &precursor_sub_score,
                                float &a_ion_sub_score) const
  {
    total_loss_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                           exp_spectrum, total_loss_spectrum);

    // bad score, likely wihout any single matching peak
    if (total_loss_score < 0.01) { return; }

    immonium_sub_score = 0;
    precursor_sub_score = 0;
    a_ion_sub_score = 0;
    auto const & tl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                       fragment_mass_tolerance_unit_ppm,
                                                       exp_spectrum,
                                                       total_loss_spectrum);

    tlss_MIC = tl_sub_scores.TIC != 0 ? tl_sub_scores.MIC / tl_sub_scores.TIC : 0;
    tlss_err = tl_sub_scores.err;
    tlss_Morph = tl_sub_scores.score;

    if (!immonium_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             immonium_sub_score_spectrum);
      immonium_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    if (!precursor_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             precursor_sub_score_spectrum);
      precursor_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    if (!a_ion_sub_score_spectrum.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             a_ion_sub_score_spectrum);
      a_ion_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
  }

  void scorePartialLossFragments_(const PeakSpectrum &exp_spectrum,
                                  double fragment_mass_tolerance,
                                  bool fragment_mass_tolerance_unit_ppm,
                                  const PeakSpectrum &partial_loss_spectrum_z1,
                                  const PeakSpectrum &partial_loss_spectrum_z2,
                                  const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                                  float &partial_loss_sub_score,
                                  float &marker_ions_sub_score,
                                  float &plss_MIC, float &plss_err, float &plss_Morph) const
  {
    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    plss_MIC = 0;
    plss_err = 0;
    plss_Morph = 0;


    if (!marker_ions_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             marker_ions_sub_score_spectrum_z1);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }
    //TODO: these are currently empty
/*
              if (!shifted_immonium_ions_sub_score_spectrum.empty())
              {
                shifted_immonium_ions_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, shifted_immonium_ions_sub_score_spectrum);
                               shifted_immonium_ions_sub_score(0),
              }
*/
    if (!partial_loss_spectrum_z1.empty()) // check if we generated partial loss spectra
    {
      if (exp_pc_charge < 3)
      {
        partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                                     exp_spectrum, partial_loss_spectrum_z1);
        auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                           fragment_mass_tolerance_unit_ppm,
                                                           exp_spectrum,
                                                           partial_loss_spectrum_z1);

        plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
        plss_err = pl_sub_scores.err;
        plss_Morph = pl_sub_scores.score;
      }
      else //if (exp_pc_charge >= 3)
      {
        partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
                                                     exp_spectrum, partial_loss_spectrum_z2);
        auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                           fragment_mass_tolerance_unit_ppm,
                                                           exp_spectrum,
                                                           partial_loss_spectrum_z2);

        plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
        plss_err = pl_sub_scores.err;
        plss_Morph = pl_sub_scores.score;
      }
    }
#ifdef DEBUG_RNPXLSEARCH
    LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
  // TODO: cap plss_err
  float ft_da = fragment_mass_tolerance_unit_ppm ? fragment_mass_tolerance * 1e-6 * 1000.0 : fragment_mass_tolerance;
  if (plss_err > ft_da) plss_err = ft_da;
  }

};

RNPxlSearch::RNPxlParameterParsing::PrecursorsToMS2Adducts
RNPxlSearch::RNPxlParameterParsing::getAllFeasibleFragmentAdducts(
  const RNPxlModificationMassesResult &precursor_adducts,
  const RNPxlSearch::RNPxlParameterParsing::NucleotideToFragmentAdductMap &nucleotide_to_fragment_adducts,
  const set<char> &can_xl) 
{
  PrecursorsToMS2Adducts all_pc_all_feasible_adducts;

  // for all possible precursor adducts
  for (auto const & pa : precursor_adducts.mod_masses)
  {
    const String& ef = pa.first;
    const set<String>& ambiguities = precursor_adducts.mod_combinations.at(ef);

    // for each ambiguous precursor adduct (stored as nucleotide formula e.g.: "AU-H2O")
    for (auto const & pc_adduct : ambiguities)
    {
      // calculate feasible fragment adducts and store them for lookup
      auto feasible_adducts = getFeasibleFragmentAdducts(pc_adduct, ef, nucleotide_to_fragment_adducts, can_xl);
      // TODO: check if needed anymore - std::sort(feasible_adducts.begin(), feasible_adducts.end());
      all_pc_all_feasible_adducts[pc_adduct] = feasible_adducts;
      break; // only store one precursor adduct for multiple ambiguities (e.g. AUG, AGU, UAG..)
    }
  }

  // print feasible fragment adducts and marker ions
  for (auto const & fa : all_pc_all_feasible_adducts)
  {
    LOG_DEBUG << "Precursor adduct: " << fa.first << "\n";

    for (auto const & ffa : fa.second.feasible_adducts)
    {
      const char & nucleotide = ffa.first;
      LOG_DEBUG << "  Cross-linkable nucleotide '" << nucleotide << "' and feasible fragment adducts:" << endl;
      for (auto const & a : ffa.second)
      {
        LOG_DEBUG << "    " << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
      }
    }

    LOG_DEBUG << "  Marker ions." << endl;
    for (auto const & ffa : fa.second.marker_ions)
    {
      LOG_DEBUG << "    "  << ffa.name << "\t" << ffa.formula.toString() << "\t" << ffa.mass << "\n";
    }
  }
  LOG_DEBUG << endl;

  return all_pc_all_feasible_adducts;
}

RNPxlSearch::RNPxlParameterParsing::NucleotideToFragmentAdductMap
RNPxlSearch::RNPxlParameterParsing::getTargetNucleotideToFragmentAdducts(StringList fragment_adducts)
{
  NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts;
  for (String t : fragment_adducts)
  {
    t.removeWhitespaces();

    EmpiricalFormula formula;
    String name;

    // format is: target_nucletide:formula;name
    char target_nucleotide = t[0];
    if (t[1] != ':')
    {
      LOG_WARN << "Missing ':'. Wrong format of fragment_adduct string: " << t << endl;
      return NucleotideToFragmentAdductMap();
    }

    // remove target nucleotide and : character from t
    t = t.substr(2);

    // split into formula and name
    vector<String> fs;
    t.split(";", fs);
    if (fs.size() == 1) // no name provided so we just take the formula as name
    {
      formula = EmpiricalFormula(fs[0]);
      name = formula.toString();
    }
    else if (fs.size() == 2)
    {
      formula = EmpiricalFormula(fs[0]);
      name = fs[1];
    }
    else
    {
      LOG_WARN << "Wrong format of fragment_adduct string: " << t << endl;
      return NucleotideToFragmentAdductMap();
    }

    FragmentAdductDefinition_ fad;
    fad.name = name;
    fad.formula = formula;
    fad.mass = formula.getMonoWeight();

    nucleotide_to_fragment_adducts[target_nucleotide].insert(fad);

    // register all fragment adducts as N- and C-terminal modification (if not already registered)
    if (!ModificationsDB::getInstance()->has(name))
    {
      ResidueModification* c_term = new ResidueModification();
      c_term->setId(name);
      c_term->setName(name);
      c_term->setFullId(name + " (C-term)");
      c_term->setTermSpecificity(ResidueModification::C_TERM);
      c_term->setDiffMonoMass(fad.mass);
      ModificationsDB::getInstance()->addModification(c_term);

      ResidueModification* n_term = new ResidueModification();
      n_term->setId(name);
      n_term->setName(name);
      n_term->setFullId(name + " (N-term)");
      n_term->setTermSpecificity(ResidueModification::N_TERM);
      n_term->setDiffMonoMass(fad.mass);
      ModificationsDB::getInstance()->addModification(n_term);
    }
  }

#ifdef DEBUG_RNPXLSEARCH
  for (auto const & p2fas : precursor_to_fragment_adducts)
    {
      for (auto const & p2fa : p2fas.second)
      {
        cout << "nucleotide:" << p2fas.first
             << " fragment adduct:" << p2fa.formula.toString()
             << " fragment adduct mass:" << p2fa.mass
             << " name:" <<  p2fa.name << endl;
      }
    }
#endif

  return nucleotide_to_fragment_adducts;
}

RNPxlSearch::MS2AdductsOfSinglePrecursorAdduct
RNPxlSearch::RNPxlParameterParsing::getFeasibleFragmentAdducts(const String &exp_pc_adduct,
                                                               const String &exp_pc_formula,
                                                               const RNPxlSearch::RNPxlParameterParsing::NucleotideToFragmentAdductMap &nucleotide_to_fragment_adducts,
                                                               const set<char> &can_xl)
{
  LOG_DEBUG << "Generating fragment adducts for precursor adduct: '" << exp_pc_adduct << "'" << endl;

  MS2AdductsOfSinglePrecursorAdduct ret;

  // no precursor adduct? 
  if (exp_pc_formula.empty()) { return ret; } // no fragment adducts or marker ions are expected!

  // count nucleotides in precursor adduct (e.g.: "TCA-H2O" yields map: T->1, C->1, A->1)
  // and determine the set of cross-linkable nucleotides in the precursor adduct
  size_t nt_count(0);
  map<char, Size> exp_pc_nucleotide_count;
  set<char> exp_pc_xl_nts;
  String::const_iterator exp_pc_it = exp_pc_adduct.begin();
  for (; exp_pc_it != exp_pc_adduct.end(); ++exp_pc_it, ++nt_count)
  {
    // we are finished with nucleotides in string if first loss/gain is encountered
    if (*exp_pc_it == '+' || *exp_pc_it == '-') break;

    // count occurence of nucleotide
    if (exp_pc_nucleotide_count.count(*exp_pc_it) == 0)
    {
      exp_pc_nucleotide_count[*exp_pc_it] = 1;
      if (can_xl.count(*exp_pc_it)) { exp_pc_xl_nts.insert(*exp_pc_it); };
    }
    else
    {
      exp_pc_nucleotide_count[*exp_pc_it]++;
    }
  }

  // check if at least one nucleotide present that can cross link
  bool has_xl_nt(false);
  for (auto const & m : exp_pc_nucleotide_count) { if (can_xl.count(m.first)) { has_xl_nt = true; break; } }

  LOG_DEBUG << "\t" << exp_pc_adduct << " has cross-linkable nucleotide (0 = false, 1 = true): " << has_xl_nt << endl;

  // no cross-linkable nt contained in the precursor adduct? Return an empty fragment adduct definition set
  if (!has_xl_nt) { return ret; }

  // HERE: at least one cross-linkable nt present in precursor adduct

  // extract loss string from precursor adduct (e.g.: "-H2O")
  // String exp_pc_loss_string(exp_pc_it, exp_pc_adduct.end());

  LOG_DEBUG << "\t" << exp_pc_adduct << " is monomer (1 = true, >1 = false): " << nt_count << endl;

  // Handle the cases of monomer or oligo nucleotide bound to the precursor.
  // This distinction is made because potential losses on the precursor only allows us to reduce the set of chemical feasible fragment adducts if they are on a monomer.
  // In the case of an oligo we can't be sure if the cross-linked amino acid or any other in the oligo had the loss.
  if (nt_count > 1)  // No monomer? For every nucleotide that can be cross-linked: Create all fragment adducts without restriction by losses (no restriction as the loss could be on the other nts)
  {
    // for each nucleotide and potential set of fragment adducts
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // the nucleotide without any associated loss
      const set<FragmentAdductDefinition_>& fragment_adducts = n2fa.second; // all potential fragment adducts that may arise from the unmodified nucleotide

      // check if nucleotide is cross-linkable and part of the precursor adduct
      if (exp_pc_xl_nts.find(nucleotide) != exp_pc_xl_nts.end())
      {
        LOG_DEBUG << "\t" << exp_pc_adduct << " found nucleotide: " << String(nucleotide) << " in precursor RNA." << endl;
        LOG_DEBUG << "\t" << exp_pc_adduct << " nucleotide: " << String(nucleotide) << " has fragment_adducts: " << fragment_adducts.size() << endl;

        // store feasible adducts associated with a cross-link with character nucleotide
        vector<FragmentAdductDefinition_> faa;
        std::copy(fragment_adducts.begin(), fragment_adducts.end(), back_inserter(faa));
        ret.feasible_adducts.emplace_back(make_pair(nucleotide, faa));
      }
    }

    // Create set of marker ions for all nucleotides contained in the precursor (including those that do not cross-link.)
    // Note: The non-cross-linked nt in the precursor adduct are more likely to produce the marker ions (=more fragile than the cross-linked nt).
    set<FragmentAdductDefinition_> marker_ion_set;
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // the nucleotide without any associated loss
      if (exp_pc_nucleotide_count.find(nucleotide) != exp_pc_nucleotide_count.end()) { marker_ion_set.insert(n2fa.second.begin(), n2fa.second.end()); }
    }
    std::move(std::begin(marker_ion_set), std::end(marker_ion_set), std::back_inserter(ret.marker_ions));
  }
  else // nt_count == 1: monomer. We need to check if the neutral loss reduces the set of feasible (e.g., chemically sound) fragment adducts
  {
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // one letter code of the nt
      set<FragmentAdductDefinition_> fas = n2fa.second; // all potential fragment adducts that may arise from nt (if no losses are considered)

      // check if nucleotide is cross-linkable and part of the precursor adduct
      if (exp_pc_xl_nts.find(nucleotide) != exp_pc_xl_nts.end())
      {
        LOG_DEBUG << "\t" << exp_pc_adduct << " found nucleotide: " << String(nucleotide) << " in precursor RNA." << endl;
        LOG_DEBUG << "\t" << exp_pc_adduct << " nucleotide: " << String(nucleotide) << " has fragment_adducts: " << fas.size() << endl;

        // check chemical feasibility by checking if subtraction of adduct would result in negative elemental composition
        for (auto it = fas.begin(); it != fas.end(); )
        {
          bool negative_elements = (EmpiricalFormula(exp_pc_formula) - it->formula).toString().hasSubstring("-");

          if (negative_elements) // fragment adduct can't be subformula of precursor adduct
          {
            it = fas.erase(it);
          }
          else
          {
            ++it; // STL erase idiom (mind the pre-increment)
          }
        }

        // store feasible adducts associated with a cross-link with character nucleotide[0]
        vector<FragmentAdductDefinition_> faa;
        std::copy(fas.begin(), fas.end(), back_inserter(faa));
        ret.feasible_adducts.emplace_back(make_pair(nucleotide, faa));

        // store feasible marker ions (in this case (monomer) we also restrict them on the nt on the precursor)
        // Note that these peaks are likely missing or of very low intensity,
        std::copy(std::begin(fas), std::end(fas), std::back_inserter(ret.marker_ions));
      }
    }
  }

  // Because, e.g., ribose might be a feasible fragment of any nucleotide, we keep only one version
  // Note: sort by formula and (as tie breaker) the name
  std::sort(ret.marker_ions.begin(), ret.marker_ions.end(),
    [](FragmentAdductDefinition_ const & a, FragmentAdductDefinition_ const & b)
    {
      const String as = a.formula.toString();
      const String bs = b.formula.toString();
      return std::tie(as, a.name) < std::tie(bs, b.name);
    }
  );
  // Note: for uniqueness, we only rely on the formula (in case of tie: keeping the first = shortest name)
  auto it = std::unique(ret.marker_ions.begin(), ret.marker_ions.end(),
    [](FragmentAdductDefinition_ const & a, FragmentAdductDefinition_ const & b)
    {
      return a.formula == b.formula;
    }
  );
  ret.marker_ions.resize(std::distance(ret.marker_ions.begin(), it));

  // print feasible fragment adducts
  for (auto const & ffa : ret.feasible_adducts)
  {
    const char & nucleotide = ffa.first;
    LOG_DEBUG << "  Cross-linkable nucleotide '" << nucleotide << "' and feasible fragment adducts:" << endl;
    for (auto const & a : ffa.second)
    {
      LOG_DEBUG << "\t" << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
    }
  }

  // print marker ions
  LOG_DEBUG << "  Marker ions:" << endl;
  for (auto const & a : ret.marker_ions)
  {
    LOG_DEBUG << "\t" << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
  }

  return ret;
}

void RNPxlSearch::RNPxlFragmentIonGenerator::addMS2MarkerIons(
  const vector<RNPxlSearch::FragmentAdductDefinition_> &marker_ions, 
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
  for (auto const & m : marker_ions)
  {
    const double mz = m.mass + Constants::PROTON_MASS_U;

    spectrum.emplace_back(mz, 1.0);
    spectrum_charge.emplace_back(1);
    spectrum_annotation.emplace_back(ANNOTATIONS_MARKER_ION_PREFIX + m.name);  // add name (e.g., MI:U-H2O)
  }
}

void RNPxlSearch::RNPxlFragmentIonGenerator::addSpecialLysImmonumIons(
  const String& unmodified_sequence,
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
   if (unmodified_sequence.has('K'))
   {
      const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight(); 
      spectrum.emplace_back(immonium_ion2_mz, 1.0);
      spectrum_charge.emplace_back(1);
      spectrum_annotation.emplace_back(String("iK(C5H10N1)"));

      // usually only observed without shift (A. Stuetzer)
      const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight(); 
      spectrum.emplace_back(immonium_ion3_mz, 1.0);
      spectrum_charge.emplace_back(1);
      spectrum_annotation.emplace_back(String("iK(C6H13N2O)"));
    }
}

void RNPxlSearch::RNPxlFragmentIonGenerator::addShiftedImmoniumIons(const String &unmodified_sequence,
                                                                    const String &fragment_shift_name,
                                                                    const double fragment_shift_mass,
                                                                    PeakSpectrum &partial_loss_spectrum,
                                                                    PeakSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
                                                                    PeakSpectrum::StringDataArray &partial_loss_spectrum_annotation) 
{

  if (unmodified_sequence.hasSubstring("Y"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10NO").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('Y', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("W"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C10H11N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('W', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("F"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('F', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("H"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H8N3").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('H', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("C"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C2H6NS").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('C', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("P"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C4H8N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('P', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("L") || unmodified_sequence.hasSubstring("I"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H12N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('L', fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("K"))
  {
    // classical immonium ion
    const double immonium_ion_mz = EmpiricalFormula("C5H13N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('K', fragment_shift_name));

    // TODO: check if only DNA specific and if also other shifts are observed
    // according to A. Stuetzer mainly observed with C‘-NH3 (94.0167 Da)
    const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion2_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C5H10N1)" + fragment_shift_name));

    // usually only observed without shift (A. Stuetzer)
    const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion3_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C6H13N2O)" + fragment_shift_name));
  }
  else if (unmodified_sequence.hasSubstring("M"))
  {
    const double immonium_ion_mz = 104.05285 + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(RNPxlFragmentAnnotationHelper::getAnnotatedImmoniumIon('M', fragment_shift_name));
  }
}

  /* 
  * Add peaks with shifts induced by the RNA/DNA:
  *   - Precursor with complete NA-oligo for charge 1..z
  *   - Partial shifts (without complete precursor adduct)
  *     - Add shifted immonium ions for charge 1 only
  *     and create shifted shifted b,y,a ions + precursors for charge 1..z (adding the unshifted version and performing the shift)
  *     based on the total_loss_spectrum provide to the method
  */
void RNPxlSearch::RNPxlFragmentIonGenerator::generatePartialLossSpectrum(const String &unmodified_sequence,
                                                                         const double &fixed_and_variable_modified_peptide_weight,
                                                                         const String &precursor_rna_adduct,
                                                                         const double &precursor_rna_weight,
                                                                         const int &precursor_charge,
                                                                         const vector<RNPxlSearch::FragmentAdductDefinition_> &partial_loss_modification,
                                                                         const PeakSpectrum& partial_loss_template_z1,
                                                                         const PeakSpectrum& partial_loss_template_z2,
                                                                         const PeakSpectrum& partial_loss_template_z3,
                                                                         PeakSpectrum &partial_loss_spectrum)
{
  partial_loss_spectrum.getIntegerDataArrays().resize(1);
  PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge = partial_loss_spectrum.getIntegerDataArrays()[0];

  partial_loss_spectrum.getStringDataArrays().resize(1); // annotation
  PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation = partial_loss_spectrum.getStringDataArrays()[0];

  // ADD: (mainly for ETD) MS2 precursor peaks of the MS1 adduct (total RNA) carrying peptide for all z <= precursor charge
  for (int charge = 1; charge <= static_cast<int>(precursor_charge); ++charge)
  {
    addPrecursorWithCompleteRNA_(fixed_and_variable_modified_peptide_weight,
                                 precursor_rna_adduct,
                                 precursor_rna_weight,
                                 charge,
                                 partial_loss_spectrum,
                                 partial_loss_spectrum_charge,
                                 partial_loss_spectrum_annotation);
  }

  // for all observable MS2 adducts ...
  for (Size i = 0; i != partial_loss_modification.size(); ++i)
  {
    // get name and mass of fragment adduct
    const String& fragment_shift_name = partial_loss_modification[i].name; // e.g. U-H2O
    const double fragment_shift_mass = partial_loss_modification[i].mass;

    // ADD: shifted immonium ion peaks of charge 1 (if the amino acid is present in the sequence)
    RNPxlFragmentIonGenerator::addShiftedImmoniumIons(
      unmodified_sequence,
      fragment_shift_name,
      fragment_shift_mass,
      partial_loss_spectrum,
      partial_loss_spectrum_charge,
      partial_loss_spectrum_annotation);

    // annotate generated a,b,y ions with fragment shift name
    PeakSpectrum shifted_series_peaks;
    shifted_series_peaks.getStringDataArrays().resize(1); // annotation
    shifted_series_peaks.getIntegerDataArrays().resize(1); // charge

    PeakSpectrum::StringDataArray& shifted_series_annotations = shifted_series_peaks.getStringDataArrays()[0];
    PeakSpectrum::IntegerDataArray& shifted_series_charges = shifted_series_peaks.getIntegerDataArrays()[0];

    // For every charge state
    for (int z = 1; z <= precursor_charge; ++z)
    {
      // 1. add shifted peaks 
      if (z == 1)
      {
        for (Size i = 0; i != partial_loss_template_z1.size(); ++i) 
        { 
          Peak1D p = partial_loss_template_z1[i];
          p.setMZ(p.getMZ() + fragment_shift_mass / static_cast<double>(z));         
          shifted_series_peaks.push_back(p);
          shifted_series_annotations.push_back(partial_loss_template_z1.getStringDataArrays()[0][i]);
          shifted_series_charges.push_back(partial_loss_template_z1.getIntegerDataArrays()[0][i]);
        } 
      }
      else if (z == 2)
      {
        for (Size i = 0; i != partial_loss_template_z2.size(); ++i) 
        { 
          Peak1D p = partial_loss_template_z2[i];
          p.setMZ(p.getMZ() + fragment_shift_mass / static_cast<double>(z));         
          shifted_series_peaks.push_back(p);
          shifted_series_annotations.push_back(partial_loss_template_z2.getStringDataArrays()[0][i]);
          shifted_series_charges.push_back(partial_loss_template_z2.getIntegerDataArrays()[0][i]);
        } 
      }
      else if (z == 3)
      {
        for (Size i = 0; i != partial_loss_template_z3.size(); ++i) 
        { 
          Peak1D p = partial_loss_template_z3[i];
          p.setMZ(p.getMZ() + fragment_shift_mass / static_cast<double>(z));         
          shifted_series_peaks.push_back(p);
          shifted_series_annotations.push_back(partial_loss_template_z3.getStringDataArrays()[0][i]);
          shifted_series_charges.push_back(partial_loss_template_z3.getIntegerDataArrays()[0][i]);
        } 
      }
      else // don't consider fragment ions with charge >= 4 

      { 
        break; 
      }    
    }

    // 2. add fragment shift name to annotation of shifted peaks
    for (Size j = 0; j != shifted_series_annotations.size(); ++j)
    {
      shifted_series_annotations[j] += " " + fragment_shift_name;
    }

    // append shifted and annotated ion series to partial loss spectrum
    partial_loss_spectrum.insert(partial_loss_spectrum.end(), shifted_series_peaks.begin(), shifted_series_peaks.end());
    // std::move strings during insert
    partial_loss_spectrum_annotation.insert(
      partial_loss_spectrum_annotation.end(),
      make_move_iterator(shifted_series_annotations.begin()),
      make_move_iterator(shifted_series_annotations.end())
    );
    partial_loss_spectrum.getIntegerDataArrays()[0].insert(
      partial_loss_spectrum_charge.end(),
      shifted_series_charges.begin(),
      shifted_series_charges.end()
    );
  }

  partial_loss_spectrum.sortByPosition();
}

void RNPxlSearch::RNPxlFragmentIonGenerator::addPrecursorWithCompleteRNA_(
  const double fixed_and_variable_modified_peptide_weight, 
  const String &precursor_rna_adduct,
  const double precursor_rna_weight, 
  const int charge, 
  PeakSpectrum &partial_loss_spectrum,
  MSSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
  MSSpectrum::StringDataArray &partial_loss_spectrum_annotation)
{
  const double xl_mz = (fixed_and_variable_modified_peptide_weight + precursor_rna_weight +
                  static_cast<double>(charge) * Constants::PROTON_MASS_U)
                 / static_cast<double>(charge);
  partial_loss_spectrum.push_back(Peak1D(xl_mz, 1.0));
  partial_loss_spectrum_charge.push_back(charge);
  partial_loss_spectrum_annotation.push_back(String("[M+") + precursor_rna_adduct + "]");
}

int main(int argc, const char** argv)
{
  RNPxlSearch tool;
  return tool.main(argc, argv);
}

