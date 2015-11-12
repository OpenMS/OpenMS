// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>
#include <OpenMS/ANALYSIS/RNPXL/PScore.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <QtCore/QProcess>

#include <OpenMS/FILTERING/ID/IDFilter.h>

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
using namespace std;
/*
  TODO:
  // check for common annotation error in unimod
  if ((origin == "C-term" || origin == "N-term") && term_specifity == ResidueModification::ANYWHERE)
        {
          continue;
        }
    move predicate member functions to class
*/



struct PeptideHitSequenceLessComparator
{
  bool operator()(const PeptideHit& a, const PeptideHit& b)
  {
    if (a.getSequence().toString() < b.getSequence().toString()) return true;

    return false;
  }

};

class RNPxlSearch :
  public TOPPBase
{
public:
  RNPxlSearch() :
    TOPPBase("RNPxlSearch", "Annotate RNA to peptide crosslinks in MS/MS spectra.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("database", "<file>", "", "input file ");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_tsv", "<file>", "", "tsv output file");
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance (+/- around precursor m/z)", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, false);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance (+/- around fragment m/z)", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

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
    registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

    StringList all_enzymes;
    EnzymesDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);

    // RNPxl specific
    registerTOPPSubsection_("RNPxl", "RNPxl Options");
    registerIntOption_("RNPxl:length", "", 2, "Oligonucleotide maximum length. 0 = disable search for RNA variants.", false);

    registerStringOption_("RNPxl:sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)", false);

    StringList target_nucleotides;
    target_nucleotides.push_back("A=C10H14N5O7P");
    target_nucleotides.push_back("C=C9H14N3O8P");
    target_nucleotides.push_back("G=C10H14N5O8P");
    target_nucleotides.push_back("U=C9H13N2O9P");

    registerStringList_("RNPxl:target_nucleotides", "", target_nucleotides, "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG", false, false);

    StringList mapping;
    mapping.push_back("A->A");
    mapping.push_back("C->C");
    mapping.push_back("G->G");
    mapping.push_back("U->U");

    registerStringList_("RNPxl:mapping", "", mapping, "format: source->target e.g. A->A, ..., U->U, U->X", false, false);

    StringList restrictions;
    restrictions.push_back("A=0");
    restrictions.push_back("C=0");
    restrictions.push_back("U=1");
    restrictions.push_back("G=0");

    registerStringList_("RNPxl:restrictions", "", restrictions, "format: target nucleotide=min_count: e.g U=1 if at least one U must be in the generated sequence.", false, false);

    StringList modifications;
    modifications.push_back("");
    modifications.push_back("-H2O");
    modifications.push_back("-H2O-HPO3");
    modifications.push_back("-HPO3");
    modifications.push_back("-H2O+HPO3");
    modifications.push_back("+HPO3");

    registerStringList_("RNPxl:modifications", "", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3", false, false);

    registerFlag_("RNPxl:CysteineAdduct", "Use this flag if the +152 adduct is expected.");
    registerFlag_("RNPxl:filter_fractional_mass", "Use this flag to filter non-crosslinks by fractional mass.");
    registerFlag_("RNPxl:localization", "Use this flag to perform crosslink localization by partial loss scoring as post-analysis.");
    registerFlag_("RNPxl:carbon_labeled_fragments", "Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).");
    registerDoubleOption_("RNPxl:filter_small_peptide_mass", "<threshold>", 600.0, "Filter precursor that can only correspond to non-crosslinks by mass.", false, true);
    registerDoubleOption_("RNPxl:marker_ions_tolerance", "<tolerance>", 0.05, "Tolerance used to determine marker ions (Da).", false, true);
  }

  // Slimmer structure to store a string representation
  struct IndexedString
  {
    String::const_iterator begin;
    String::const_iterator end; // one after last character in substring

    bool operator<(const IndexedString& other) const
    {
      if (end - begin < other.end - other.begin) return true;

      if (end - begin > other.end - other.begin) return false;

      // same size
      String::const_iterator b = begin;
      String::const_iterator bo = other.begin;

      for (; b != end; ++b, ++bo)
      {
        if (*b < *bo) return true;

        if (*b > *bo) return false;
      }

      return false;
    }

    inline String getString() const
    {
      return String(begin, end);
    }

  };

  // Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
  struct AnnotatedHit
  {
    IndexedString sequence;
    SignedSize peptide_mod_index; // enumeration index of the non-RNA peptide modification
    Size rna_mod_index; // index of the RNA modification
    double score;
    double best_localization_score;
    String localization_scores;
    String best_localization;  
    String fragment_annotation_string;
    static bool hasBetterScore(const AnnotatedHit& a, const AnnotatedHit& b)
    {
      return a.score > b.score;
    }
  };

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      ResidueModification rm;
      if (modification.hasSubstring(" (N-term)"))
      {
        modification.substitute(" (N-term)", "");
        rm = ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::N_TERM);
      }
      else if (modification.hasSubstring(" (C-term)"))
      {
        modification.substitute(" (C-term)", "");
        rm = ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::C_TERM);
      }
      else
      {
        rm = ModificationsDB::getInstance()->getModification(modification);
      }
      modifications.push_back(rm);

      // attempt to register modified residue in the single thread context (no locking required) and obtain thread safety this way
      ResidueDB::getInstance()->getModifiedResidue(modification);
    }

    return modifications;
  }

  // check for minimum size
  class HasInvalidPeptideLengthPredicate
  {
public:
    explicit HasInvalidPeptideLengthPredicate(Size min_size)
      : min_size_(min_size)
    {
    }

    bool operator()(const AASequence& aas)
    {
      return aas.size() < min_size_;
    }

private:
    Size min_size_;
  };

  // spectrum must not contain 0 intensity peaks and must be sorted by m/z
  template <typename SpectrumType>
  void deisotopeAndSingleChargeMSSpectrum(SpectrumType& in, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = true)
  {
    if (in.empty())
    {
      return;
    }

    SpectrumType old_spectrum = in;

    // determine charge seeds and extend them
    vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
    vector<Int> features(old_spectrum.size(), -1);
    Int feature_number = 0;

    for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
    {
      double current_mz = old_spectrum[current_peak].getPosition()[0];

      for (Int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1) // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton) // test for missing peak
            {
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            }
            else
            {
              // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
              Size n_extensions = extensions.size();
              if (n_extensions != 0)
              {
                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
                {
                  if (i < min_isopeaks)
                  {
                    has_min_isopeaks = false;
                  }
                  break;
                }
              }

              // averagine check passed
              extensions.push_back(p);
            }
          }

          if (has_min_isopeaks)
          {
            //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
            mono_isotopic_peak[current_peak] = q;
            for (Size i = 0; i != extensions.size(); ++i)
            {
              features[extensions[i]] = feature_number;
            }
            feature_number++;
          }
        }
      }
    }

    in.clear(false);
    for (Size i = 0; i != old_spectrum.size(); ++i)
    {
      Int z = mono_isotopic_peak[i];
      if (keep_only_deisotoped)
      {
        if (z == 0)
        {
          continue;
        }

        // if already single charged or no decharging selected keep peak as it is
        if (!make_single_charged)
        {
          in.push_back(old_spectrum[i]);
        }
        else
        {
          Peak1D p = old_spectrum[i];
          p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
          in.push_back(p);
        }
      }
      else
      {
        // keep all unassigned peaks
        if (features[i] < 0)
        {
          in.push_back(old_spectrum[i]);
          continue;
        }

        // convert mono-isotopic peak with charge assigned by deisotoping
        if (z != 0)
        {
          if (!make_single_charged)
          {
            in.push_back(old_spectrum[i]);
          }
          else
          {
            Peak1D p = old_spectrum[i];
            p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            in.push_back(p);
          }
        }
      }
    }

    in.sortByPosition();
  }

  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool single_charge_spectra)
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
      deisotopeAndSingleChargeMSSpectrum(exp[exp_index], 1, 3, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 10, single_charge_spectra);

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
    }
  }

  struct FragmentAnnotationDetail_
  {
    String shift;
    int charge;
    double mz;
    double intensity;

    bool operator<(const FragmentAnnotationDetail_& other) const
    {
      if (charge < other.charge) 
      { 
        return true;
      } 
      else if (charge > other.charge)
      {
        return false;
      }

      if (shift < other.shift)
      { 
        return true;
      } 
      else if (shift > other.shift)
      {
        return false;
      }

      if (mz < other.mz)
      {
        return true;
      }
      else if (mz > other.mz)
      {
        return false;
      }

      if (intensity < other.intensity)
      {
        return true;
      }
      else if (intensity > other.intensity)
      {
        return false;
      }

      return false;
    }

    bool operator==(const FragmentAnnotationDetail_& other) const
    {
      double mz_diff = fabs(mz - other.mz);
      double intensity_diff = fabs(intensity - other.intensity);
      return (shift == other.shift && mz_diff < 1e-6 && intensity_diff < 1e-6); // mz and intensity difference comparison actually not needed but kept for completeness
    }
  };

  RichPeak1D getAnnotatedImmoniumIon_(char c, double immonium_ion_mz, const String& fragment_shift_name)
  {
    RichPeak1D RNA_fragment_peak;
    RNA_fragment_peak.setIntensity(1.0);
    RNA_fragment_peak.setMZ(immonium_ion_mz); 
    RNA_fragment_peak.setMetaValue("IonName", String("i") + c + "+" + fragment_shift_name);
    return RNA_fragment_peak;
  }
 
  String fragmentAnnotationDetailsToString_(const String& ion_type, map<Size, vector<FragmentAnnotationDetail_> > ion_annotation_details)
  {
    String fas;
    for (map<Size, vector<FragmentAnnotationDetail_> >::const_iterator ait = ion_annotation_details.begin(); ait != ion_annotation_details.end(); ++ait)
    {
      for (vector<FragmentAnnotationDetail_>::const_iterator sit = ait->second.begin(); sit != ait->second.end(); ++sit)
      {
        if (ait != ion_annotation_details.begin() || sit != ait->second.begin())
        {
          fas += "|";
        }

        String annotation_text;
        annotation_text = sit->shift.empty() ? "[" + ion_type + String(ait->first) + "]" + String(sit->charge, '+') : "[" + ion_type + String(ait->first) + "+" + sit->shift + "]" + String(sit->charge, '+'); // e.g.: [b3]+ and  [y3+H3PO4]++
        // e.g.: (343.5,99.5,"[b2-H2O]+")
        fas += "(" + String::number(sit->mz, 3) + "," + String::number(100.0 * sit->intensity, 1) + "," + "\"" + annotation_text+ "\")";
      }
    }
    return fas;
  }

  String shiftedImmoniumIonsToString_(const map<String, set<pair<String, double > > >& shifted_immonium_ions)
  {
    String fas;
    for (map<String, set<pair<String, double > > >::const_iterator ait = shifted_immonium_ions.begin(); ait != shifted_immonium_ions.end(); ++ait)
    {
      for (set<pair<String, double> >::const_iterator sit = ait->second.begin(); sit != ait->second.end(); ++sit)
      {
        if (ait != shifted_immonium_ions.begin() || sit != ait->second.begin())
        {
          fas += "|";
        }
        fas += sit->first;
      }
    }
    return fas;
  }

  String shiftedMarkerIonsToString_(const map<String, set<String> >& annotated_marker_ions)
  {
    String fas;
    for (map<String, set<String> >::const_iterator ait = annotated_marker_ions.begin(); ait != annotated_marker_ions.end(); ++ait)
    {
      for (set<String>::const_iterator sit = ait->second.begin(); sit != ait->second.end(); ++sit)
      {
        if (ait != annotated_marker_ions.begin() || sit != ait->second.begin())
        {
          fas += "|";
        }
        fas += *sit;
      }
    }
    return fas;
  }


  void postScoreHits_(const PeakMap& exp, vector<vector<AnnotatedHit> >& annotated_hits, Size top_hits, const RNPxlModificationMassesResult& mm, const vector<ResidueModification>& fixed_modifications, const vector<ResidueModification>& variable_modifications, Size max_variable_mods_per_peptide, TheoreticalSpectrumGenerator spectrum_generator, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool carbon_is_labeled)
  {
    // for pscore calculation only
    vector<vector<Size> > rank_map = PScore::calculateRankMap(exp);

    Param ps = spectrum_generator.getParameters();
    ps.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    spectrum_generator.setParameters(ps);

    SpectrumAlignment spectrum_aligner;
    Param pa = spectrum_aligner.getParameters();
    pa.setValue("tolerance", (double)fragment_mass_tolerance, "Defines the absolute (in Da) or relative (in ppm) tolerance in the alignment");
    if (fragment_mass_tolerance_unit_ppm)
    {
      pa.setValue("is_relative_tolerance", "true");
    } 
    else
    {
      pa.setValue("is_relative_tolerance", "false");
    } 
  
    spectrum_aligner.setParameters(pa);

  // remove all but top n scoring for localization (usually all but the first one)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      const PeakSpectrum& exp_spectrum = exp[scan_index];
      const Size precursor_charge = exp_spectrum.getPrecursors()[0].getCharge();
      if (!annotated_hits[scan_index].empty())
      {
        for (vector<AnnotatedHit>::iterator a_it = annotated_hits[scan_index].begin(); a_it != annotated_hits[scan_index].end(); ++a_it)
        {
          // get unmodified string
          String unmodified_sequence = a_it->sequence.getString();
          AASequence aas = AASequence::fromString(unmodified_sequence);

          // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
          vector<AASequence> all_modified_peptides;
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);

          // reannotate much more memory heavy AASequence object
          AASequence fixed_and_variable_modified_peptide = all_modified_peptides[a_it->peptide_mod_index]; 

          // determine RNA on precursor from index in map
          std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, a_it->rna_mod_index);
          String precursor_rna_adduct = *mod_combinations_it->second.begin();

          // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
          // First get all possible RNA fragment shifts in the MS2 (based on the precursor RNA)
          vector<String> partial_loss_modification_names = RNPxlModificationsGenerator::getRNAFragmentModificationNames(precursor_rna_adduct, aas);

          RichPeakSpectrum partial_loss_spectrum;

          // generate total loss spectrum
          RichPeakSpectrum total_loss_spectrum;
          for (Size z = 1; z <= precursor_charge; ++z)
          {
            spectrum_generator.addPeaks(total_loss_spectrum, aas, Residue::AIon, z);
            spectrum_generator.addPeaks(total_loss_spectrum, aas, Residue::BIon, z);
            spectrum_generator.addPeaks(total_loss_spectrum, aas, Residue::YIon, z);
          }

          total_loss_spectrum.sortByPosition();

          //TODO generate unshifted immonium ions

          // Add peaks for marker ions A', G', C' marker ions (presence of these are determined by precursor RNA)
          RichPeak1D RNA_fragment_peak;
          RNA_fragment_peak.setIntensity(1.0);
          if (precursor_rna_adduct.hasSubstring("A"))
          {
            RNA_fragment_peak.setMZ(136.0623); // C5H6N5
            RNA_fragment_peak.setMetaValue("IonName", "A'");
            partial_loss_spectrum.push_back(RNA_fragment_peak);
          }

          if (precursor_rna_adduct.hasSubstring("G"))
          {
            RNA_fragment_peak.setMZ(152.0572); //C5H6N5O
            RNA_fragment_peak.setMetaValue("IonName", "G'");
            partial_loss_spectrum.push_back(RNA_fragment_peak);
          }

          if (precursor_rna_adduct.hasSubstring("C"))
          {
            RNA_fragment_peak.setMZ(112.0510); // C4H6N3O
            RNA_fragment_peak.setMetaValue("IonName", "C'");
            partial_loss_spectrum.push_back(RNA_fragment_peak);
          }

          for (Size i = 0; i != partial_loss_modification_names.size(); ++i)
          {
            const String& fragment_shift_name = partial_loss_modification_names[i]; // e.g. U-H2O

            ResidueModification fragment_modification = ModificationsDB::getInstance()->getTerminalModification(fragment_shift_name, ResidueModification::N_TERM);

            // full labeling of nucleotide if flag is set
            if (carbon_is_labeled)
            {
              const Element* carbon = ElementDB::getInstance()->getElement("Carbon");
              fragment_modification.setDiffMonoMass(fragment_modification.getDiffMonoMass() + fragment_modification.getDiffFormula().getNumberOf(carbon) * Constants::C13C12_MASSDIFF_U); // replace 12C by 13C
            }

            const double fragment_shift_mass = fragment_modification.getDiffMonoMass();

            // RNA mass peak
            RichPeak1D RNA_fragment_peak;
            RNA_fragment_peak.setIntensity(1.0);
            RNA_fragment_peak.setMZ(fragment_shift_mass + Constants::PROTON_MASS_U); // there is exactly one RNA fragment modification that we added to this partial loss spectrum. So get the modification and mass to calculate the RNA peak mass.
            RNA_fragment_peak.setMetaValue("IonName", fragment_shift_name);  // add name (e.g. RNA:U-H2O)
            partial_loss_spectrum.push_back(RNA_fragment_peak);

            // Add shifted immonium ion peaks if the amino acid is present in the sequence
            if (unmodified_sequence.hasSubstring("Y"))
            {
              double immonium_ion_mz = EmpiricalFormula("C8H10NO").getMonoWeight() + fragment_shift_mass; 
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('Y', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("W"))
            {
              double immonium_ion_mz = EmpiricalFormula("C10H11N2").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('W', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("F"))
            {
              double immonium_ion_mz = EmpiricalFormula("C8H10N").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('F', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("H"))
            {
              double immonium_ion_mz = EmpiricalFormula("C5H8N3").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('H', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("C"))
            {
              double immonium_ion_mz = EmpiricalFormula("C2H6NS").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('C', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("P"))
            {
              double immonium_ion_mz = EmpiricalFormula("C4H8N").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('P', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("L") || unmodified_sequence.hasSubstring("I"))
            {
              double immonium_ion_mz = EmpiricalFormula("C5H12N").getMonoWeight() + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('L', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("K"))
            {
              double immonium_ion_mz = 101.10732 + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('K', immonium_ion_mz, fragment_shift_name));
            }
            else if (unmodified_sequence.hasSubstring("M"))
            {
              double immonium_ion_mz = 104.05285 + fragment_shift_mass;
              partial_loss_spectrum.push_back(getAnnotatedImmoniumIon_('M', immonium_ion_mz, fragment_shift_name));
            }

            // generate all possible shifted ion a,b,y ion peaks by putting a modification on the N or C terminus
            AASequence c_term_shifted = aas;
            c_term_shifted.setCTerminalModification(fragment_shift_name);

            AASequence n_term_shifted = aas;
            n_term_shifted.setNTerminalModification(fragment_shift_name);
 
            RichPeakSpectrum shifted_series_peaks;
            for (Size z = 1; z <= precursor_charge; ++z)
            {
              spectrum_generator.addPeaks(shifted_series_peaks, n_term_shifted, Residue::AIon, z); // N-term
              spectrum_generator.addPeaks(shifted_series_peaks, n_term_shifted, Residue::BIon, z); // N-term
              spectrum_generator.addPeaks(shifted_series_peaks, c_term_shifted, Residue::YIon, z); // C-term
            }

            // annotate generated a,b,y ions with fragment shift name
            for (Size j = 0; j != shifted_series_peaks.size(); ++j)
            {
              const String& ion_name = shifted_series_peaks[j].getMetaValue("IonName");
              shifted_series_peaks[j].setMetaValue("IonName", ion_name + " " + fragment_shift_name);
            }
            
            // add shifted and annotated ion series to partial loss spectrum
            partial_loss_spectrum.insert(partial_loss_spectrum.end(), shifted_series_peaks.begin(), shifted_series_peaks.end());            

          }

          partial_loss_spectrum.sortByPosition();

          // fill annotated spectrum information
          set<Size> peak_is_annotated;  // experimental peak index

          // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)
          map<Size, vector<FragmentAnnotationDetail_> > unshifted_b_ions, unshifted_y_ions, unshifted_a_ions, shifted_b_ions, shifted_y_ions, shifted_a_ions;
          map<String, set<pair<String, double> > > shifted_immonium_ions;
          map<String, set<String> > annotated_marker_ions;

          // first annotate total loss peaks (these give no information where the actual shift occured)
          #ifdef DEBUG_RNPXLSEARCH
            cout << "Annotating ion (total loss spectrum): " << aas.toString()  << endl;
          #endif
          std::vector<std::pair<Size, Size> > alignment;
          spectrum_aligner.getSpectrumAlignment(alignment, total_loss_spectrum, exp_spectrum);
          for (vector<std::pair<Size, Size> >::const_iterator pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
          {
            const double fragment_intensity = exp_spectrum[pair_it->second].getIntensity(); // in percent (%)
            const double fragment_mz = exp_spectrum[pair_it->second].getMZ();
            String ion_name = total_loss_spectrum[pair_it->first].getMetaValue("IonName");

            // define which ion names are annotated 
            if (ion_name.hasPrefix("y"))
            { 
                String ion_nr_string = ion_name;
                ion_nr_string.substitute("y", "");
                ion_nr_string.substitute("+", "");
                Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                const AASequence& peptide_sequence = aas.getSuffix(ion_number);
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              peak_is_annotated.insert(pair_it->second);                  

              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              FragmentAnnotationDetail_ d;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              unshifted_y_ions[ion_number].push_back(d);
            }
            else if (ion_name.hasPrefix("b"))
            { 
                String ion_nr_string = ion_name;
                ion_nr_string.substitute("b", "");
                ion_nr_string.substitute("+", "");
                Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                const AASequence& peptide_sequence = aas.getPrefix(ion_number);
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              peak_is_annotated.insert(pair_it->second);                  

              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              FragmentAnnotationDetail_ d;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              unshifted_b_ions[ion_number].push_back(d);
            }
            else if (ion_name.hasPrefix("a"))
            { 
                String ion_nr_string = ion_name;
                ion_nr_string.substitute("a", "");
                ion_nr_string.substitute("+", "");
                Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                const AASequence& peptide_sequence = aas.getPrefix(ion_number);
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
              #endif
              peak_is_annotated.insert(pair_it->second);                  

              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              FragmentAnnotationDetail_ d;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              unshifted_a_ions[ion_number].push_back(d);
            }
          }

          vector<double> sites_sum_score(aas.size(), 0);

          // loss spectrum to the experimental measured one
          alignment.clear();
          spectrum_aligner.getSpectrumAlignment(alignment, partial_loss_spectrum, exp_spectrum);

          if (alignment.empty()) continue;
          
          set<String> observed_immonium_ions;

          for (vector<std::pair<Size, Size> >::const_iterator pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
          {
            // only annotate experimental peak if not annotated as complete loss peak (e.g. b and y ions without shift)
            if (peak_is_annotated.find(pair_it->second) != peak_is_annotated.end())
            {
              continue;
            }

            const double fragment_intensity = exp_spectrum[pair_it->second].getIntensity();
            const double fragment_mz = exp_spectrum[pair_it->second].getMZ();

            String ion_name = partial_loss_spectrum[pair_it->first].getMetaValue("IonName");

            vector<String> f;

            ion_name.split(' ', f);  // e.g. "y3+ RNA:C3O" or just "y2"
            String fragment_shift_name;
            if (f.size() == 2) 
            {
              fragment_shift_name = f[1];
            }

            // define which ion names are annotated 
            if (ion_name.hasPrefix("y"))
            { 
              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("y", "");
              ion_nr_string.substitute("+", "");
              Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
              #endif
              // remove RNA: substring for nicer annotation of ion
              String annotation = fragment_shift_name;
              annotation.substitute("RNA:", "");

              FragmentAnnotationDetail_ d;
              d.shift = annotation;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              shifted_y_ions[ion_number].push_back(d);
            }
            else if (ion_name.hasPrefix("b"))
            { 
              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("b", "");
              ion_nr_string.substitute("+", "");
              Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
              #endif
              // remove RNA: substring for nicer annotation of ion
              String annotation = fragment_shift_name;
              annotation.substitute("RNA:", "");

              FragmentAnnotationDetail_ d;
              d.shift = annotation;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              shifted_b_ions[ion_number].push_back(d);
            }
            else if (ion_name.hasPrefix("a"))
            { 
              Size charge = std::count(ion_name.begin(), ion_name.end(), '+');
              String ion_nr_string = ion_name;
              ion_nr_string.substitute("a", "");
              ion_nr_string.substitute("+", "");
              Size ion_number = (Size)ion_nr_string.toInt();
              #ifdef DEBUG_RNPXLSEARCH
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
              #endif
              // remove RNA: substring for nicer annotation of ion
              String annotation = fragment_shift_name;
              annotation.substitute("RNA:", "");
              FragmentAnnotationDetail_ d;
              d.shift = annotation;
              d.charge = charge;
              d.mz = fragment_mz;
              d.intensity = fragment_intensity;
              shifted_a_ions[ion_number].push_back(d);
            }
            else if (ion_name.hasPrefix("RNA:"))
            {
              #ifdef DEBUG_RNPXLSEARCH
                cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " intensity: " << fragment_intensity << endl;                
              #endif
              // remove RNA prefix from string and annotate ion
              String annotation = ion_name;
              annotation.substitute("RNA:", "");
              annotated_marker_ions[annotation].insert("(" + String::number(fragment_mz, 3) + "," + String::number(100.0 * fragment_intensity, 1) + ",\"" + annotation + "\")");
            }
            else if (ion_name.hasPrefix("i"))
            {
              String ion_nr_string = ion_name;
              String origin = ion_name[1];  // type of immonium ion
              observed_immonium_ions.insert(origin);
            #ifdef DEBUG_RNPXLSEARCH
              cout << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " intensity: " << fragment_intensity << endl;                
            #endif
              // remove RNA: substring for nicer annotation of ion
              String annotation = ion_name;
              annotation.substitute("RNA:", "");
              shifted_immonium_ions[origin].insert(make_pair("(" + String::number(fragment_mz, 3) + "," + String::number(100.0 * fragment_intensity, 1) + ",\"" + annotation + "\")", fragment_intensity));
            }
          }

          for (Size i = 0; i != sites_sum_score.size(); ++i) 
          {            
            sites_sum_score[i] = 0.0;
            // caclculate how many ions are explained by this position
            //
            // b-ions
            for (Size bi = 1; bi <= sites_sum_score.size(); ++bi) 
            {      
              double distance = fabs((static_cast<double>(bi - 1) - static_cast<double>(i)) / static_cast<double>(sites_sum_score.size() - 1.0));
              if ((bi - 1) < i) // penalize contradicting observations
              {
                if (shifted_b_ions.find(bi) != shifted_b_ions.end())
                {
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_b_ions[bi].begin(); k != shifted_b_ions[bi].end(); ++k)
                  {
                    sites_sum_score[i] -= 2.0 * (1.0 - distance) * k->intensity;
                  }
                }
              } 
              else // add supporting observations
              {                            
                if (shifted_b_ions.find(bi) != shifted_b_ions.end())
                {  
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_b_ions[bi].begin(); k != shifted_b_ions[bi].end(); ++k)
                  {
                    sites_sum_score[i] += (1.0 - distance) * k->intensity; // weight by distance
                  }
                }
              } 
            }
            
            // a-ions
            for (Size ai = 1; ai <= sites_sum_score.size(); ++ai) 
            {            
              double distance = fabs((static_cast<double>(ai - 1) - static_cast<double>(i)) / static_cast<double>(sites_sum_score.size() - 1.0));
              if ((ai - 1) < i) // penalize contradicting observations
              {
                if (shifted_a_ions.find(ai) != shifted_a_ions.end())
                {
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_a_ions[ai].begin(); k != shifted_a_ions[ai].end(); ++k)
                  {
                    sites_sum_score[i] -= 2.0 * (1.0 - distance) * k->intensity;
                  }
                }
              } 
              else // add supporting observations
              {                            
                if (shifted_a_ions.find(ai) != shifted_a_ions.end())
                {  
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_a_ions[ai].begin(); k != shifted_a_ions[ai].end(); ++k)
                  {
                    sites_sum_score[i] += (1.0 - distance) * k->intensity;
                  }
                }
              } 
            }
            
            // y-ions
            for (Size yi = 1; yi <= sites_sum_score.size(); ++yi) 
            {
              Size position = sites_sum_score.size() - yi;
              double distance = fabs((static_cast<double>(i) - static_cast<double>(position)) / static_cast<double>(sites_sum_score.size() - 1.0));
              if (position > i) // penalize contradicting observations
              {
                if (shifted_y_ions.find(yi) != shifted_y_ions.end())
                {
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_y_ions[yi].begin(); k != shifted_y_ions[yi].end(); ++k)
                  {
                    sites_sum_score[i] -= 2.0 * (1.0 - distance) * k->intensity;
                  }
                }
              }
              else
              {
                if (shifted_y_ions.find(yi) != shifted_y_ions.end())
                {
                  for (vector<FragmentAnnotationDetail_>::const_iterator k = shifted_y_ions[yi].begin(); k != shifted_y_ions[yi].end(); ++k)
                  {
                    sites_sum_score[i] += (1.0 - distance) * k->intensity;
                  }
                }
              }              
            }
          }

          #ifdef DEBUG_RNPXLSEARCH
            cout << "Localisation based on immonium ions: ";
          #endif
          String aas_unmodified = aas.toUnmodifiedString();
          for (Size i = 0; i != aas_unmodified.size(); ++i)
          {
            String origin = String(aas_unmodified[i]);
            if (shifted_immonium_ions.find(origin) != shifted_immonium_ions.end())
            {                                
              #ifdef DEBUG_RNPXLSEARCH
                cout << i+1 << " ";
              #endif
              for (set<pair<String, double> >::const_iterator k = shifted_immonium_ions[origin].begin(); k != shifted_immonium_ions[origin].end(); ++k)
              {
                sites_sum_score[i] += k->second; 
              }
            }
          }

          String best_localization = unmodified_sequence;
          double best_localization_score = 0;
          String localization_scores;
          for (Size i = 0; i != sites_sum_score.size(); ++i)
          {
            if (sites_sum_score[i] > best_localization_score) best_localization_score = sites_sum_score[i];
          }

          for (Size i = 0; i != sites_sum_score.size(); ++i)
          {
            if (i != 0) localization_scores += ' ';
            if (sites_sum_score[i] > 0 )
            {
              localization_scores += String::number(100.0 * sites_sum_score[i], 2);
            }
            else
            {
              localization_scores += "0";
            }

            if (best_localization_score > 0.0 && sites_sum_score[i] >= best_localization_score - 1e-6) best_localization[i] = tolower(best_localization[i]);
          }

          // store score of best localization(s)
          a_it->localization_scores = localization_scores;
          a_it->best_localization = best_localization;
          a_it->best_localization_score = best_localization_score;

          StringList fragment_annotations;
          String ub = fragmentAnnotationDetailsToString_("b", unshifted_b_ions);
          if (!ub.empty()) fragment_annotations.push_back(ub);
          String uy = fragmentAnnotationDetailsToString_("y", unshifted_y_ions);
          if (!uy.empty()) fragment_annotations.push_back(uy);
          String ua = fragmentAnnotationDetailsToString_("a", unshifted_a_ions);
          if (!ua.empty()) fragment_annotations.push_back(ua);
          String sb = fragmentAnnotationDetailsToString_("b", shifted_b_ions);
          if (!sb.empty()) fragment_annotations.push_back(sb);
          String sy = fragmentAnnotationDetailsToString_("y", shifted_y_ions);
          if (!sy.empty()) fragment_annotations.push_back(sy);
          String sa = fragmentAnnotationDetailsToString_("a", shifted_a_ions);
          if (!sa.empty()) fragment_annotations.push_back(sa);
          String sii = shiftedImmoniumIonsToString_(shifted_immonium_ions);
          if (!sii.empty()) fragment_annotations.push_back(sii);
          String smi = shiftedMarkerIonsToString_(annotated_marker_ions);
          if (!smi.empty()) fragment_annotations.push_back(smi);

          a_it->fragment_annotation_string = ListUtils::concatenate(fragment_annotations, "|");

          #ifdef DEBUG_RNPXLSEARCH
            cout << "Ion centric annotation: " << endl;
            cout << "unshifted b ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("b", unshifted_b_ions) << endl;
            cout << "unshifted y ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("y", unshifted_y_ions) << endl;
            cout << "unshifted a ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("a", unshifted_a_ions) << endl;
            cout << "shifted b ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("b", shifted_b_ions) << endl;
            cout << "shifted y ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("y", shifted_y_ions) << endl;
            cout << "shifted a ions: " << endl;
            cout << fragmentAnnotationDetailsToString_("a", shifted_a_ions) << endl;
            cout << "shifted immonium ions: " << endl;
            cout << shiftedImmoniumIonsToString_(shifted_immonium_ions) << endl;
            cout << "shifted marker ions: " << endl;
            cout << shiftedMarkerIonsToString_(annotated_marker_ions) << endl;
            cout << "Localization scores: ";
            cout << localization_scores << endl;
            cout << "Localisation based on ion series and immonium ions of all observed fragments: ";
            cout << best_localization << endl;
          #endif
        }
      }
    }
  }

  void postProcessHits_(const PeakMap& exp, vector<vector<AnnotatedHit> >& annotated_hits, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids, Size top_hits, const RNPxlModificationMassesResult& mm, const vector<ResidueModification>& fixed_modifications, const vector<ResidueModification>& variable_modifications, Size max_variable_mods_per_peptide)
  {
  // remove all but top n scoring
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
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
        Size charge = exp[scan_index].getPrecursors()[0].getCharge();

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        for (vector<AnnotatedHit>::const_iterator a_it = annotated_hits[scan_index].begin(); a_it != annotated_hits[scan_index].end(); ++a_it)
        {
          PeptideHit ph;
          ph.setCharge(charge);

          // get unmodified string
          AASequence aas = AASequence::fromString(a_it->sequence.getString());

          // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
          vector<AASequence> all_modified_peptides;
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);

          // reannotate much more memory heavy AASequence object
          AASequence fixed_and_variable_modified_peptide = all_modified_peptides[a_it->peptide_mod_index]; 
          ph.setScore(a_it->score);

          // determine RNA modification from index in map
          std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
          std::advance(mod_combinations_it, a_it->rna_mod_index);
          ph.setMetaValue(String("RNPxl:RNA"), *mod_combinations_it->second.begin()); // return first nucleotide formula matching the index of the empirical formula
          ph.setMetaValue(String("RNPxl:RNA_MASS_z0"), EmpiricalFormula(mod_combinations_it->first).getMonoWeight()); // RNA uncharged mass via empirical formula

          ph.setMetaValue(String("RNPxl:best_localization_score"), a_it->best_localization_score);
          ph.setMetaValue(String("RNPxl:localization_scores"), a_it->localization_scores);
          ph.setMetaValue(String("RNPxl:best_localization"), a_it->best_localization);

          if (!a_it->fragment_annotation_string.empty())
          {
            ph.setMetaValue(String("RNPxl:fragment_annotation"), a_it->fragment_annotation_string);
          }
          // set the amino acid sequence (for complete loss spectra this is just the variable and modified peptide. For partial loss spectra it additionally contains the loss induced modification)
          ph.setSequence(fixed_and_variable_modified_peptide);
          phs.push_back(ph);
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
    search_parameters.missed_cleavages = getIntOption_("peptide:missed_cleavages");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    search_parameters.precursor_tolerance = getDoubleOption_("precursor:mass_tolerance");
    protein_ids[0].setSearchParameters(search_parameters);
  }


  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in_mzml = getStringOption_("in");
    String in_db = getStringOption_("database");
    String out_idxml = getStringOption_("out");
    String out_csv = getStringOption_("out_tsv");

    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");

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

    vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
    vector<ResidueModification> variable_modifications = getModifications_(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    Int report_top_hits = getIntOption_("report:top_hits");

    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides = getStringList_("RNPxl:target_nucleotides");

    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings = getStringList_("RNPxl:mapping");

    // string format: target,min_count: e.g "X=1" if at least one tU must be in the generated sequence.
    // All target nucleotides must be included. X=0 -> disable restriction
    StringList restrictions = getStringList_("RNPxl:restrictions");

    StringList modifications = getStringList_("RNPxl:modifications");

    String sequence_restriction = getStringOption_("RNPxl:sequence");

    Int max_nucleotide_length = getIntOption_("RNPxl:length");

    bool cysteine_adduct = getFlag_("RNPxl:CysteineAdduct");

    bool localization = getFlag_("RNPxl:localization");

    bool carbon_is_labeled = getFlag_("RNPxl:carbon_labeled_fragments");

    RNPxlModificationMassesResult mm;

    if (max_nucleotide_length != 0)
    {
      mm = RNPxlModificationsGenerator::initModificationMassesRNA(target_nucleotides, mappings, restrictions, modifications, sequence_restriction, cysteine_adduct, max_nucleotide_length);
    }

    mm.mod_masses[""] = 0; // insert "null" modification otherwise peptides without RNA will not be searched
    mm.mod_combinations[""].insert("none");

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

    progresslogger.startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, true);
    progresslogger.endProgress();

    // build multimap of precursor mass to scan index
    Size fractional_mass_filtered(0);
    Size small_peptide_mass_filtered(0);
    multimap<double, Size> multimap_mass_2_scan_index;
    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();

        if (precursor_charge < min_precursor_charge || precursor_charge > max_precursor_charge)
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();
        double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

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

        multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
      }
    }

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;

    vector<vector<AnnotatedHit> > annotated_hits(spectra.size(), vector<AnnotatedHit>());

    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_db, fasta_db);
    progresslogger.endProgress();

    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    EnzymaticDigestion digestor;
    digestor.setEnzyme(getStringOption_("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

    progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    set<IndexedString> processed_petides;

    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    Size count_proteins = 0;
    Size count_peptides = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++count_proteins;

      IF_MASTERTHREAD
      {
        progresslogger.setProgress((SignedSize)fasta_index * NUMBER_OF_THREADS);
      }

      vector<pair<String::const_iterator, String::const_iterator> > current_digest;
      digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<pair<String::const_iterator, String::const_iterator> >::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
        IndexedString string_idx;
        string_idx.begin = cit->first;
        string_idx.end = cit->second;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          if (processed_petides.find(string_idx) != processed_petides.end())
          {
            // peptide (and all modified variants) already processed so skip it
            already_processed = true;
          }
        }

        if (already_processed)
        {
          continue;
        }

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          processed_petides.insert(string_idx);
        }

#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;
        vector<AASequence> all_modified_peptides;

        // no critial section is needed despite ResidueDB not beeing thread sage.
        // It is only written to on introduction of novel modified residues. These residues have been already added above (single thread context).
        {
          AASequence aas = AASequence::fromString(String(cit->first, cit->second));
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass_without_RNA = candidate.getMonoWeight();

          //create empty theoretical spectrum
          MSSpectrum<RichPeak1D> complete_loss_spectrum = MSSpectrum<RichPeak1D>();

          // iterate over all RNA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          Size rna_mod_index = 0;
          for (std::map<String, double>::const_iterator rna_mod_it = mm.mod_masses.begin(); rna_mod_it != mm.mod_masses.end(); ++rna_mod_it, ++rna_mod_index)
          {
            double current_peptide_mass = current_peptide_mass_without_RNA + rna_mod_it->second; // add RNA mass

            // determine MS2 precursors that match to the current peptide mass
            multimap<double, Size>::const_iterator low_it;
            multimap<double, Size>::const_iterator up_it;

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

            if (low_it == up_it) continue; // no matching precursor in data

            //add peaks for b and y ions with charge 1
            if (complete_loss_spectrum.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              spectrum_generator.getSpectrum(complete_loss_spectrum, candidate, 1);
              complete_loss_spectrum.sortByPosition(); //sort by mz
            }

            for (; low_it != up_it; ++low_it)
            {
              const Size& scan_index = low_it->second;
              const PeakSpectrum& exp_spectrum = spectra[scan_index];

              double score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, complete_loss_spectrum);

              // no good hit
              if (score < 1.0)
              {
                continue;
              }

              // add peptide hit
              AnnotatedHit ah;
              ah.sequence.begin = cit->first;
              ah.sequence.end = cit->second;
              ah.peptide_mod_index = mod_pep_idx;
              ah.score = score;
              ah.rna_mod_index = rna_mod_index;

              #ifdef DEBUG_RNPXLSEARCH
                cout << "best score in pre-score: " << score << endl;
              #endif

#ifdef _OPENMP
#pragma omp critical (annotated_hits_access)
#endif
              {

                annotated_hits[scan_index].push_back(ah);
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    cout << "Proteins: " << count_proteins << endl;
    cout << "Peptides: " << count_peptides << endl;
    cout << "Processed peptides: " << processed_petides.size() << endl;

    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;
    progresslogger.startProgress(0, 1, "Post-processing PSMs...");

    if (localization)
    {
      // reload spectra from disc
      spectra.clear(true);
      f.load(in_mzml, spectra);
      spectra.sortSpectra(true);    
      preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false); // for post scoring don't convert fragments to single charge as we need this information
      progresslogger.startProgress(0, 1, "localization...");
      postScoreHits_(spectra, annotated_hits, report_top_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, spectrum_generator, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, carbon_is_labeled);
    }

    progresslogger.startProgress(0, 1, "annotation...");
    postProcessHits_(spectra, annotated_hits, protein_ids, peptide_ids, report_top_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide);
    progresslogger.endProgress();

    // annotate RNPxl related information to hits and create report
    vector<RNPxlReportRow> csv_rows = RNPxlReport::annotate(spectra, peptide_ids, marker_ions_tolerance);

    // write ProteinIdentifications and PeptideIdentifications to IdXML
    IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

    QStringList qparam;
    qparam << "-in" << out_idxml.toQString() << "-out" << out_idxml.toQString() << "-fasta" << in_db.toQString() << "-prefix" << "-enzyme:specificity" << "none" << "-missing_decoy_action" << "warn";
    Int status;
#ifdef OPENMS_WINDOWSPLATFORM
    if (db_name_contains_space)
    {
      // for some reason QProcess doesn't handle escaped " in arguments properly so we use a system call
      String call_string = "PeptideIndexer " + ListUtils::concatenate(qparam, " ");
      writeDebug_(call_string, 5);
      status = system(call_string.c_str());
    }
    else
    {
      status = QProcess::execute("PeptideIndexer", qparam);
    }
#else
    status = QProcess::execute("PeptideIndexer", qparam);
#endif

    if (status != 0)
    {
      writeLog_("Error: Calling PeptideIndexer resulted in an error.");
      return EXTERNAL_PROGRAM_ERROR;
    }

    IdXMLFile().load(out_idxml, protein_ids, peptide_ids);
    csv_rows = RNPxlReport::annotate(spectra, peptide_ids, marker_ions_tolerance);

    // save report
    TextFile csv_file;
    csv_file.addLine(RNPxlReportRowHeader().getString("\t"));
    for (Size i = 0; i != csv_rows.size(); ++i)
    {
      csv_file.addLine(csv_rows[i].getString("\t"));
    }
    csv_file.store(out_csv);
    
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  RNPxlSearch tool;
  return tool.main(argc, argv);
}

