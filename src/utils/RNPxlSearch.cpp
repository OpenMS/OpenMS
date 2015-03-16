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

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
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


#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <map>
#include <algorithm>
 
#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

using namespace OpenMS;
using namespace std;

/*
  TODO:
    proper C-term N-term handling of terminal modifications that can be at every amino acid

        // should be something like this: check if AA of modification and peptide match
        if (origin != aa_seq[pos].getOneLetterCode() && origin != "C-term" && origin != "N-term")
        {
          continue;
        }

  // check for common annotation error in unimod
  if ((origin == "C-term" || origin == "N-term") && term_specifity == ResidueModification::ANYWHERE)
        {
          continue;
        }
    move predicate member functions to class
*/

struct MarkerIonExtractor
{
  typedef map<String, vector<pair<double, double> > > MarkerIonsType;
  static MarkerIonsType extractMarkerIons(const PeakSpectrum& s, const double marker_tolerance)
  {
    MarkerIonsType marker_ions;
    marker_ions["A"].push_back(make_pair(136.06231, 0.0));
    marker_ions["A"].push_back(make_pair(330.06033, 0.0));
    marker_ions["C"].push_back(make_pair(112.05108, 0.0));
    marker_ions["C"].push_back(make_pair(306.04910, 0.0));
    marker_ions["G"].push_back(make_pair(152.05723, 0.0));
    marker_ions["G"].push_back(make_pair(346.05525, 0.0));
    marker_ions["U"].push_back(make_pair(113.03509, 0.0));
    marker_ions["U"].push_back(make_pair(307.03311, 0.0));

    PeakSpectrum spec(s);
    Normalizer normalizer;
    normalizer.filterSpectrum(spec);
    spec.sortByPosition();

    // for each nucleotide with marker ions
    for (Map<String, vector<pair<double, double> > >::iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      // for each marker ion of the current nucleotide
      for (Size i = 0; i != it->second.size(); ++i)
      {
        double mz = it->second[i].first;
        double max_intensity = 0;
        for (PeakSpectrum::ConstIterator sit = spec.begin(); sit != spec.end(); ++sit)
        {
          if (sit->getMZ() + marker_tolerance < mz)
          {
            continue;
          }
          if (mz < sit->getMZ() - marker_tolerance)
          {
            break;
          }
          if (fabs(mz - sit->getMZ()) < marker_tolerance)
          {
            if (max_intensity < sit->getIntensity())
            {
              max_intensity = sit->getIntensity();
            }
          }
        }
        it->second[i].second = max_intensity;
      }
    }
    return marker_ions;
  }
};

struct RNPxlReportRow
{
  bool no_id;
  double rt;
  double original_mz;
  String accessions;
  String RNA;
  String peptide;
  Int charge;
  double score;
  double peptide_weight;
  double RNA_weight;
  double xl_weight;
  double abs_prec_error;
  double rel_prec_error;
  MarkerIonExtractor::MarkerIonsType marker_ions;
  double m_H;
  double m_2H;
  double m_3H;
  double m_4H;

  String getString(String separator)
  {
    StringList sl;

    // rt mz
    sl << String::number(rt, 2) << String::number(original_mz, 4);

    // id if available
    if (no_id)
    {
      sl << "" << "" << "" << "" << "" << "" << "" << "";
    }
    else
    {
      sl << accessions << RNA << peptide << String(charge) << String(score)
         << String::number(peptide_weight, 4) << String::number(RNA_weight, 4) << String::number(peptide_weight + RNA_weight, 4);
    }

    // marker ions
    for (Map<String, vector<pair<double, double> > >::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String::number(it->second[i].second * 100.0, 2);
      }
    }

    // id error and multiple charged mass
    if (no_id)
    {
      sl << "" << ""
         << "" << "" << "" << "";
    }
    else
    {
      // error
      sl << String::number(abs_prec_error, 4)
         << String::number(rel_prec_error, 1);

      // weight
      sl << String::number(m_H, 4)
         << String::number(m_2H, 4)
         << String::number(m_3H, 4)
         << String::number(m_4H, 4);
    }
    return ListUtils::concatenate(sl, separator);
  }

};



struct RNPxlReportRowHeader
{
  String getString(String separator)
  {
    StringList sl;
    sl << "#RT" << "original m/z" << "proteins" << "RNA" << "peptide" << "charge" << "score"
       << "peptide weight" << "RNA weight" << "cross-link weight";

    // marker ion fields
    MarkerIonExtractor::MarkerIonsType marker_ions = MarkerIonExtractor::extractMarkerIons(PeakSpectrum(), 0.0); // call only to generate header entries
    for (MarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String(it->first + "_" + it->second[i].first);
      }
    }
    sl << "abs prec. error Da" << "rel. prec. error ppm" << "M+H" << "M+2H" << "M+3H" << "M+4H";
    return ListUtils::concatenate(sl, separator);
  }

};

// create report
vector<RNPxlReportRow> annotateRNPxlInformation_(const PeakMap& spectra, vector<PeptideIdentification>& peptide_ids, double marker_ions_tolerance)
{
  map<Size, Size> map_spectra_to_id;
  for (Size i = 0; i != peptide_ids.size(); ++i)
  {
    OPENMS_PRECONDITION(!peptide_ids[i].getHits().empty(), "Error: no empty peptide ids allowed.");
    Size scan_index = (unsigned int)peptide_ids[i].getMetaValue("scan_index");
    map_spectra_to_id[scan_index] = i;
  }

  vector<RNPxlReportRow> csv_rows;

  for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
  {
    int scan_index = s_it - spectra.begin();
    vector<Precursor> precursor = s_it->getPrecursors();

    // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
    if (s_it->getMSLevel() == 2 && precursor.size() == 1)
    {
      Size charge = precursor[0].getCharge();
      double mz = precursor[0].getMZ();
      MarkerIonExtractor::MarkerIonsType marker_ions = MarkerIonExtractor::extractMarkerIons(*s_it, marker_ions_tolerance);

      double rt = s_it->getRT();

      RNPxlReportRow row;

      // case 1: no peptide identification: store rt, mz, charge and marker ion intensities
      if (map_spectra_to_id.find(scan_index) == map_spectra_to_id.end())
      {	   
        row.no_id = true;
        row.rt = rt;
        row.original_mz = mz;
	row.charge = charge;
        row.marker_ions = marker_ions;
        csv_rows.push_back(row);
        continue;
      }

      PeptideIdentification& pi = peptide_ids[map_spectra_to_id[scan_index]];
      vector<PeptideHit>& phs = pi.getHits();

      // case 2: identification data present for spectrum
      PeptideHit& ph = phs[0];
      const AASequence& sequence = ph.getSequence();
      double peptide_weight = sequence.getMonoWeight();
      String rna_name = ph.getMetaValue("RNPxl:RNA");
      double rna_weight = ph.getMetaValue("RNPxl:RNA_MASS_z0");

      // crosslink weight for different charge states
      double weight_z1 = (peptide_weight + rna_weight + 1.0 * Constants::PROTON_MASS_U);
      double weight_z2 = (peptide_weight + rna_weight + 2.0 * Constants::PROTON_MASS_U) / 2.0;
      double weight_z3 = (peptide_weight + rna_weight + 3.0 * Constants::PROTON_MASS_U) / 3.0;
      double weight_z4 = (peptide_weight + rna_weight + 4.0 * Constants::PROTON_MASS_U) / 4.0;

      double xl_weight = peptide_weight + rna_weight;
      double theo_mz = (xl_weight + static_cast<double>(charge) * Constants::PROTON_MASS_U) / (double)charge;
      double absolute_difference = theo_mz - mz;
      double ppm_difference =  absolute_difference / theo_mz * 1e6;

      String protein_accessions;
      set<String> accs = ph.extractProteinAccessions();

      // concatenate set into String
      for (set<String>::const_iterator a_it = accs.begin(); a_it != accs.end(); ++a_it)
      {
        if (a_it != accs.begin())
        {
          protein_accessions += ",";
        }
        protein_accessions += *a_it;
      }

      row.no_id = false;
      row.rt = rt;
      row.original_mz = mz;
      row.accessions = protein_accessions;
      row.RNA = rna_name;
      row.peptide = sequence.toString();
      row.charge = charge;
      row.score = ph.getScore();
      row.peptide_weight = peptide_weight;
      row.RNA_weight = rna_weight;
      row.xl_weight = peptide_weight + rna_weight;

      ph.setMetaValue("RNPxl:peptide_mass_z0", DataValue(peptide_weight));
      ph.setMetaValue("RNPxl:xl_mass_z0", xl_weight);

      for (MarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
      {
        for (Size i = 0; i != it->second.size(); ++i)
        {
          ph.setMetaValue(it->first + "_" + it->second[i].first, static_cast<double>(it->second[i].second * 100.0));
        }
      }

      row.marker_ions = marker_ions;
      row.abs_prec_error = absolute_difference;
      row.rel_prec_error = ppm_difference;
      row.m_H = weight_z1;
      row.m_2H = weight_z2;
      row.m_3H = weight_z3;
      row.m_4H = weight_z4;

      ph.setMetaValue("RNPxl:Da difference", (double)absolute_difference);
      ph.setMetaValue("RNPxl:ppm difference", (double)ppm_difference);
      ph.setMetaValue("RNPxl:z1 mass", (double)weight_z1);
      ph.setMetaValue("RNPxl:z2 mass", (double)weight_z2);
      ph.setMetaValue("RNPxl:z3 mass", (double)weight_z3);
      ph.setMetaValue("RNPxl:z4 mass", (double)weight_z4);

      csv_rows.push_back(row);

    }
  }

  return csv_rows;
}

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

    registerOutputFile_("out_csv", "<file>", "", "csv output file");
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, false);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance", false);

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
    registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);

    // RNPxl specific
    registerTOPPSubsection_("RNPxl", "RNPxl Options");
    registerIntOption_("RNPxl:length", "", 4, "Oligonucleotide maximum length.", false);

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
    restrictions.push_back("U=0");
    restrictions.push_back("G=0");

    registerStringList_("RNPxl:restrictions", "", restrictions, "format: target nucleotide=min_count: e.g U=1 if at least one U must be in the generated sequence.", false, false);

    StringList modifications;
    modifications.push_back("-H2O");
    modifications.push_back("");
    modifications.push_back("-H2O-HPO3");
    modifications.push_back("-HPO3");
    modifications.push_back("-H2O+HPO3");
    modifications.push_back("+HPO3");

    registerStringList_("RNPxl:modifications", "", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3", false, false);

    registerFlag_("RNPxl:CysteineAdduct", "Use this flag if the +152 adduct is expected.");
    registerFlag_("RNPxl:filter_fractional_mass", "Use this flag to filter non-crosslinks by fractional mass.");
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
      if (end-begin < other.end-other.begin) return true;

      if (end-begin > other.end-other.begin) return false;

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
      ResidueModification rm = ModificationsDB::getInstance()->getModification(modification);
      modifications.push_back(rm);
      // attempt to register modified residue in the single thread context (no locking required) and obtain thread safety this way
      ResidueDB::getInstance()->getModifiedResidue(modification);
    }

    return modifications;
  }

  // check if for minimum size
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

      for (Int q = max_charge; q >= min_charge; --q)     // important: test charge hypothesis from high to low
      {
        // try to extend isotopes from mono-isotopic peak
        // if extension larger then min_isopeaks possible:
        //   - save charge q in mono_isotopic_peak[]
        //   - annotate all isotopic peaks with feature number
        if (features[current_peak] == -1)     // only process peaks which have no assigned feature number
        {
          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton)     // test for missing peak
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

  double logfactorial(UInt x)
  {
    UInt y;

    if (x < 2)
      return 1;
    else
    {
      double z = 0;
      for (y = 2; y <= x; y++)
      {
        z = log((double)y) + z;
      }

      return z;
    }
  }

  double computeHyperScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const MSSpectrum<Peak1D>& exp_spectrum, const MSSpectrum<RichPeak1D>& theo_spectrum)
  {
    double dot_product = 0.0;
    UInt y_ion_count = 0;
    UInt b_ion_count = 0;

    for (MSSpectrum<RichPeak1D>::ConstIterator theo_peak_it = theo_spectrum.begin(); theo_peak_it != theo_spectrum.end(); ++theo_peak_it)
    {
      const double& theo_mz = theo_peak_it->getMZ();

      double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
      Size index = exp_spectrum.findNearest(theo_mz);
      double exp_mz = exp_spectrum[index].getMZ();

      // found peak match
      if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
      {
        dot_product += exp_spectrum[index].getIntensity();
        if (theo_peak_it->getMetaValue("IonName").toString()[0] == 'y')
        {
          ++y_ion_count;
        }
        else
        {
          ++b_ion_count;
        }
      }
    }

    if (dot_product > 1e-1)
    {
      double yFact = logfactorial(y_ion_count);
      double bFact = logfactorial(b_ion_count);
      double hyperScore = log(dot_product) + yFact + bFact;
      return hyperScore;
    }
    else
    {
      return 0;
    }
  }

  void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
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
      deisotopeAndSingleChargeMSSpectrum(exp[exp_index], 1, 3, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 10, true);

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
    }
  }

  void postProcessHits_(const PeakMap& exp, vector<vector<AnnotatedHit> >& annotated_hits, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids, Size top_hits, const RNPxlModificationMassesResult& mm, const vector<ResidueModification>& fixed_modifications, const vector<ResidueModification>& variable_modifications, Size max_variable_mods_per_peptide)
  {
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
	  ph.setSequence(all_modified_peptides[a_it->peptide_mod_index]);
	  ph.setScore(a_it->score);

	  // determine RNA modification from index in map
    std::map<String, std::set<String> >::const_iterator mod_combinations_it = mm.mod_combinations.begin();
	  std::advance(mod_combinations_it, a_it->rna_mod_index);
    ph.setMetaValue(String("RNPxl:RNA"), *mod_combinations_it->second.begin()); // return first nucleotide formula matching the index of the empirical formula
    ph.setMetaValue(String("RNPxl:RNA_MASS_z0"), EmpiricalFormula(mod_combinations_it->first).getMonoWeight()); // RNA uncharged mass via empirical formula
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
  }


  ExitCodes main_(int, const char**)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    String in_mzml = getStringOption_("in");
    String in_db = getStringOption_("database");
    String out_idxml = getStringOption_("out");
    String out_csv = getStringOption_("out_csv");

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

    RNPxlModificationMassesResult mm = RNPxlModificationsGenerator::initModificationMassesRNA(target_nucleotides, mappings, restrictions, modifications, sequence_restriction, cysteine_adduct, max_nucleotide_length);
    mm.mod_masses[""] = 0;    // insert "null" modification otherwise unmodified peptide will not be searched
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
    preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
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
    digestor.setEnzyme(EnzymaticDigestion::ENZYME_TRYPSIN);
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
          MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();

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
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * current_peptide_mass * precursor_mass_tolerance * 1e-6);
            }
            else    // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * precursor_mass_tolerance);
            }

            if (low_it == up_it) continue; // no matching precursor in data

            //add peaks for b and y ions with charge 1
            if (theo_spectrum.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              spectrum_generator.getSpectrum(theo_spectrum, candidate, 1);
              theo_spectrum.sortByPosition(); //sort by mz
            }

            for (; low_it != up_it; ++low_it)
            {
              const Size& scan_index = low_it->second;
              const MSSpectrum<Peak1D>& exp_spectrum = spectra[scan_index];

              double score = computeHyperScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);

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
    postProcessHits_(spectra, annotated_hits, protein_ids, peptide_ids, report_top_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide);
    progresslogger.endProgress();

    // annotate RNPxl related information to hits and create report
    vector<RNPxlReportRow> csv_rows = annotateRNPxlInformation_(spectra, peptide_ids, marker_ions_tolerance);

    // save report
    TextFile csv_file;
    csv_file.addLine(RNPxlReportRowHeader().getString("\t"));
    for (Size i = 0; i != csv_rows.size(); ++i)
    {
      csv_file.addLine(csv_rows[i].getString("\t"));
    }
    csv_file.store(out_csv);

    // write ProteinIdentifications and PeptideIdentifications to IdXML
    IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  RNPxlSearch tool;
  return tool.main(argc, argv);
}
