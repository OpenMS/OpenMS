// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>

#include <map>
#include <algorithm>

#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

using namespace OpenMS;
using namespace std;
using namespace boost::unordered;

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

struct PeptideHitSequenceLessComparator
{
    bool operator()(const PeptideHit& a, const PeptideHit& b)
    {
      if (a.getSequence().toString() < b.getSequence().toString()) return true;
      return false;
    }
};

class RNPxlSearch
    : public TOPPBase
{
  public:
    RNPxlSearch() :
      TOPPBase("RNPxlSearch", "Annotates MS/MS spectra using RNPxlSearch.", false)
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

      // mass tolerance
      registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor mass tolerance", false);

      StringList precursor_mass_tolerance_unit_valid_strings;
      precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
      precursor_mass_tolerance_unit_valid_strings.push_back("Da");
      registerStringOption_("precursor_mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
      setValidStrings_("precursor_mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

      registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance", false);

      StringList fragment_mass_tolerance_unit_valid_strings;
      fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
      fragment_mass_tolerance_unit_valid_strings.push_back("Da");

      registerStringOption_("fragment_mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
      setValidStrings_("fragment_mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

      // modifications
      vector<String> all_mods;
      ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
      registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
      setValidStrings_("fixed_modifications", all_mods);

      registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
      setValidStrings_("variable_modifications", all_mods);

      registerIntOption_("peptide_max_var_mods", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);
      registerIntOption_("top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);

      registerIntOption_("missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

      // RNPxl specific
      registerIntOption_("length", "", 4, "Oligonucleotide maximum length.", false);

      registerStringOption_("sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)", false);

      StringList target_nucleotides;
      target_nucleotides.push_back("A=C10H14N5O7P");
      target_nucleotides.push_back("C=C9H14N3O8P");
      target_nucleotides.push_back("G=C10H14N5O8P");
      target_nucleotides.push_back("U=C9H13N2O9P");

      registerStringList_("target_nucleotides", "", target_nucleotides, "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG", false, false);

      StringList mapping;
      mapping.push_back("A->A");
      mapping.push_back("C->C");
      mapping.push_back("G->G");
      mapping.push_back("U->U");

      registerStringList_("mapping", "", mapping, "format: source->target e.g. A->A, ..., U->U, U->X", false, false);

      StringList restrictions;
      restrictions.push_back("A=0");
      restrictions.push_back("C=0");
      restrictions.push_back("U=0");
      restrictions.push_back("G=0");

      registerStringList_("restrictions", "", restrictions, "format: target nucleotide=min_count: e.g U=1 if at least one U must be in the generated sequence.", false, false);

      StringList modifications;
      modifications.push_back("-H2O");
      modifications.push_back("");
      modifications.push_back("-H2O-HPO3");
      modifications.push_back("-HPO3");
      modifications.push_back("-H2O+HPO3");
      modifications.push_back("+HPO3");

      registerStringList_("modifications", "", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3", false, false);

      registerFlag_("CysteineAdduct", "Use this flag if the +152 adduct is expected.");
    }

    vector<ResidueModification> getModifications_(StringList modNames)
    {
      vector<ResidueModification> modifications;

      // iterate over modification names and add to vector
      for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
      {
        String modification(*mod_it);
        modifications.push_back( ModificationsDB::getInstance()->getModification(modification));
      }

      return modifications;
    }

    static bool HasInvalidPeptideLengthPredicate_(const AASequence& aas)
    {
      bool has_invalid_length = aas.size() < 7;
      return has_invalid_length;
    }

    // spectrum must not contain 0 intensity peaks and must be sorted by mz
    template <typename SpectrumType>
    void deisotopeMSSpectrum(SpectrumType& in, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 2, Size max_isopeaks = 10)
    {
      if (in.empty())
      {
        return;
      }

      static const Int UNASSIGNED_FEAUTURE_ID = -1;
      SpectrumType old_spectrum = in;

      // annotate every peak not assigned to a feature
      vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
      vector<Int> features(old_spectrum.size(), UNASSIGNED_FEAUTURE_ID);
      Int feature_number = 0;

      // determine charge charge seeds and extend them
      for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
      {
        double current_mz = old_spectrum[current_peak].getPosition()[0];

        for (Int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
        {
          // skip peaks already assigned to a feature
          if (features[current_peak] != UNASSIGNED_FEAUTURE_ID)
          {
            continue;
          }

          bool has_min_isopeaks = true;
          vector<Size> extensions;
          for (Size i = 0; i < max_isopeaks; ++i)
          {
            double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
            Size p = old_spectrum.findNearest(expected_mz);
            double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
            if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton) // test for missing peak
            {
              // if current isotope is missing and smaller than expected number of observed isotopic peaks then mark as incomplete
              if (i < min_isopeaks)
              {
                has_min_isopeaks = false;
              }
              break;
            } else
            {
              // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis that intensities are decreasing
              Size n_extensions = extensions.size();
              if (n_extensions != 0)
              {
                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions-1]].getIntensity())
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

      in.clear(false);
      for (Size i = 0; i != old_spectrum.size(); ++i)
      {
        if (keep_only_deisotoped)
        {
          if (mono_isotopic_peak[i] != 0)
          {
            in.push_back(old_spectrum[i]);
          }
        } else  // deisotoped and unassigned
        {
          if (mono_isotopic_peak[i] != 0 || features[i] < 0)
          {
            in.push_back(old_spectrum[i]);
          }
        }
      }

      /*
    cout << "monoisotopic peak charges: " << endl;
    for (Size i = 0; i != in.size(); ++i)
    {
        cout << mono_isotopic_peak[i];
    }
    cout <<  endl;

    cout << "feature numbers for each peak: " << endl;
    for (Size i = 0; i != in.size(); ++i)
    {
        if (features[i] != -1)
        {
            cout << features[i] << " " << in[i].getMZ() << " . ";
        } else
        {
            cout << ". ";
        }
    }
    cout <<  endl;
    */

      in.sortByPosition();
    }

    double logfactorial(UInt x)
    {
      int y;

      if (x < 2)
        return 1;
      else
      {
        double z = 0;
        for (y = 2; y<=x; y++)
        {
          z = log(y)+z;
        }

        return (z);
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
          } else
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
      } else
      {
        return 0;
      }
    }

    ExitCodes main_(int, const char**)
    {
      String in_mzml = getStringOption_("in");
      String in_db = getStringOption_("database");
      String out_idxml = getStringOption_("out");

      double precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
      double fragment_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");

      bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor_mass_tolerance_unit") == "ppm");
      bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment_mass_tolerance_unit") == "ppm");

      StringList fixedModNames = getStringList_("fixed_modifications");
      set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

      if (fixed_unique.size() != fixedModNames.size())
      {
        cout << "duplicate fixed modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      StringList varModNames = getStringList_("variable_modifications");
      set<String> var_unique(varModNames.begin(), varModNames.end());
      if (var_unique.size() != varModNames.size())
      {
        cout << "duplicate variable modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
      StringList target_nucleotides = getStringList_("target_nucleotides");

      // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
      StringList mappings = getStringList_("mapping");

      // string format: target,min_count: e.g "X=1" if at least one tU must be in the generated sequence.
      // All target nucleotides must be included. X=0 -> disable restriction
      StringList restrictions = getStringList_("restrictions");

      StringList modifications = getStringList_("modifications");

      String sequence_restriction = getStringOption_("sequence");

      Int max_nucleotide_length = getIntOption_("length");

      bool cysteine_adduct = getFlag_("CysteineAdduct");

      RNPxlModificationMassesResult mm = RNPxlModificationsGenerator::initModificationMassesRNA(target_nucleotides, mappings, restrictions, modifications, sequence_restriction, cysteine_adduct, max_nucleotide_length);
      mm.mod_masses[""] = 0;  // insert "null" modification otherwise unmodified peptide will not be searched
      mm.mod_combinations[""].insert("none");

      // load MS2 map
      PeakMap exp;
      MzMLFile f;
      f.setLogType(log_type_);

      PeakFileOptions options;
      options.clearMSLevels();
      options.addMSLevel(2);
      f.getOptions() = options;
      f.load(in_mzml, exp);
      exp.sortSpectra(true);

      cout << "filtering spectra ... " << endl;
      // filter MS2 map
      // remove 0 intensities
      ThresholdMower threshold_mower_filter;
      threshold_mower_filter.filterPeakMap(exp);

      Normalizer normalizer;
      normalizer.filterPeakMap(exp);

      exp.sortSpectra(true);
      for (PeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
      {
        deisotopeMSSpectrum(*it, 1, 3, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 2, 10);
      }

      if (getIntOption_("debug") >= 1)
      {
        f.store(String("deisotoped_" + String(in_mzml)), exp);
      }

      // remove noise
      WindowMower window_mower_filter;
      Param filter_param = window_mower_filter.getParameters();
      filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
      filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
      filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
      window_mower_filter.setParameters(filter_param);
      window_mower_filter.filterPeakMap(exp);

      // conservativly restrict number of peaks
      NLargest nlargest_filter = NLargest(400);
      nlargest_filter.filterPeakMap(exp);
      exp.sortSpectra(true);

      // build multimap of precursor mass to scan index
      multimap<double, Size> multimap_mass_2_scan_index;
      for (PeakMap::ConstIterator s_it = exp.begin(); s_it != exp.end(); ++s_it)
      {
        int scan_index = s_it - exp.begin();
        vector<Precursor> precursor = s_it->getPrecursors();

        if (precursor.size() == 1 && s_it->size() > 5)  // there should only one precursor and MS2 should contain at least 5 peaks
        {
          int precursor_charge = precursor[0].getCharge();
          double precursor_mz = precursor[0].getMZ();
          double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;
          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        }
      }

      // create spectrum generator
      TheoreticalSpectrumGenerator spectrum_generator;
      vector<PeptideIdentification> peptide_ids;
      vector<vector<PeptideHit> > peptide_hits(exp.size(), vector<PeptideHit>());

      // load database from FASTA file
      cout << "loading database ... " << endl;
      FASTAFile fastaFile;
      vector<FASTAFile::FASTAEntry> fasta_db;
      fastaFile.load(in_db, fasta_db);

      // digest database
      const Size missed_cleavages = getIntOption_("missed_cleavages");
      EnzymaticDigestion digestor;
      digestor.setEnzyme(EnzymaticDigestion::ENZYME_TRYPSIN);
      digestor.setMissedCleavages(missed_cleavages);

      cout << "digesting database, apply modifications and score peptides against spectra... " << endl;
      vector<ResidueModification> fixedMods = getModifications_(fixedModNames);
      vector<ResidueModification> varMods = getModifications_(varModNames);
      Size max_variable_mods_per_peptide = getIntOption_("peptide_max_var_mods");

      boost::unordered_set<Size> cached_peptides;
      for (vector<FASTAFile::FASTAEntry>::const_iterator fasta_it = fasta_db.begin(); fasta_it != fasta_db.end(); ++fasta_it)
      {
        vector<AASequence> current_digest;
        digestor.digest(AASequence::fromString(fasta_it->sequence), current_digest);
        // c++ STL pattern for deleting entries from vector based on predicate evaluation
        current_digest.erase(std::remove_if(current_digest.begin(), current_digest.end(), HasInvalidPeptideLengthPredicate_), current_digest.end());
        // make unique
        set<AASequence> tmp(current_digest.begin(), current_digest.end());

        // look up if peptide has already been searched
        for (set<AASequence>::iterator set_it = tmp.begin(); set_it != tmp.end(); )  // note that this is the STL pattern for deleting while iterating a set so don't change operator order
        {
          Size string_hash = boost::hash<std::string>()(set_it->toUnmodifiedString());
          if (cached_peptides.find(string_hash) != cached_peptides.end())
          {
            tmp.erase(set_it++);
          } else
          {
            ++set_it;
            cached_peptides.insert(string_hash);  // ok to add already here as s is unique for this digested protein
          }
        }

        // assign unprocessed peptides
        current_digest.assign(tmp.begin(), tmp.end());

        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixedMods.begin(), fixedMods.end(), current_digest.begin(), current_digest.end());
        ModifiedPeptideGenerator::applyVariableModifications(varMods.begin(), varMods.end(), current_digest.begin(), current_digest.end(), max_variable_mods_per_peptide, all_modified_peptides);

        for (vector<AASequence>::const_iterator mod_peps_it = all_modified_peptides.begin(); mod_peps_it != all_modified_peptides.end(); ++mod_peps_it)
        {
          const AASequence& candidate = *mod_peps_it;
          double current_peptide_mass_without_RNA = candidate.getMonoWeight();

          //create empty theoretical spectrum
          MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();

          // iterate over all RNA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          for (std::map<String, double>::const_iterator rna_mod_it = mm.mod_masses.begin(); rna_mod_it != mm.mod_masses.end(); ++rna_mod_it)
          {
            double current_peptide_mass = current_peptide_mass_without_RNA + rna_mod_it->second; // add RNA mass

            // determine MS2 precursors that match to the current peptide mass
            multimap<double, Size>::const_iterator low_it;
            multimap<double, Size>::const_iterator up_it;

            if (precursor_mass_tolerance_unit_ppm)  // ppm
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + current_peptide_mass * precursor_mass_tolerance * 1e-6);
            } else  // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + precursor_mass_tolerance);
            }

            if (low_it == up_it) continue; // no matching precursor in data

            //add peaks for b and y ions with charge 1
            if (theo_spectrum.empty())  // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              spectrum_generator.getSpectrum(theo_spectrum, candidate, 1);
              theo_spectrum.sortByPosition(); //sort by mz
            }

            for (; low_it != up_it; ++low_it)
            {
              const Size& scan_index = low_it->second;
              const MSSpectrum<Peak1D>& exp_spectrum = exp[scan_index];

              double score = computeHyperScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);

              // no hit
              if (score < 1e-16)
              {
                continue;
              }

              // add peptide hit
              PeptideHit hit;
              hit.setSequence(candidate);
              hit.setMetaValue(String("RNA"), *mm.mod_combinations.at(rna_mod_it->first).begin());  // return first nucleotide formula matching current empirical formula and mass
              hit.setCharge(exp_spectrum.getPrecursors()[0].getCharge());
              hit.setScore(score);
              peptide_hits[scan_index].push_back(hit);
            }
          }
        }
      }

      // annotate with RT and MZ and sort by rank
      for (vector< vector<PeptideHit> >::const_iterator pit = peptide_hits.begin(); pit != peptide_hits.end(); ++pit)
      {
        if (!pit->empty())
        {
          Size scan_index = pit - peptide_hits.begin();

          // create empty PeptideIdentification object and fill meta data
          PeptideIdentification pi;
          pi.setScoreType("hyperscore");
          pi.setHigherScoreBetter(true);
          pi.setMetaValue("RT", exp[scan_index].getRT());
          pi.setMetaValue("MZ", exp[scan_index].getPrecursors()[0].getMZ());
          pi.setHits(*pit);

          // sort by score
          pi.assignRanks();

          peptide_ids.push_back(pi);
        }
      }

      // calculate delta ion score
      for (vector<PeptideIdentification>::iterator pids_it = peptide_ids.begin(); pids_it != peptide_ids.end(); ++pids_it)
      {
        vector<PeptideHit>& pit = pids_it->getHits();
        if (pit.size() == 1)
        {
          PeptideHit& best = pids_it->getHits()[0];
          best.setMetaValue(String("D-Score"), 1000.0);  // set to maximum if there is no second peptide
        } else if (pit.size() >= 2)
        {
          PeptideHit& best = pids_it->getHits()[0];
          PeptideHit& second = pids_it->getHits()[1];
          best.setMetaValue(String("D-Score"), best.getScore() - second.getScore());
        }
      }

      // filter for top n hits
      IDFilter filter;
      Int top_hits = getIntOption_("top_hits");

      for (vector<PeptideIdentification>::iterator pids_it = peptide_ids.begin(); pids_it != peptide_ids.end(); ++pids_it)
      {
        PeptideIdentification& pi = *pids_it;
        PeptideIdentification temp_identification = pi;
        filter.filterIdentificationsByBestNHits(temp_identification, top_hits, pi);
      }

      // protein identifications (leave as is...)
      vector<ProteinIdentification> protein_ids(1);
      protein_ids[0].setDateTime(DateTime::now());
      protein_ids[0].setSearchEngine("RNPxlSearch");
      protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

      return EXECUTION_OK;
    }
};

int main(int argc, const char ** argv)
{
  RNPxlSearch tool;
  return tool.main(argc, argv);
}

