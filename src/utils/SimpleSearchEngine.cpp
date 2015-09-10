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
#include <OpenMS/CHEMISTRY/EnzymesDB.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
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

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

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
   - proper C-term N-term handling of terminal modifications that can be at every amino acid

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

   - make hyperscore scoring linear
   - single and multiple neutral loss spectra creation
*/

class SimpleSearchEngine :
    public TOPPBase
{
  public:
    SimpleSearchEngine() :
      TOPPBase("SimpleSearchEngine", "Annotates MS/MS spectra using SimpleSearchEngine.", false)
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

      registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
      registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

      StringList precursor_mass_tolerance_unit_valid_strings;
      precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
      precursor_mass_tolerance_unit_valid_strings.push_back("Da");

      registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
      setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

      registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, true);
      registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, true);

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

      vector<String> all_enzymes;
      EnzymesDB::getInstance()->getAllNames(all_enzymes);
      registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
      setValidStrings_("enzyme", all_enzymes);

      registerTOPPSubsection_("peptide", "Peptide Options");
      registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
      registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);

      registerTOPPSubsection_("report", "Reporting Options");
      registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);
    }

    vector<ResidueModification> getModifications_(StringList modNames)
    {
      vector<ResidueModification> modifications;

      // iterate over modification names and add to vector
      for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
      {
        String modification(*mod_it);
        modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
      }

      return modifications;
    }

    // check if for minimum size
    class HasInvalidPeptideLengthPredicate
    {
        public:
          explicit HasInvalidPeptideLengthPredicate(Size min_size)
            :min_size_(min_size)
          {
          }

          bool operator()(const AASequence& aas)
          {
            return (aas.size() < min_size_);
          }
      private:
          Size min_size_;
    };

    // spectrum must not contain 0 intensity peaks and must be sorted by m/z
    template <typename SpectrumType>
    static void deisotopeAndSingleChargeMSSpectrum(SpectrumType& in, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = true)
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

        for (Int q = max_charge; q >= min_charge; --q)   // important: test charge hypothesis from high to low
        {
          // try to extend isotopes from mono-isotopic peak
          // if extension larger then min_isopeaks possible:
          //   - save charge q in mono_isotopic_peak[]
          //   - annotate all isotopic peaks with feature number
          if (features[current_peak] == -1)   // only process peaks which have no assigned feature number
          {
            bool has_min_isopeaks = true;
            vector<Size> extensions;
            for (Size i = 0; i < max_isopeaks; ++i)
            {
              double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
              Size p = old_spectrum.findNearest(expected_mz);
              double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
              if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton)   // test for missing peak
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

    void postProcessHits_(const PeakMap& exp, const vector<vector<PeptideHit> >& peptide_hits, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids, Size top_hits)
    {
      for (vector<vector<PeptideHit> >::const_iterator pit = peptide_hits.begin(); pit != peptide_hits.end(); ++pit)
      {
        if (!pit->empty())
        {
          Size scan_index = pit - peptide_hits.begin();

          // create empty PeptideIdentification object and fill meta data
          PeptideIdentification pi;
          pi.setScoreType("hyperscore");
          pi.setHigherScoreBetter(true);
          pi.setRT(exp[scan_index].getRT());
          pi.setMZ(exp[scan_index].getPrecursors()[0].getMZ());
          pi.setHits(*pit);
          pi.assignRanks();
          peptide_ids.push_back(pi);
        }
      }

      IDFilter filter;

      // only store top n hits
      for (vector<PeptideIdentification>::iterator pids_it = peptide_ids.begin(); pids_it != peptide_ids.end(); ++pids_it)
      {
        PeptideIdentification& pi = *pids_it;
        PeptideIdentification temp_identification = pi;
        filter.filterIdentificationsByBestNHits(temp_identification, top_hits, pi);
      }

      // protein identifications (leave as is...)
      protein_ids = vector<ProteinIdentification>(1);
      protein_ids[0].setDateTime(DateTime::now());
      protein_ids[0].setSearchEngine("SimpleSearchEngine");
      protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    }

    ExitCodes main_(int, const char**)
    {
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);
      String in_mzml = getStringOption_("in");
      String in_db = getStringOption_("database");
      String out_idxml = getStringOption_("out");

      Int min_precursor_charge = getIntOption_("precursor:min_charge");
      Int max_precursor_charge = getIntOption_("precursor:max_charge");
      double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
      bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");

      double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
      bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

      StringList fixedModNames = getStringList_("modifications:fixed");
      set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

      Size peptide_min_size = getIntOption_("peptide:min_size");

      if (fixed_unique.size() != fixedModNames.size())
      {
        cout << "duplicate fixed modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      StringList varModNames = getStringList_("modifications:variable");
      set<String> var_unique(varModNames.begin(), varModNames.end());
      if (var_unique.size() != varModNames.size())
      {
        cout << "duplicate variable modification provided." << endl;
        return ILLEGAL_PARAMETERS;
      }

      vector<ResidueModification> fixedMods = getModifications_(fixedModNames);
      vector<ResidueModification> varMods = getModifications_(varModNames);
      Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

      Int report_top_hits = getIntOption_("report:top_hits");

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
          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
        }
      }

      // create spectrum generator
      TheoreticalSpectrumGenerator spectrum_generator;

      vector<vector<PeptideHit> > peptide_hits(spectra.size(), vector<PeptideHit>());

      progresslogger.startProgress(0, 1, "Load database from FASTA file...");
      FASTAFile fastaFile;
      vector<FASTAFile::FASTAEntry> fasta_db;
      fastaFile.load(in_db, fasta_db);
      progresslogger.endProgress();

      const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
      EnzymaticDigestion digestor;
      digestor.setEnzyme(getStringOption_("enzyme"));
      digestor.setMissedCleavages(missed_cleavages);

      progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

      // lookup for processed peptides. must be defined outside of omp section and synchronized
      set<std::string> processed_petides;

      // set minimum size of peptide after digestion
      HasInvalidPeptideLengthPredicate has_invalid_length(getIntOption_("peptide:min_size"));

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
      {
        IF_MASTERTHREAD
        {
          progresslogger.setProgress((SignedSize)fasta_index * NUMBER_OF_THREADS);
        }

        const AASequence& seq = AASequence::fromString(fasta_db[fasta_index].sequence);

        vector<AASequence> current_digest;
        digestor.digest(seq, current_digest);

        // c++ STL pattern for deleting entries from vector based on predicate evaluation
        current_digest.erase(std::remove_if(current_digest.begin(), current_digest.end(), has_invalid_length), current_digest.end());

        for (vector<AASequence>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
        {
          const std::string& s = cit->toUnmodifiedString();

          bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
          {
            if (processed_petides.find(s) != processed_petides.end())
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
            processed_petides.insert(s);
          }

          vector<AASequence> all_modified_peptides;

          // this critial section is because ResidueDB is not thread safe and new residues are created based on the PTMs
#ifdef _OPENMP
#pragma omp critical (residuedb_access)
#endif
          {
            ModifiedPeptideGenerator::applyFixedModifications(fixedMods.begin(), fixedMods.end(), *cit);
            ModifiedPeptideGenerator::applyVariableModifications(varMods.begin(), varMods.end(), *cit, max_variable_mods_per_peptide, all_modified_peptides);
          }

          for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
          {
            const AASequence& candidate = all_modified_peptides[mod_pep_idx];
            double current_peptide_mass = candidate.getMonoWeight();

            // determine MS2 precursors that match to the current peptide mass
            multimap<double, Size>::const_iterator low_it;
            multimap<double, Size>::const_iterator up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * current_peptide_mass * precursor_mass_tolerance * 1e-6);
            }
            else // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - 0.5 * precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + 0.5 * precursor_mass_tolerance);
            }

            if (low_it == up_it)
            {
              continue;     // no matching precursor in data
            }

            //create theoretical spectrum
            MSSpectrum<RichPeak1D> theo_spectrum = MSSpectrum<RichPeak1D>();

            //add peaks for b and y ions with charge 1
            spectrum_generator.getSpectrum(theo_spectrum, candidate, 1);

            //sort by mz
            theo_spectrum.sortByPosition();

            for (; low_it != up_it; ++low_it)
            {
              const Size& scan_index = low_it->second;
              const MSSpectrum<Peak1D>& exp_spectrum = spectra[scan_index];

              double score = computeHyperScore(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);

              // no hit
              if (score < 1e-16)
              {
                continue;
              }

              PeptideHit hit;
              hit.setSequence(candidate);
              hit.setCharge(exp_spectrum.getPrecursors()[0].getCharge());
              hit.setScore(score);
#ifdef _OPENMP
#pragma omp critical (peptide_hits_access)
#endif
              {
                peptide_hits[scan_index].push_back(hit);
              }
            }
          }
        }
      }
      progresslogger.endProgress();

      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      progresslogger.startProgress(0, 1, "Post-processing PSMs...");
      postProcessHits_(spectra, peptide_hits, protein_ids, peptide_ids, report_top_hits);
      progresslogger.endProgress();

      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

      return EXECUTION_OK;
    }

};

int main(int argc, const char** argv)
{
  SimpleSearchEngine tool;
  return tool.main(argc, argv);
}
