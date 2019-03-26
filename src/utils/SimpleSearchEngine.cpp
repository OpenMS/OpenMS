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

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/HyperScore.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/METADATA/SpectrumSettings.h>

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

class SimpleSearchEngine :
    public TOPPBase
{
  /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
  struct AnnotatedHit
  {
    StringView sequence;
    SignedSize peptide_mod_index; // enumeration index of the non-RNA peptide modification
    double score = 0; // main score
    std::vector<PeptideHit::PeakAnnotation> fragment_annotations;

    static bool hasBetterScore(const AnnotatedHit& a, const AnnotatedHit& b)
    {
      return a.score > b.score;
    }
  };

  public:
    SimpleSearchEngine() :
      TOPPBase("SimpleSearchEngine", 
        "Annotates MS/MS spectra using SimpleSearchEngine.", 
        false)
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

      registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
      registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

      StringList precursor_mass_tolerance_unit_valid_strings;
      precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
      precursor_mass_tolerance_unit_valid_strings.push_back("Da");

      registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
      setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

      registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, true);
      registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, true);

      // consider one before annotated monoisotopic peak and the annotated one
      IntList isotopes = {0, 1};
      registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)", false, false);

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
      ProteaseDB::getInstance()->getAllNames(all_enzymes);
      registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
      setValidStrings_("enzyme", all_enzymes);

      registerTOPPSubsection_("peptide", "Peptide Options");
      registerIntOption_("peptide:min_size", "<num>", 7, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
      registerIntOption_("peptide:max_size", "<num>", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).", false, true);
      registerIntOption_("peptide:missed_cleavages", "<num>", 1, "Number of missed cleavages.", false, false);
      registerStringOption_("peptide:motif", "<regex>", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.", false);

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
        modifications.push_back(*ModificationsDB::getInstance()->getModification(modification));
      }

      return modifications;
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
        Deisotoper::deisotopeAndSingleCharge(exp[exp_index], 
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
          1, 3,   // min / max charge 
          false,  // keep only deisotoped
          3, 10,  // min / max isopeaks 
          true);  // convert fragment m/z to mono-charge

        // remove noise
        window_mower_filter.filterPeakSpectrum(exp[exp_index]);
        nlargest_filter.filterPeakSpectrum(exp[exp_index]);

        // sort (nlargest changes order)
        exp[exp_index].sortByPosition();
      }
    }

    void postProcessHits_(const PeakMap& exp, 
      vector<vector<AnnotatedHit> >& annotated_hits, 
      vector<ProteinIdentification>& protein_ids, 
      vector<PeptideIdentification>& peptide_ids, 
      Size top_hits,
      const vector<ResidueModification>& fixed_modifications, 
      const vector<ResidueModification>& variable_modifications, 
      Size max_variable_mods_per_peptide)
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
    protein_ids[0].setSearchEngine("SimpleSearchEngine");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    search_parameters.charges = String(getIntOption_("precursor:min_charge")) + ":" + String(getIntOption_("precursor:max_charge"));

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = getStringList_("modifications:fixed");
    search_parameters.variable_modifications = getStringList_("modifications:variable");
    search_parameters.missed_cleavages = getIntOption_("peptide:missed_cleavages");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme"));
    protein_ids[0].setSearchParameters(search_parameters);
  }

    ExitCodes main_(int, const char**) override
    {
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);
      String in_mzml = getStringOption_("in");
      String in_db = getStringOption_("database");
      String out_idxml = getStringOption_("out");
      const String peptide_motif = getStringOption_("peptide:motif");      
      boost::regex peptide_motif_regex(peptide_motif);

      Int min_precursor_charge = getIntOption_("precursor:min_charge");
      Int max_precursor_charge = getIntOption_("precursor:max_charge");
      double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
      bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
      IntList precursor_isotopes = getIntList_("precursor:isotopes");

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

      vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
      vector<ResidueModification> variable_modifications = getModifications_(varModNames);
      Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

      size_t top_hits = static_cast<size_t>(getIntOption_("report:top_hits"));

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

          // calculate precursor mass (optionally corrected for misassignment) and map it to MS scan index
          for (int isotope_number : precursor_isotopes)
          {
            double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

            // correct for monoisotopic misassignments of the precursor annotation
            if (isotope_number != 0) { precursor_mass -= isotope_number * Constants::C13C12_MASSDIFF_U; }

            multimap_mass_2_scan_index.insert(make_pair(precursor_mass, scan_index));
          }
        }
      }

      // create spectrum generator
      TheoreticalSpectrumGenerator spectrum_generator;
      Param param(spectrum_generator.getParameters());
      param.setValue("add_first_prefix_ion", "true");
      param.setValue("add_metainfo", "true");
      spectrum_generator.setParameters(param);

      // preallocate storage for PSMs
      vector<vector<AnnotatedHit> > annotated_hits(spectra.size(), vector<AnnotatedHit>());
      for (auto & a : annotated_hits) { a.reserve(2 * top_hits); }

#ifdef _OPENMP
      // we want to do locking at the spectrum level so we get good parallelisation 
      vector<omp_lock_t> annotated_hits_lock(annotated_hits.size());
      for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_init_lock(&(annotated_hits_lock[i])); }
#endif

      progresslogger.startProgress(0, 1, "Load database from FASTA file...");
      FASTAFile fastaFile;
      vector<FASTAFile::FASTAEntry> fasta_db;
      fastaFile.load(in_db, fasta_db);
      progresslogger.endProgress();

      const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
      ProteaseDigestion digestor;
      digestor.setEnzyme(getStringOption_("enzyme"));
      digestor.setMissedCleavages(missed_cleavages);

      progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

      // lookup for processed peptides. must be defined outside of omp section and synchronized
      set<StringView> processed_petides;

      // set minimum / maximum size of peptide after digestion
      Size min_peptide_length = getIntOption_("peptide:min_size");
      Size max_peptide_length = getIntOption_("peptide:max_size");
      Size count_proteins(0), count_peptides(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
      {
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_proteins;

        IF_MASTERTHREAD
        {
          progresslogger.setProgress(count_proteins);
        }

        vector<StringView> current_digest;
        digestor.digestUnmodified(fasta_db[fasta_index].sequence, current_digest, min_peptide_length, max_peptide_length);

        for (auto const & c : current_digest)
        { 
          const String current_peptide = c.getString();
          if (current_peptide.find_first_of("XBZ") != std::string::npos) { continue; }

          // if a peptide motif is provided skip all peptides without match
          if (!peptide_motif.empty() && !boost::regex_match(current_peptide, peptide_motif_regex)) { continue; }          
        
          bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
          {
            // peptide (and all modified variants) already processed so skip it
            if (processed_petides.find(c) != processed_petides.end())
            {
              already_processed = true;
            }
          }

          // skip peptides that have already been processed
          if (already_processed) { continue; }

#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
          {
            processed_petides.insert(c);
          }

#ifdef _OPENMP
#pragma omp atomic
#endif
          ++count_peptides;

          vector<AASequence> all_modified_peptides;

          // this critial section is because ResidueDB is not thread safe and new residues are created based on the PTMs
#ifdef _OPENMP
#pragma omp critical (residuedb_access)
#endif
          {
            AASequence aas = AASequence::fromString(current_peptide);
            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
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

            // no matching precursor in data
            if (low_it == up_it) { continue; }

            // create theoretical spectrum
            PeakSpectrum theo_spectrum;

            // add peaks for b and y ions with charge 1
            spectrum_generator.getSpectrum(theo_spectrum, candidate, 1, 1);

            // sort by mz
            theo_spectrum.sortByPosition();

            for (; low_it != up_it; ++low_it)
            {
              const Size& scan_index = low_it->second;
              const PeakSpectrum& exp_spectrum = spectra[scan_index];
              // const int& charge = exp_spectrum.getPrecursors()[0].getCharge();
              const double& score = HyperScore::compute(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum);

              if (score == 0) { continue; } // no hit?

              // add peptide hit
              AnnotatedHit ah;
              ah.sequence = c;
              ah.peptide_mod_index = mod_pep_idx;
              ah.score = score;

#ifdef _OPENMP
              omp_set_lock(&(annotated_hits_lock[scan_index]));
              {
#endif
                annotated_hits[scan_index].push_back(ah);

                // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                if (annotated_hits[scan_index].size() >= 2 * top_hits)
                {
                  std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + top_hits, annotated_hits[scan_index].end(), AnnotatedHit::hasBetterScore);
                  annotated_hits[scan_index].resize(top_hits); 
                }
#ifdef _OPENMP
              }
              omp_unset_lock(&(annotated_hits_lock[scan_index]));
#endif
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
      postProcessHits_(spectra, 
        annotated_hits, 
        protein_ids, 
        peptide_ids, 
        top_hits,
        fixed_modifications, 
        variable_modifications, 
        max_variable_mods_per_peptide
        );
      progresslogger.endProgress();

      // add meta data on spectra file
      StringList ms_runs;
      spectra.getPrimaryMSRunPath(ms_runs);
      protein_ids[0].setPrimaryMSRunPath(ms_runs);

      // reindex peptides to proteins
      PeptideIndexing indexer;
      Param param_pi = indexer.getParameters();
      param_pi.setValue("decoy_string", "DECOY_");
      param_pi.setValue("decoy_string_position", "prefix");
      param_pi.setValue("enzyme:name", getStringOption_("enzyme"));
      param_pi.setValue("enzyme:specificity", "full");
      param_pi.setValue("missing_decoy_action", "silent");
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

      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

#ifdef _OPENMP
      // free locks
      for (size_t i = 0; i != annotated_hits_lock.size(); i++) { omp_destroy_lock(&(annotated_hits_lock[i])); }
#endif

      return EXECUTION_OK;
    }

};

int main(int argc, const char** argv)
{
  SimpleSearchEngine tool;
  return tool.main(argc, argv);
}
