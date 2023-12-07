// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/PeptideSearchEngineAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

#include <algorithm>
#include <map>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  PeptideSearchEngineAlgorithm::PeptideSearchEngineAlgorithm() :
    DefaultParamHandler("PeptideSearchEngineAlgorithm"),
    ProgressLogger()
  {
    defaults_.setValue("precursor:mass_tolerance", 10.0, "+/- tolerance for precursor mass.");

    std::vector<std::string> precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered.");

    defaults_.setSectionDescription("precursor", "Precursor (Parent Ion) Options");

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    defaults_.setValue("precursor:isotopes", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)");

    defaults_.setValue("fragment:mass_tolerance", 10.0, "Fragment mass tolerance");

    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);


    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");    

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));

    defaults_.setValue("decoys", "false", "Should decoys be generated?");
    defaults_.setValidStrings("decoys", {"true","false"} );

    defaults_.setValue("annotate:PSM",  std::vector<std::string>{"ALL"}, "Annotations added to each PSM.");
    defaults_.setValidStrings("annotate:PSM", 
      std::vector<std::string>{
        "ALL",
        Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM, 
        Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
        Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
        Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION}
      );

    defaults_.setSectionDescription("annotate", "Annotation Options");

    defaults_.setValue("peptide:min_size", 7, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");
    defaults_.setValue("peptide:motif", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.");
    defaults_.setSectionDescription("peptide", "Peptide Options");

    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setSectionDescription("report", "Reporting Options");

    defaultsToParam_();
  }

  void PeptideSearchEngineAlgorithm::updateMembers_()
  {
    precursor_mass_tolerance_ = param_.getValue("precursor:mass_tolerance");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();

    precursor_min_charge_ = param_.getValue("precursor:min_charge");
    precursor_max_charge_ = param_.getValue("precursor:max_charge");

    precursor_isotopes_ = param_.getValue("precursor:isotopes");

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");

    fragment_mass_tolerance_unit_ = param_.getValue("fragment:mass_tolerance_unit").toString();

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    set<String> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());
    if (fixed_unique.size() != modifications_fixed_.size())
    {
      OPENMS_LOG_WARN << "Duplicate fixed modification provided. Making them unique." << endl;
      modifications_fixed_.assign(fixed_unique.begin(), fixed_unique.end());
    }    

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<String> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      OPENMS_LOG_WARN << "Duplicate variable modification provided. Making them unique." << endl;
      modifications_variable_.assign(var_unique.begin(), var_unique.end());
    }

    modifications_max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    enzyme_ = param_.getValue("enzyme").toString();

    peptide_min_size_ = param_.getValue("peptide:min_size");
    peptide_max_size_ = param_.getValue("peptide:max_size");
    peptide_missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_motif_ = param_.getValue("peptide:motif").toString(); // TODO: remove unused parameters

    report_top_hits_ = param_.getValue("report:top_hits");

    decoys_ = param_.getValue("decoys") == "true";
    annotate_psm_ = ListUtils::toStringList<std::string>(param_.getValue("annotate:PSM"));

  }

  // static
  void PeptideSearchEngineAlgorithm::preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
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

#pragma omp parallel for default(none) shared(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, window_mower_filter, nlargest_filter)
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

void PeptideSearchEngineAlgorithm::postProcessHits_(const PeakMap& exp, 
      std::vector<std::vector<PeptideSearchEngineAlgorithm::AnnotatedHit_> >& annotated_hits, 
      std::vector<ProteinIdentification>& protein_ids, 
      std::vector<PeptideIdentification>& peptide_ids, 
      Size top_hits,
      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
      Size max_variable_mods_per_peptide,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const String& precursor_mass_tolerance_unit_ppm,
      const String& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const String& enzyme,
      const String& database_name) const
  {
    // remove all but top n scoring TODO: use two parameters to distinguish between number of reported peptides and number of pre-scored peptides
#pragma omp parallel for default(none) shared(annotated_hits, top_hits)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
      annotated_hits.shrink_to_fit();
    }

    bool annotation_precursor_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_fragment_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_prefix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_suffix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION) != annotate_psm_.end();

    // "ALL" adds all annotations
    if (std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end())
    {
      annotation_precursor_error_ppm = true;
      annotation_fragment_error_ppm = true;
      annotation_prefix_fraction = true;
      annotation_suffix_fraction = true;
    }

#pragma omp parallel for
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spec = exp[scan_index];
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setSpectrumReference( spec.getNativeID());
        pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        pi.setScoreType("hyperscore");
        pi.setHigherScoreBetter(true);
        double mz = spec.getPrecursors()[0].getMZ();
        pi.setRT(spec.getRT());
        pi.setMZ(mz);
        Size charge = spec.getPrecursors()[0].getCharge();

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        for (const auto& ah : annotated_hits[scan_index])
        {
          PeptideHit ph;
          ph.setCharge(charge);
          ph.setScore(ah.score);
          ph.setSequence(ah.sequence);

          if (annotation_fragment_error_ppm)
          {
            TheoreticalSpectrumGenerator tsg;
            vector<pair<Size, Size> > alignment;
            MSSpectrum theoretical_spec;
            tsg.getSpectrum(theoretical_spec, ah.sequence, 1, std::min((int)charge - 1, 2));
            SpectrumAlignment sa;
            sa.getSpectrumAlignment(alignment, theoretical_spec, spec);

            vector<double> err;
            for (const auto& match : alignment)
            {
              double fragment_error = fabs(Math::getPPM(spec[match.second].getMZ(), theoretical_spec[match.first].getMZ()));
              err.push_back(fragment_error);
            }
            double median_ppm_error(0);
            if (!err.empty()) { median_ppm_error = Math::median(err.begin(), err.end(), false); }
            ph.setMetaValue(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM, median_ppm_error);
          }

          if (annotation_precursor_error_ppm)
          {
            double theo_mz = ah.sequence.getMZ(charge);
            double ppm_difference = Math::getPPM(mz, theo_mz);
            ph.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, ppm_difference);
          }

          if (annotation_prefix_fraction)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION, ah.prefix_fraction);
          }

          if (annotation_suffix_fraction)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION, ah.suffix_fraction);
          }

          // store PSM
          phs.push_back(ph);
        }
        pi.setHits(phs);
        pi.assignRanks();

#pragma omp critical (peptide_ids_access)
        {
          //clang-tidy: seems to be a false-positive in combination with omp
          peptide_ids.push_back(std::move(pi));
        }
      }
    }

#ifdef _OPENMP
    // we need to sort the peptide_ids by scan_index in order to have the same output in the idXML-file
    if (omp_get_max_threads() > 1)
    {
      std::sort(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return a.getMetaValue("scan_index") < b.getMetaValue("scan_index");
      });
    }
#endif

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("PeptideSearchEngine");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    DateTime now = DateTime::now();
    String identifier("PSE_" + now.get());
    protein_ids[0].setIdentifier(identifier);
    for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges = String(precursor_min_charge) + ":" + String(precursor_max_charge);

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
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
    // note: precursor error is calculated by percolator itself
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));

    search_parameters.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    protein_ids[0].setSearchParameters(std::move(search_parameters));
  }

  PeptideSearchEngineAlgorithm::ExitCodes PeptideSearchEngineAlgorithm::search(const String& in_mzML, const String& in_db, vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids) const
  {
    bool precursor_mass_tolerance_unit_ppm = (precursor_mass_tolerance_unit_ == "ppm");
    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    // load MS2 map
    PeakMap spectra;
    FileHandler f;

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.loadExperiment(in_mzML, spectra, {FileTypes::MZML});
    spectra.sortSpectra(true);

    startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
    endProgress();

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param(spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    spectrum_generator.setParameters(param);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(report_top_hits_); }

    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(in_db, fasta_db);

    // generate decoy protein sequences by reversing them
    if (decoys_)
    {
      startProgress(0, 1, "Generate decoys...");
      DecoyGenerator decoy_generator;

      // append decoy proteins
      const size_t old_size = fasta_db.size();
      fasta_db.reserve(fasta_db.size() * 2);
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = fasta_db[i];
        e.sequence = decoy_generator.reversePeptides(AASequence::fromString(e.sequence), enzyme_).toString();
        e.identifier = "DECOY_" + e.identifier;
        fasta_db.push_back(std::move(e));
      }
      // randomize order of targets and decoys to introduce no global bias in the case that
      // many targets have the same score as their decoy. (As we always take the first best scoring one)
      Math::RandomShuffler shuffler;
      shuffler.portable_random_shuffle(fasta_db.begin(), fasta_db.end());
      endProgress();
    }
    
    // build fragment index
    //TODO: Pass all the other parameters from this class to FragmentIndex
    //TODO: Can we do it with p.setValue or is there a more sophisticated way?
    startProgress(0, 1, "Building fragment index...");    
    FragmentIndex fragment_index_;
    auto p = fragment_index_.getParameters();
    p.setValue("max_processed_hits", report_top_hits_);
    p.setValue("fragment_min_mz", param_.getValue("fragment:min_mz"));
    p.setValue("fragment_max_mz", param_.getValue("fragment:max_mz"));

    fragment_index_.setParameters(p);
    fragment_index_.build(fasta_db);
    endProgress();

    startProgress(0, spectra.size(), "Scoring peptide models against spectra...");
    size_t count_spectra{};

#pragma omp parallel for schedule(static) default(none) shared(annotated_hits, count_spectra, fragment_index_, spectrum_generator, fasta_db, precursor_mass_tolerance_unit_ppm, fragment_mass_tolerance_unit_ppm, spectra)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {

      #pragma omp atomic
      ++count_spectra;

      IF_MASTERTHREAD
      {
        setProgress(count_spectra);
      }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      FragmentIndex::SpectrumMatchesTopN top_sms;
      fragment_index_.querySpectrum(exp_spectrum, top_sms); // TODO: expose top N as argument here and use report_top_hits_

      for (const auto& sms : top_sms.hits_)
      {
        FragmentIndex::Peptide sms_pep = fragment_index_.getPeptides()[sms.peptide_idx_];
        pair<size_t, size_t> candidate_snippet = sms_pep.sequence_;
        AASequence unmod_candidate = AASequence::fromString(fasta_db[sms_pep.protein_idx].sequence.substr(candidate_snippet.first, candidate_snippet.second));
        AASequence mod_candidate;
        //reapply modifications.
        if (!(modifications_variable_.empty() && modifications_fixed_.empty()))
        {
          vector<AASequence> mod_candidates;
          ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
          ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, unmod_candidate);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, unmod_candidate, modifications_max_variable_mods_per_peptide_, mod_candidates);
          mod_candidate = mod_candidates[sms_pep.modification_idx_];
        }
        else
        {
          mod_candidate = unmod_candidate;
        }


        // create theoretical spectrum
        PeakSpectrum theo_spectrum;

        // add peaks for b and y ions with charge 1
        spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, 1);

        // sort by mz
        theo_spectrum.sortByPosition();

        // const int& charge = exp_spectrum.getPrecursors()[0].getCharge();
        HyperScore::PSMDetail detail;
        const double& score = HyperScore::computeWithDetail(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);

        if (score == 0)
        { 
          continue; // no hit?
        }

        // add peptide hit
        AnnotatedHit_ ah;
        ah.sequence = std::move(mod_candidate);
        ah.score = score;
        double seq_length = (double)ah.sequence.size();
        ah.prefix_fraction = (double)detail.matched_b_ions/seq_length;
        ah.suffix_fraction = (double)detail.matched_y_ions/seq_length;
        ah.mean_error = detail.mean_error;

        annotated_hits[scan_index].push_back(std::move(ah)); 
      }
    }

    endProgress();

    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);

    startProgress(0, 1, "Post-processing PSMs...");
    PeptideSearchEngineAlgorithm::postProcessHits_(spectra, 
      annotated_hits, 
      protein_ids, 
      peptide_ids, 
      report_top_hits_,
      fixed_modifications, 
      variable_modifications, 
      modifications_max_variable_mods_per_peptide_,
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      precursor_mass_tolerance_,
      fragment_mass_tolerance_,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      in_db
      );
    endProgress();

    // add meta data on spectra file
    protein_ids[0].setPrimaryMSRunPath({in_mzML}, spectra);

    // reindex peptides to proteins
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string", "DECOY_");
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", enzyme_);
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
      {
        return ExitCodes::INPUT_FILE_EMPTY;       
      }
      else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
      {
        return ExitCodes::UNEXPECTED_RESULT;
      }
      else
      {
        return ExitCodes::UNKNOWN_ERROR;
      }
    }

    return ExitCodes::EXECUTION_OK;
  }

} // namespace OpenMS

