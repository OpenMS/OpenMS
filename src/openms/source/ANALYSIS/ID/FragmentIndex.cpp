// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#include <functional>


using namespace std;



namespace OpenMS
{

#ifdef DEBUG_FRAGMENT_INDEX
    static void print_slice(const std::vector<FragmentIndex::Fragment>& slice, size_t low, size_t high)
    {
      cout << "Slice: ";
      for(size_t i = low; i <= high; i++)
      {
        cout << slice[i].fragment_mz_ << " ";
      }
      cout << endl;
    }


  void FragmentIndex::addSpecialPeptide( OpenMS::AASequence& peptide, Size source_idx)
  {
    float temp_mono = peptide.getMonoWeight();
    fi_peptides_.push_back({AASequence(std::move(peptide)), source_idx,temp_mono});
  }
#endif

  void FragmentIndex::clear()
  {
    fi_fragments_.clear();
    fi_peptides_.clear();
    bucket_min_mz_.clear();
    is_build_ = false;
  }


  // TODO: check if it makes sense to stream fasta from disc (we have some code for that... would safe memory but might have other drawbacks)
  void FragmentIndex::generatePeptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {

      fi_peptides_.reserve(fasta_entries.size() * 5 * modifications_variable_.size()); //TODO: Calculate the average cleavage site number for the most important model organisms
      ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
      ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
      size_t skipped_peptides = 0;

      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);

      std::cout << "Generating peptides..." << std::endl;
      
      vector<pair<size_t, size_t>> digested_peptides; // every thread gets it own copy that is only cleared, not destructed (prevents frequent reallocations)
      #pragma omp parallel for private(digested_peptides)
      for (SignedSize i = 0; i < fasta_entries.size(); ++i)
      {
        digested_peptides.clear();
        const FASTAFile::FASTAEntry& protein = fasta_entries[i];
        /// DIGEST (if bottom-up)
        digestor.digestUnmodified(StringView(protein.sequence), digested_peptides, peptide_min_length_, peptide_max_length_);

        for (const pair<size_t, size_t>& digested_peptide : digested_peptides)
        {
          //remove peptides containing unknown AA
          if (protein.sequence.substr(digested_peptide.first, digested_peptide.second).find('X') != string::npos)
          {
            #pragma omp atomic
            skipped_peptides++;
            continue;
          }

          /// MODIFY (if modifications are specified)
          AASequence unmod_peptide = AASequence::fromString(protein.sequence.substr(digested_peptide.first, digested_peptide.second));
          float unmodified_mz = unmod_peptide.getMZ(1);

          if (!(modifications_fixed_.empty() && modifications_variable_.empty()))
          {
            vector<AASequence> modified_peptides;
            AASequence mod_peptide = AASequence(unmod_peptide); // copy the peptide

            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, modified_peptides);

            UInt32 modification_idx = 0;
            for (const AASequence& modified_peptide : modified_peptides)
            {
              float modified_mz = modified_peptide.getMZ(1);
              if (modified_mz < peptide_min_mass_ || modified_mz > peptide_max_mass_) // exclude peptides that are not in the min-max window
              {
                continue;   //TODO: Integrate this check in ModPepGen from #6859
              }
              #pragma omp critical (FIIndex)
              {
                fi_peptides_.emplace_back(static_cast<UInt32>(i),
                                        modification_idx,
                                        digested_peptide,
                                        modified_mz);
                ++modification_idx;
              }
            }
          }
          else
          {
            if (peptide_min_mass_ <= unmodified_mz && unmodified_mz <= peptide_max_mass_)
            {
              #pragma omp critical (FIIndex)
              {
                fi_peptides_.emplace_back(static_cast<UInt32>(i), 0, digested_peptide, unmodified_mz);
              }
            }
          }
        }
      }
      if (skipped_peptides > 0)
      {
        OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unkown AA \n";
      }
      std::cout << "Sorting peptides..." << std::endl;        
      // sort the peptide vector, critical for following steps
      sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b)
           {
        return std::tie(a.precursor_mz_, a.protein_idx) < std::tie(b.precursor_mz_, b.protein_idx);
           });
      std::cout << "done." << std::endl;        
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      // reserve some memory for fi_fragments_: Each Peptide can approx give rise to up to #AA*2 fragments
      fi_fragments_.reserve(fi_peptides_.size() * 2 * peptide_min_length_); //TODO: Does this make senese?

      // get the spectrum generator and set the ion-types
      //TheoreticalSpectrumGenerator tsg;
      SimpleTSGXLMS tsg;

      auto tsg_params = tsg.getParameters();
      auto this_params = getParameters();
      tsg_params.setValue("add_a_ions", this_params.getValue("add_a_ions"));
      tsg_params.setValue("add_b_ions", this_params.getValue("add_b_ions"));
      tsg_params.setValue("add_c_ions", this_params.getValue("add_c_ions"));
      tsg_params.setValue("add_x_ions", this_params.getValue("add_x_ions"));
      tsg_params.setValue("add_y_ions", this_params.getValue("add_y_ions"));
      tsg_params.setValue("add_z_ions", this_params.getValue("add_z_ions"));
      tsg_params.setValue("add_first_prefix_ion", "true");
      tsg.setParameters(tsg_params);


      /// generate all Peptides
      generatePeptides(fasta_entries);

      /// Since we (the new) Peptide struct does not store the AASequence, we must reconstruct the modified ones
      /// therefore we need the modificationGenerators:
      ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
      ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);


      vector<AASequence> mod_peptides;
      std::vector<SimpleTSGXLMS::SimplePeak> b_y_ions;

      OPENMS_LOG_INFO << "Generating fragments..." << std::endl;

      #pragma omp parallel for private(mod_peptides, b_y_ions)
      for(size_t peptide_idx = 0; peptide_idx < fi_peptides_.size(); peptide_idx++)
      {
        const Peptide& pep = fi_peptides_[peptide_idx];
        mod_peptides.clear();
        b_y_ions.clear();
        AASequence unmod_peptide = AASequence::fromString(fasta_entries[pep.protein_idx].sequence.substr(pep.sequence_.first, pep.sequence_.second));

        if (!(modifications_fixed_.empty() && modifications_variable_.empty()))
        {
          AASequence mod_peptide = AASequence(unmod_peptide); // copy the peptide
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, mod_peptides);
          //tsg.getSpectrum(b_y_ions, mod_peptides[pep.modification_idx_], 1,1);
          tsg.getLinearIonSpectrum(b_y_ions, mod_peptides[pep.modification_idx_], mod_peptides[pep.modification_idx_].size(), 1, 1);
        }
        else
        {
         //tsg.getSpectrum(b_y_ions, unmod_peptide, 1,1);
          tsg.getLinearIonSpectrum(b_y_ions, unmod_peptide, unmod_peptide.size(), 1, 1);
        }

        for (const SimpleTSGXLMS::SimplePeak& frag : b_y_ions)
        {
          if (fragment_min_mz_ > frag.mz || frag.mz > fragment_max_mz_  ) continue;

          #pragma omp critical (CreateFragment)
          fi_fragments_.emplace_back(static_cast<UInt32>(peptide_idx),(float) frag.mz);
        }        
      }

      std::cout << "Sorting fragments..." << std::endl;

      /// 1.) First all Fragments are sorted by their own mass!
      sort(fi_fragments_.begin(), fi_fragments_.end(), [](const Fragment& a, const Fragment& b)
      {
        return std::tie(a.fragment_mz_, a.peptide_idx_) < std::tie(b.fragment_mz_, b.peptide_idx_);
      });

      /// Calculate the bucket size
      bucketsize_ = sqrt(fi_fragments_.size()); //Todo: MSFragger uses a different approach, which might be better
      OPENMS_LOG_INFO << "Creating DB with bucket_size " << bucketsize_ << endl;

      /// 2.) next sort after precursor mass and save the min_mz of each bucket
      #pragma omp parallel for
      for (size_t i = 0; i < fi_fragments_.size(); i += bucketsize_)
      {

        #pragma omp critical
        bucket_min_mz_.emplace_back(fi_fragments_[i].fragment_mz_);

        auto bucket_start = fi_fragments_.begin() + i;
        auto bucket_end = (i + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

//TODO: is this thread safe????
        sort(bucket_start, bucket_end, [](const Fragment& a, const Fragment& b) {
          return a.peptide_idx_ < b.peptide_idx_; // we don´t need a tie, because the idx are unique
        });
      }
      OPENMS_LOG_INFO << "Sorting by bucket min m/z:" << bucketsize_ << endl;
      //Resort in case the parallelization block above messed something up TODO: check if this can happen
      std::sort( bucket_min_mz_.begin(), bucket_min_mz_.end());
      is_build_ = true;
      OPENMS_LOG_INFO << "Fragment index built!" << endl;
  }

  std::pair<size_t, size_t > FragmentIndex::getPeptidesInPrecursorRange(float precursor_mass, const std::pair<float, float>& window)
  {
      float prec_tol = precursor_mz_tolerance_unit_ppm_ ? Math::ppmToMass(precursor_mz_tolerance_, precursor_mass) : precursor_mz_tolerance_ ;

      auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass - prec_tol + window.first, [](const Peptide& a, float b) { return a.precursor_mz_ < b;});
      auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass + prec_tol + window.second, [](float b, const Peptide& a) { return b < a.precursor_mz_;});
      return make_pair(std::distance(fi_peptides_.begin(), left_it), std::distance(fi_peptides_.begin(), right_it));
  }

  vector<FragmentIndex::Hit> FragmentIndex::query(const OpenMS::Peak1D& peak, const pair<size_t, size_t>& peptide_idx_range, uint16_t peak_charge)
  {
      float adjusted_mass = peak.getMZ() * (float)peak_charge -((peak_charge-1) * Constants::PROTON_MASS_U);

      float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass) : fragment_mz_tolerance_;

      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass - frag_tol);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass + frag_tol);

      if (left_it != bucket_min_mz_.begin()) --left_it;

      auto in_range_buckets = make_pair(std::distance(bucket_min_mz_.begin(), left_it), std::distance(bucket_min_mz_.begin(), right_it));

      vector<FragmentIndex::Hit> hits;
      hits.reserve(peptide_idx_range.second - peptide_idx_range.first);

      for (UInt32 j = in_range_buckets.first; j < in_range_buckets.second; j++)
      {
        auto slice_begin = fi_fragments_.begin() + (j*bucketsize_);
        auto slice_end = ((j+1) * bucketsize_) >= fi_fragments_.size() ? fi_fragments_.end() : (fi_fragments_.begin() + ((j+1) * bucketsize_)) ;

        auto left_iter = std::lower_bound(slice_begin, slice_end, peptide_idx_range.first, [](Fragment a, UInt32 b) { return a.peptide_idx_ < b;} );

        while (left_iter != slice_end) // sequential scan
        {
          if(left_iter->peptide_idx_ > peptide_idx_range.second) break;

          if ((adjusted_mass >= left_iter->fragment_mz_ - frag_tol ) && adjusted_mass <= (left_iter->fragment_mz_+ frag_tol))
          {
            hits.emplace_back(left_iter->peptide_idx_, left_iter->fragment_mz_);
            #ifdef DEBUG_FRAGMENT_INDEX
            if (left_iter->peptide_idx_ < peptide_idx_range.first || left_iter->peptide_idx_ > peptide_idx_range.second)
              OPENMS_LOG_WARN << "idx out of range" << endl;
            #endif
          }
          ++left_iter;
        }
      }

      return hits;
  }

  void FragmentIndex::queryPeak(OpenMS::FragmentIndex::SpectrumMatchesTopN& candidates,
                                const OpenMS::Peak1D& peak,
                                const std::pair<size_t, size_t>& candidates_range,
                                const int16_t isotope_error,
                                const uint16_t precursor_charge)
  {
      vector<Hit> query_hits;
      uint16_t actual_max = std::min(precursor_charge, max_fragment_charge_);
      for (uint16_t fragment_charge = 1; fragment_charge <= actual_max; fragment_charge++)
      {
        query_hits = query(peak, candidates_range, fragment_charge);

        for (const auto& hit : query_hits)
        {
          size_t idx = hit.peptide_idx - candidates_range.first;

          auto& source = candidates.hits_[idx];
          if (source.num_matched_ == 0)
          {
            ++candidates.scored_candidates_;
            source.precursor_charge_ = precursor_charge;
            source.peptide_idx_ = hit.peptide_idx;
            source.isotope_error_ = isotope_error;
          }
          ++source.num_matched_;
          ++candidates.matched_peaks_;
        }
      }
  }

  void FragmentIndex::trimHits(OpenMS::FragmentIndex::SpectrumMatchesTopN& init_hits) const
  {
      if (init_hits.hits_.size() > max_processed_hits_)
      {
        std::partial_sort(init_hits.hits_.begin(), init_hits.hits_.begin() + max_processed_hits_, init_hits.hits_.end(), [](const SpectrumMatch& a,const SpectrumMatch& b){
          if (a.num_matched_ != b.num_matched_)
          {
            return a.num_matched_ > b.num_matched_;
          }
          else
          {
            return std::tie(a.isotope_error_, a.precursor_charge_) < std::tie(b.isotope_error_, b.precursor_charge_);
          }
        });

        init_hits.hits_.resize(max_processed_hits_);
      }
      for (auto hit_iter = init_hits.hits_.rbegin(); hit_iter != init_hits.hits_.rend(); hit_iter++)
      {
        if (hit_iter->num_matched_ >= min_matched_peaks_)           // search for the first element that should be included
        {
          init_hits.hits_.resize(init_hits.hits_.size() - (distance(init_hits.hits_.rbegin(), hit_iter)));
          break;
        }
        if (hit_iter == init_hits.hits_.rend() -1) // we reached the last element without activating the previous if statement -> no hits at all
        {
          init_hits.hits_.resize(0);
        }
      }
  }

  void FragmentIndex::searchDifferentPrecursorRanges(const MSSpectrum& spectrum, float precursor_mass, SpectrumMatchesTopN& sms, uint16_t charge)
  {
      int16_t  min_isotope_error_applied;
      int16_t  max_isotope_error_applied;
      float precursor_window_upper_applied;
      float precursor_window_lower_applied;
      if (open_search)
      {
        min_isotope_error_applied = 0;
        max_isotope_error_applied = 0;
        precursor_window_upper_applied = open_precursor_window_upper_;
        precursor_window_lower_applied = open_precursor_window_lower_;
      }
      else
      {
        min_isotope_error_applied = min_isotope_error_;
        max_isotope_error_applied = max_isotope_error_;
        precursor_window_upper_applied = 0;
        precursor_window_lower_applied = 0;
      }
      for (int16_t isotope_error = min_isotope_error_applied; isotope_error <= max_isotope_error_applied; isotope_error++)
      {
        SpectrumMatchesTopN candidates_iso_error;
        float precursor_mass_isotope_error = precursor_mass + ((float)isotope_error * (float)Constants::C13C12_MASSDIFF_U);
        auto candidates_range = getPeptidesInPrecursorRange(precursor_mass_isotope_error, {precursor_window_lower_applied, precursor_window_upper_applied}); // for the simple search we do not apply any modification window!!
        candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first + 1);

        for (const Peak1D& peak : spectrum)
        {
          queryPeak(candidates_iso_error, peak, candidates_range, isotope_error, charge);
        }
        // take only top 50 hits
        trimHits(candidates_iso_error);
        sms += candidates_iso_error;
      }
      trimHits(sms);
  }

  void FragmentIndex::querySpectrum(const OpenMS::MSSpectrum& spectrum, OpenMS::FragmentIndex::SpectrumMatchesTopN& sms)
  {
      if (!isBuild())
      {
        OPENMS_LOG_WARN << "FragmentIndex not yet build \n";
        return;
      }

      if (spectrum.empty() || (spectrum.getMSLevel() != 2))
      {
        return;
      }

      const auto& precursor = spectrum.getPrecursors();
      if (precursor.size() != 1)
      {
        OPENMS_LOG_WARN << "Number of precursors is not equal 1 \n";
        return;
      }

      // two posible modes. Precursor has a charge or we test all possible charges
      vector<size_t> charges;
      //cout << "precursor charge = " << precursor[0].getCharge() << endl;
      if (precursor[0].getCharge())
      {
        //cout << "precursor charge found" << endl;
        charges.push_back(precursor[0].getCharge());
      }
      else
      {
        for (uint16_t i = min_precursor_charge_; i <= max_precursor_charge_; i++)
        {
          charges.push_back(i);
        }
      }
      // loop over all PRECURSOR-charges

      for (uint16_t charge : charges)
      {
        SpectrumMatchesTopN candidates_charge;
        float mz;
        mz = (float)precursor[0].getMZ() * charge - ((charge-1) * Constants::PROTON_MASS_U);
        searchDifferentPrecursorRanges(spectrum, mz, candidates_charge, charge);
        sms += candidates_charge;
      }
      trimHits(sms);
  }

  FragmentIndex::FragmentIndex() : DefaultParamHandler("FragmentIndex")
  {
    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", {"true","false"});

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", {"true","false"});

    defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", {"true","false"});

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", {"true","false"});

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", {"true","false"});

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", {"true","false"});

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    vector<string> tolerance_units{"Da", "ppm"}; // TODO: check if same string in other OpenMS files
    defaults_.setValue("enzyme", "Trypsin", "Enzyme for digestion");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));
    defaults_.setValue("peptide:missed_cleavages", 1, "Missed cleavages for digestion");
    defaults_.setValue("peptide:min_mass", 100, "Minimal peptide mass for database");
    defaults_.setValue("peptide:max_mass", 9000, "Maximal peptide mass for database"); //Todo: set unlimited option
    defaults_.setValue("peptide:min_size", 7, "Minimal peptide length for database");
    defaults_.setValue("peptide:max_size", 40, "Maximal peptide length for database");
    defaults_.setValue("fragment_min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment_max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("precursor:mass_tolerance", 10, "Tolerance for precursor-m/z in search");
    defaults_.setValue("fragment:mass_tolerance", 10, "Tolerance for fragment-m/z in search");
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of tolerance for precursor-m/z");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", tolerance_units);
    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of tolerance for fragment-m/z");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", tolerance_units);
    defaults_.setValue("modifications_fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications_fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications_variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications_variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("max_variable_mods_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    is_build_ = false; // TODO: remove this and build on construction

    //Search-related params

    defaults_.setValue("min_matched_peaks", 5, "Minimal number of matched ions to report a PSM");
    defaults_.setValue("min_isotope_error", -1, "Precursor isotope error");
    defaults_.setValue("max_isotope_error", 1, "precursor isotope error");
    defaults_.setValue("precursor:min_charge", 2, "min precursor charge");
    defaults_.setValue("precursor:max_charge", 5, "max precursor charge");
    defaults_.setValue("fragment:max_charge", 4, "max fragment charge");
    defaults_.setValue("max_processed_hits", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue("open_search", "false", "Open or standard search");
    defaults_.setValue("open_precursor_window_lower_", -100.0, "lower bound of the open precursor window");
    defaults_.setValue("open_precursor_window_upper_", 200.0, "upper bound of the open precursor window");

    defaultsToParam_();
}

  void FragmentIndex::updateMembers_()
  {
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    digestion_enzyme_ = param_.getValue("enzyme").toString();
    missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_min_mass_ = param_.getValue("peptide:min_mass");
    peptide_max_mass_ = param_.getValue("peptide:max_mass");
    peptide_min_length_ = param_.getValue("peptide:min_size");
    peptide_max_length_ = param_.getValue("peptide:max_size");
    fragment_min_mz_ = param_.getValue("fragment_min_mz");
    fragment_max_mz_ = param_.getValue("fragment_max_mz");

    precursor_mz_tolerance_ = param_.getValue("precursor:mass_tolerance");
    fragment_mz_tolerance_ = param_.getValue("fragment:mass_tolerance");
    precursor_mz_tolerance_unit_ppm_ = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment:mass_tolerance_unit").toString() == "ppm";

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_variable"));
    max_variable_mods_per_peptide_ = param_.getValue("max_variable_mods_per_peptide");

    min_matched_peaks_ = param_.getValue("min_matched_peaks");
    min_isotope_error_ = param_.getValue("min_isotope_error");
    max_isotope_error_ = param_.getValue("max_isotope_error");
    min_precursor_charge_ = param_.getValue("precursor:min_charge");
    max_precursor_charge_ = param_.getValue("precursor:max_charge");
    max_fragment_charge_ = param_.getValue("fragment:max_charge");
    max_processed_hits_ = param_.getValue("max_processed_hits");
    open_search = param_.getValue("open_search").toBool();
    open_precursor_window_lower_ = param_.getValue("open_precursor_window_lower_");
    open_precursor_window_upper_ = param_.getValue("open_precursor_window_upper_");
  }
 
  bool FragmentIndex::isBuild() const
  {
    return is_build_;
  }

  const vector<FragmentIndex::Peptide>& FragmentIndex::getPeptides() const
  {
    return fi_peptides_;
  }

}
