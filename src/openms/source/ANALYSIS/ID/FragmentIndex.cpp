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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
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
        cout << slice[i].fragment_mz << " ";
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
    is_build_ = false;
  }

  FragmentIndex::IonTypes FragmentIndex::getIonTypes() const
  {
    return ion_types_;
  }

  // TODO: check if it makes sense to stream fasta from disc (we have some code for that... would safe memory but might have other drawbacks)
  void FragmentIndex::generate_peptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  { 
      ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
      ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
      size_t skipped_peptides = 0;

      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);

      size_t protein_idx = 0;
      
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
            skipped_peptides++;
            continue;
          }

          /// MODIFY (if modifications are specified)
          AASequence unmod_peptide = AASequence::fromString(protein.sequence.substr(digested_peptide.first, digested_peptide.second));
          float unmodified_mz = unmod_peptide.getMZ(1); // TODO: What is getMonoWeight??

          if (!(modifications_fixed_.empty() && modifications_variable_.empty()))
          {
            vector<AASequence> modified_peptides;
            AASequence mod_peptide = AASequence(unmod_peptide); // copy the peptide

            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, modified_peptides);

            for (const AASequence& modified_peptide : modified_peptides)
            {
              float modified_mz = modified_peptide.getMZ(1);
              if (modified_mz < peptide_min_mass_ || modified_mz > peptide_max_mass_) // exclude peptides that are not in the min-max window
              {
                continue;
              }
                
              fi_peptides_.push_back({modified_peptide, protein_idx, modified_mz});
            }
          }
          else
          {
            if (peptide_min_mass_ < unmodified_mz && unmodified_mz < peptide_max_mass_) // TODO: <= instad of <?
            {
              fi_peptides_.push_back({unmod_peptide, protein_idx, unmodified_mz});
            }              
          }
        }

        protein_idx++;
      }
      if (skipped_peptides > 0)
      {
        OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unkown AA \n";
      }        
      // sort the peptide vector, critical for following steps
      sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b){ return a.precursor_mz < b.precursor_mz; });
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      // get the spectrum generator and set the ion-types
      TheoreticalSpectrumGenerator tsg;
      auto tsg_params = tsg.getParameters();
      tsg_params.setValue("add_a_ions", add_a_ions_);
      tsg_params.setValue("add_b_ions", add_b_ions_);
      tsg_params.setValue("add_c_ions", add_c_ions_);
      tsg_params.setValue("add_x_ions", add_x_ions_);
      tsg_params.setValue("add_y_ions", add_y_ions_);
      tsg_params.setValue("add_z_ions", add_z_ions_);
      tsg.setParameters(tsg_params); // TODO: this was missing before!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PeakSpectrum b_y_ions;

      /// generate all Peptides
      generate_peptides(fasta_entries);

      size_t peptide_idx = 0;

      /// For each Peptides get all theoretical b and y ions // TODO: include other fragmentation methods
      for (const Peptide& pep: fi_peptides_)
      {
        tsg.getSpectrum(b_y_ions, pep.sequence, 1, 1);
        for (Peak1D frag : b_y_ions)
        {
          if (fragment_min_mz_ > frag.getMZ() || frag.getMZ() > fragment_max_mz_  ) continue;
          fi_fragments_.push_back({peptide_idx,frag.getMZ()});
        }
        peptide_idx ++;
        b_y_ions.clear(true);
      }
      /// 1.) First all Fragments are sorted by their own mass!
      sort(fi_fragments_.begin(), fi_fragments_.end(), [](const Fragment& a, const Fragment& b) 
      {
        return a.fragment_mz < b.fragment_mz;
      });

      /// Calculate the bucket size
      bucketsize_ = sqrt(fi_fragments_.size()); //Todo: MSFragger uses a different approach, which might be better
      cout << "creating DB with bucket_size " << bucketsize_ << endl;

      /// 2.) next sort after precursor mass and save the min_mz of each bucket
      for (size_t i = 0; i < fi_fragments_.size(); i += bucketsize_)
      {
        bucket_min_mz_.emplace_back(fi_fragments_[i].fragment_mz);

        auto bucket_start = fi_fragments_.begin() + i;
        auto bucket_end = (i + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

        sort(bucket_start, bucket_end, [](const Fragment& a, const Fragment& b) 
        {
          return a.peptide_idx < b.peptide_idx;
        });
      }
      is_build_ = true;
  }

  std::pair<size_t, size_t > FragmentIndex::getPeptidesInPrecursorRange(float precursor_mass, std::pair<float, float> window)
  {
      float prec_tol = precursor_mz_tolerance_unit_ppm_ ? Math::ppmToMass(precursor_mz_tolerance_, precursor_mass) : precursor_mz_tolerance_ ;  //TODO??? apply ppm to mass * charge or not

      auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass - prec_tol + window.first, [](const Peptide& a, float b) { return a.precursor_mz < b;});
      auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass + prec_tol + window.second, [](float b, const Peptide& a) { return b < a.precursor_mz;});
      return make_pair(std::distance(fi_peptides_.begin(), left_it), std::distance(fi_peptides_.begin(), right_it));
  }

  vector<FragmentIndex::Hit> FragmentIndex::query(OpenMS::Peak1D peak, pair<size_t, size_t> peptide_idx_range, uint16_t peak_charge)
  {
      float adjusted_mass = peak.getMZ() * peak_charge;

      float frag_tol = fragment_mz_tolerance_unit_ppm_ ? Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass) : fragment_mz_tolerance_;  //TODO??? apply ppm to mass * charge or not

      auto left_it = std::lower_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass - frag_tol);
      auto right_it = std::upper_bound(bucket_min_mz_.begin(), bucket_min_mz_.end(), adjusted_mass + frag_tol);      
      if (left_it != bucket_min_mz_.begin()) --left_it;

      auto in_range_buckets = make_pair(std::distance(bucket_min_mz_.begin(), left_it), std::distance(bucket_min_mz_.begin(), right_it));

      vector<FragmentIndex::Hit> hits;
      hits.reserve(64);

      for (size_t j = in_range_buckets.first; j <= in_range_buckets.second; j++)
      {
        auto slice_begin = fi_fragments_.begin() + (j*bucketsize_);
        auto slice_end = ((j+1) * bucketsize_) >= fi_fragments_.size() ? fi_fragments_.end() : fi_fragments_.begin() + ((j+1) * bucketsize_) + 1;

        auto left_it = std::lower_bound(slice_begin, slice_end, peptide_idx_range.first, [](Fragment a, UInt32 b) { return a.peptide_idx < b;} ); 

        while (left_it != slice_end && left_it->peptide_idx < peptide_idx_range.second) // sequential scan
        {
          if (left_it->fragment_mz >= (adjusted_mass - frag_tol ) && left_it->fragment_mz <= (adjusted_mass + frag_tol))
          {
            hits.emplace_back(left_it->peptide_idx, left_it->fragment_mz);
            #ifdef DEBUG_FRAGMENT_INDEX
            if (left_it->peptide_idx < peptide_idx_range.first || left_it->peptide_idx > peptide_idx_range.second)
              OPENMS_LOG_WARN << "idx out of range" << endl;
            #endif
          }
          ++left_it;
        }

      }

      return hits;
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

    vector<string> tolerance_units{"DA", "PPM"}; // TODO: check if same string in other OpenMS files
    defaults_.setValue("digestor_enzyme", "Trypsin", "Enzyme for digestion");
    defaults_.setValidStrings("digestor_enzyme", ListUtils::create<std::string>(all_enzymes));
    defaults_.setValue("missed_cleavages", 0, "Missed cleavages for digestion");
    defaults_.setValue("peptide_min_mass", 1, "Minimal peptide mass for database");
    defaults_.setValue("peptide_max_mass", 500000, "Maximal peptide mass for database"); //Todo: set unlimited option
    defaults_.setValue("peptide_min_length", 10, "Minimal peptide length for database");
    defaults_.setValue("peptide_max_length", 5000000, "Maximal peptide length for database");
    defaults_.setValue("fragment_min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment_max_mz", 500000, "Maximal fragment mz for database");
    defaults_.setValue("precursor_mz_tolerance", 2.0, "Tolerance for precursor-m/z in search");
    defaults_.setValue("fragment_mz_tolerance", 0.05, "Tolerance for fragment-m/z in search");
    defaults_.setValue("precursor_mz_tolerance_unit", "DA", "Unit of tolerance for precursor-m/z");
    defaults_.setValidStrings("precursor_mz_tolerance_unit", tolerance_units);
    defaults_.setValue("fragment_mz_tolerance_unit", "DA", "Unit of tolerance for fragment-m/z");
    defaults_.setValidStrings("fragment_mz_tolerance_unit", tolerance_units);
    defaults_.setValue("modifications_fixed", std::vector<std::string>{}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications_fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications_variable", std::vector<std::string>{}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications_variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("max_variable_mods_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    is_build_ = false; // TODO: remove this and build on construction
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
    digestion_enzyme_ = param_.getValue("digestor_enzyme").toString();
    missed_cleavages_ = param_.getValue("missed_cleavages");
    peptide_min_mass_ = param_.getValue("peptide_min_mass");
    peptide_max_mass_ = param_.getValue("peptide_max_mass");
    peptide_min_length_ = param_.getValue("peptide_min_length");
    peptide_max_length_ = param_.getValue("peptide_max_length");
    fragment_min_mz_ = param_.getValue("fragment_min_mz");
    fragment_max_mz_ = param_.getValue("fragment_max_mz");

    precursor_mz_tolerance_ = param_.getValue("precursor_mz_tolerance");
    fragment_mz_tolerance_ = param_.getValue("fragment_mz_tolerance");
    precursor_mz_tolerance_unit_ppm_ = param_.getValue("precursor_mz_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment_mz_tolerance_unit").toString() == "ppm";

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_variable"));
    max_variable_mods_per_peptide_ = param_.getValue("max_variable_mods_per_peptide");
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
