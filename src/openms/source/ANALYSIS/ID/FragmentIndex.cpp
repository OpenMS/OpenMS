//
// Created by trapho on 10/5/23.
//
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

  void FragmentIndex::addSpecialPeptide( OpenMS::AASequence& peptide, Size source_idx)

  {
    double temp_mono = peptide.getMonoWeight();
    fi_peptides_.push_back({exp_type_, AASequence(std::move(peptide)), source_idx,temp_mono});
  }

  void FragmentIndex::clear()
  {
    fi_fragments_.clear();
    fi_peptides_.clear();
  }

  void FragmentIndex::generate_peptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  { //TODO: Multithreading
    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(modifications_fixed_);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(modifications_variable_);
    size_t skipped_peptides = 0;


      ProteaseDigestion digestor;
      digestor.setEnzyme(digestion_enzyme_);
      digestor.setMissedCleavages(missed_cleavages_);


      size_t protein_idx = 0;
      for(FASTAFile::FASTAEntry protein: fasta_entries){
        vector<pair<size_t ,size_t >> digested_peptides;
        /// DIGEST (if bottom-up)
        if(exp_type_)
        {
          digestor.digestUnmodified(StringView(protein.sequence), digested_peptides, peptide_min_length_, peptide_max_length_);
        }else{
          digested_peptides.push_back(make_pair(0, protein.sequence.length())); // in case of top down the whole protein is the peptide
        }

        for(pair<size_t, size_t > digested_peptide : digested_peptides){

          //remove peptides containing unknown AA
          if (protein.sequence.substr(digested_peptide.first, digested_peptide.second).find('X') != string::npos){
            skipped_peptides++;
            continue;
          }

          /// MODIFY (if modifications are specified)
          AASequence unmod_peptide = AASequence::fromString(protein.sequence.substr(digested_peptide.first, digested_peptide.second));
          double unmodified_mz = unmod_peptide.getMZ(1); // TODO: What is getMonoWeight??

          if (!(modifications_fixed_.empty() && modifications_variable_.empty())){
            vector<AASequence> modified_peptides;
            AASequence mod_peptide = AASequence(unmod_peptide); //copy the peptide

            ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, mod_peptide);
            ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, mod_peptide, max_variable_mods_per_peptide_, modified_peptides);

            for(AASequence modified_peptide: modified_peptides){
              double modified_mz = modified_peptide.getMZ(1);
              if(peptide_min_mass_ > modified_mz && modified_mz > peptide_max_mass_) //exclude peptides that are not in the min-max window
                continue;

              fi_peptides_.push_back({exp_type_, modified_peptide, protein_idx, modified_mz});
            }
          }else{
            if(peptide_min_mass_ < unmodified_mz && unmodified_mz < peptide_max_mass_)
              fi_peptides_.push_back({exp_type_,unmod_peptide, protein_idx, unmodified_mz});
          }


        }

        protein_idx++;
      }
      if(skipped_peptides > 0)
        OPENMS_LOG_WARN << skipped_peptides << " peptides skipped due to unkown AA \n";
      //sort the peptide vector, critical for following steps
      sort(fi_peptides_.begin(), fi_peptides_.end(), [](const Peptide& a, const Peptide& b){return a.mass < b.mass;});
  }

  void FragmentIndex::build(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
  {
      TheoreticalSpectrumGenerator tsg;
      PeakSpectrum b_y_ions;
      std::vector<Fragment> all_frags;
      /// generate all Peptides
      generate_peptides(fasta_entries);

      size_t peptide_idx = 0;
      /// For each Peptides get all theoretical b and y ions // TODO: include other fragmentation methods
      for (Peptide pep: fi_peptides_){
        tsg.getSpectrum(b_y_ions, pep.sequence, 1, 1);
        for (Peak1D frag : b_y_ions){
          if (fragment_min_mz_ > frag.getMZ() && frag.getMZ() > fragment_max_mz_  ) continue;
          fi_fragments_.push_back({peptide_idx,frag.getMZ()});

        }
        peptide_idx ++;
        b_y_ions.clear(true);
      }
      /// 1.) First all Fragments are sorted by their own mass!
      sort(fi_fragments_.begin(), fi_fragments_.end(), [](const Fragment& a, const Fragment& b) {
        return a.fragment_mz < b.fragment_mz;
      });

      /// Calculate the bucket size
      bucketsize_ = sqrt(fi_fragments_.size()); //Todo: MSFragger uses a different approach, which might be better
      cout << "creating DB with bucket_size " << bucketsize_ << endl;

      /// 2.) next sort after precursor mass and save the min_mz of each bucket
      for (size_t i = 0; i < fi_fragments_.size(); i += bucketsize_){
        bucket_min_mz_.emplace_back(fi_fragments_[i].fragment_mz);

        auto bucket_start = fi_fragments_.begin()+i;
        auto bucket_end = (i + bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : bucket_start + bucketsize_;

        sort(bucket_start, bucket_end, [](const Fragment& a, const Fragment& b) {
          return a.peptide_idx < b.peptide_idx;
        });
      }
      is_build_ = true;

  }

  std::pair<size_t, size_t > FragmentIndex::getPeptideRange(double precursor_mass, std::pair<double, double> window)
  {
      double prec_tol = (precursor_mz_tolerance_unit_ == "DA") ? precursor_mz_tolerance_ : Math::ppmToMass(precursor_mz_tolerance_, precursor_mass);
      //set include to false, be we dont want the lower element
      return binary_search_slice<Peptide, double>(fi_peptides_, precursor_mass - prec_tol + window.first, precursor_mass + prec_tol + window.second, [](Peptide a){return a.mass;}, false);
  }


  vector<FragmentIndex::Hit> FragmentIndex::query(OpenMS::Peak1D peak, pair<size_t, size_t> peptide_idx_range, uint16_t peak_charge)
  {
      double adjusted_mass = peak.getMZ() * peak_charge;

      vector<FragmentIndex::Hit> hits;
      double frag_tol = (fragment_mz_tolerance_unit_ == "DA") ? fragment_mz_tolerance_ : Math::ppmToMass(fragment_mz_tolerance_, adjusted_mass);  //TODO??? apply ppm to mass * charge or not

      //include set to true, because the lowest element might contain the actuall value
      auto in_range_buckets = binary_search_slice<double,double>(bucket_min_mz_, adjusted_mass- frag_tol , adjusted_mass + frag_tol, [](double a){return a;},true);

      for (size_t j = in_range_buckets.first; j <= in_range_buckets.second; j++){

        auto slice_start = fi_fragments_.begin() + (j*bucketsize_);
        auto slice_end = ((j+1) * bucketsize_) > fi_fragments_.size() ? fi_fragments_.end() : fi_fragments_.begin()+((j+1)*bucketsize_);
        vector<Fragment> slice(slice_start, slice_end);
        //include set to false, bc we do not want lower idx
        auto in_range_fragments = binary_search_slice<Fragment, size_t >(slice, peptide_idx_range.first, peptide_idx_range.second,
                                               [](Fragment a){return a.peptide_idx;}, false);
       for(size_t candidate = in_range_fragments.first; candidate <= in_range_fragments.second; candidate++){
          if(slice[candidate].fragment_mz >= (adjusted_mass - frag_tol ) && slice[candidate].fragment_mz <= (adjusted_mass + frag_tol)){
            hits.push_back({slice[candidate].peptide_idx, slice[candidate]. fragment_mz});
            if(slice[candidate].peptide_idx < peptide_idx_range.first || slice[candidate].peptide_idx > peptide_idx_range.second)
              OPENMS_LOG_WARN << "idx out of range" << endl;
          }
       }
      }

      return hits;
  }


  template<class S, class B> std::pair<size_t, size_t> FragmentIndex::binary_search_slice(const std::vector<S>& slice, B low, B high, B (*access)(S), bool include)
  {
      pair<size_t, size_t> out;

      ///1.) find low index
      auto low_idx = slice.begin();
      auto high_idx = slice.end() - 1;
      auto mid_idx = low_idx + (high_idx - low_idx)/2;

      while((high_idx - low_idx) > 1){

        mid_idx = low_idx + (high_idx - low_idx)/2;
        if(low > access(*mid_idx))    // here is the only difference between low and high
          low_idx = mid_idx;
        else
          high_idx = mid_idx;
      }


      out.first = std::distance(slice.begin(),((access(*high_idx) >= low) || (access(*high_idx) == access(*low_idx))) ? low_idx : high_idx);
      //must be >= because of how exclude works!
      if(!include)  // we do not want to include elements that are lower, with this we also get such pairs (4,3) if the value is incorrect
                    // this automatically results later that we wont enter loops were we did not find any slices!!
        out.first++;


      ///2.) high index
      low_idx = slice.begin();
      high_idx =  slice.end() - 1;
      while((high_idx - low_idx) > 1){

        mid_idx = low_idx + (high_idx - low_idx)/2;
        if(high >= access(*mid_idx))
          low_idx = mid_idx;
        else
          high_idx = mid_idx;
      }
      out.second = std::distance(slice.begin(),(access(*high_idx) <= high) ? high_idx : low_idx);
      return out;

  }

  std::pair<size_t, size_t> FragmentIndex::binary_search_slice_double(const std::vector<double>& slice, double low, double high, bool include)
  {
      return binary_search_slice<double, double>(slice, low, high, [](double a){return a;}, include);
  }

  std::pair<size_t, size_t> FragmentIndex::binary_search_slice_mf(const std::vector<OpenMS::MultiFragment>& slice, size_t low, size_t high, size_t (*access) (OpenMS::MultiFragment), bool include)
  {
      return binary_search_slice<MultiFragment, size_t >(slice, low, high, access, include);
  }

  FragmentIndex::FragmentIndex() : DefaultParamHandler("FragmentIndex")
  {


    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    all_mods.push_back({});
    vector<string> tolerance_units{"DA", "PPM"};
    defaults_.setValue("experiment_type","false", "Bottom Up (true) or TopDown (false)");
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
    is_build_ = false;
    defaultsToParam_();
  }

  void FragmentIndex::updateMembers_()
  {
    exp_type_ = param_.getValue("experiment_type").toBool();
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
    precursor_mz_tolerance_unit_ = param_.getValue("precursor_mz_tolerance_unit").toString();
    fragment_mz_tolerance_unit_ = param_.getValue("fragment_mz_tolerance_unit").toString();

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_fixed"));
    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications_variable"));
    max_variable_mods_per_peptide_ = param_.getValue("max_variable_mods_per_peptide");
  }
  vector<AASequence> FragmentIndex::getFiPeptidesSequences() const
  {
    vector<AASequence> output;
    for(Peptide pep: fi_peptides_){
        output.push_back(AASequence(pep.sequence));
    }
    return output;
  }
  bool FragmentIndex::isBuild() const
  {
    return is_build_;
  }
  const vector<FragmentIndex::Peptide>& FragmentIndex::getFiPeptides() const
  {
    return fi_peptides_;
  }

}