// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

// #define MIN_DOUBLE_MZ 900.0

// #define SELECT_PIVOT_DEBUG

using namespace std;

namespace OpenMS
{
  CompNovoIdentificationBase::CompNovoIdentificationBase() :
    DefaultParamHandler("CompNovoIdentificationBase"),
    max_number_aa_per_decomp_(0),
    tryptic_only_(true),
    fragment_mass_tolerance_(0),
    max_number_pivot_(0),
    decomp_weights_precision_(0),
    max_mz_(2000.0),
    min_mz_(200.0),
    max_decomp_weight_(450.0),
    max_subscore_number_(30),
    max_isotope_(3)
  {
    defaults_.setValue("max_number_aa_per_decomp", 4, "maximal amino acid frequency per decomposition", ListUtils::create<String>("advanced"));
    defaults_.setValue("tryptic_only", "true", "if set to true only tryptic peptides are reported");
    defaults_.setValue("precursor_mass_tolerance", 1.5, "precursor mass tolerance");
    defaults_.setValue("fragment_mass_tolerance", 0.3, "fragment mass tolerance");
    defaults_.setValue("max_number_pivot", 9, "maximal number of pivot ions to be used", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_subscore_number", 40, "maximal number of solutions of a subsegment that are kept", ListUtils::create<String>("advanced"));
    defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!", ListUtils::create<String>("advanced"));
    defaults_.setValue("double_charged_iso_threshold", 0.6, "minimal isotope intensity correlation of doubly charged ions to be used to score the single scored ions", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_mz", 2000.0, "maximal m/z value used to calculate isotope distributions");
    defaults_.setValue("min_mz", 200.0, "minimal m/z value used to calculate the isotope distributions");
    defaults_.setValue("max_isotope_to_score", 3, "max isotope peak to be considered in the scoring", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_decomp_weight", 450.0, "maximal m/z difference used to calculate the decompositions", ListUtils::create<String>("advanced"));
    defaults_.setValue("max_isotope", 3, "max isotope used in the theoretical spectra to score", ListUtils::create<String>("advanced"));
    defaults_.setValue("missed_cleavages", 1, "maximal number of missed cleavages allowed per peptide");
    defaults_.setValue("number_of_hits", 100, "maximal number of hits which are reported per spectrum");
    defaults_.setValue("estimate_precursor_mz", "true", "If set to true, the precursor charge will be estimated, e.g. from the precursor peaks of the ETD spectrum.\n"
                                                        "The input is believed otherwise.");
    defaults_.setValidStrings("estimate_precursor_mz", ListUtils::create<String>("true,false"));
    defaults_.setValue("number_of_prescoring_hits", 250, "how many sequences are kept after first rough scoring for better scoring", ListUtils::create<String>("advanced"));

    // set all known modifications as restriction
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("fixed_modifications", ListUtils::create<String>(""), "fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
    defaults_.setValidStrings("fixed_modifications", all_mods);

    defaults_.setValue("variable_modifications", ListUtils::create<String>(""), "variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
    defaults_.setValidStrings("variable_modifications", all_mods);

    defaults_.setValue("residue_set", "Natural19WithoutI", "The predefined amino acid set that should be used, see doc of ResidueDB for possible residue sets", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  CompNovoIdentificationBase::CompNovoIdentificationBase(const CompNovoIdentificationBase & rhs) :
    DefaultParamHandler(rhs)
  {
    updateMembers_();
  }

  CompNovoIdentificationBase & CompNovoIdentificationBase::operator=(const CompNovoIdentificationBase & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
      updateMembers_();
      // TODO
    }
    return *this;
  }

  CompNovoIdentificationBase::~CompNovoIdentificationBase() = default;

  void CompNovoIdentificationBase::getCIDSpectrumLight_(PeakSpectrum & spec, const String & sequence, double prefix, double suffix)
  {
    static double h2o_mass = EmpiricalFormula("H2O").getMonoWeight();
    Peak1D p;
    double b_pos(0.0 + prefix);
    double y_pos(h2o_mass + suffix);
    for (Size i = 0; i != sequence.size() - 1; ++i)
    {
      char aa(sequence[i]);
      b_pos += aa_to_weight_[aa];

      char aa2(sequence[sequence.size() - i - 1]);
      y_pos += aa_to_weight_[aa2];

      if (b_pos > min_mz_ && b_pos < max_mz_)
      {
        p.setPosition(b_pos + Constants::PROTON_MASS_U);
        p.setIntensity(1.0f);
        spec.push_back(p);
      }

      if (y_pos > min_mz_ && y_pos < max_mz_)
      {
        p.setPosition(y_pos + Constants::PROTON_MASS_U);
        p.setIntensity(1.0f);
        spec.push_back(p);
      }
    }

    spec.sortByPosition();
  }

  void CompNovoIdentificationBase::getCIDSpectrum_(PeakSpectrum & spec, const String & sequence, Size charge, double prefix, double suffix)
  {
    if (isotope_distributions_.empty())
    {
      initIsotopeDistributions_();
    }
    static double h2o_mass = EmpiricalFormula("H2O").getMonoWeight();
    static double nh3_mass = EmpiricalFormula("NH3").getMonoWeight();
    static double co_mass = EmpiricalFormula("CO").getMonoWeight();
    Peak1D p;
    double b_pos(0 + prefix);
    double y_pos(h2o_mass + suffix);
    bool b_H2O_loss(false), b_NH3_loss(false), y_NH3_loss(false);

    for (Size i = 0; i != sequence.size() - 1; ++i)
    {
      char aa(sequence[i]);
      b_pos += aa_to_weight_[aa];

      char aa2(sequence[sequence.size() - i - 1]);
      y_pos += aa_to_weight_[aa2];
      for (Size z = 1; z <= charge && z < 3; ++z)
      {
        // b-ions
        if (b_pos >= min_mz_ && b_pos <= max_mz_)
        {
          for (Size j = 0; j != max_isotope_; ++j)
          {
            if (z == 1 /*|| b_pos > MIN_DOUBLE_MZ*/)
            {
              p.setPosition((b_pos + (double)z * Constants::PROTON_MASS_U + (double)j + Constants::NEUTRON_MASS_U) / (double)z);
              p.setIntensity(isotope_distributions_[(Size)b_pos][j] * 0.8 / (z * z));
              spec.push_back(p);
            }
          }
        }

        // b-ion losses
        if (b_pos - h2o_mass > min_mz_ && b_pos - h2o_mass < max_mz_)
        {
          if (b_H2O_loss || aa == 'S' || aa == 'T' || aa == 'E' || aa == 'D')
          {
            b_H2O_loss = true;
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - h2o_mass) / z);
            p.setIntensity(0.02 / (double)(z * z));
            if (z == 1 /* || b_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
          if (b_NH3_loss || aa == 'Q' || aa == 'N' || aa == 'R' || aa == 'K')
          {
            b_NH3_loss = true;
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - nh3_mass) / z);
            p.setIntensity(0.02 / (double)(z * z));

            if (z == 1 /* || b_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
        }

        // a-ions only for charge 1
        if (z == 1)
        {
          if (b_pos - co_mass > min_mz_ && b_pos - co_mass < max_mz_)
          {
            // a-ions
            p.setPosition((b_pos + z * Constants::PROTON_MASS_U - co_mass) / (double)z);
            p.setIntensity(0.1f);
            spec.push_back(p);
          }
        }

        if (y_pos > min_mz_ && y_pos < max_mz_)
        {
          // y-ions
          for (Size j = 0; j != max_isotope_; ++j)
          {
            if (z == 1 /* || y_pos > MIN_DOUBLE_MZ*/)
            {
              p.setPosition((y_pos + (double)z * Constants::PROTON_MASS_U + (double)j * Constants::NEUTRON_MASS_U) / (double)z);
              p.setIntensity(isotope_distributions_[(Size)y_pos][j] / (double) (z * z));
              spec.push_back(p);
            }
          }

          // H2O loss
          p.setPosition((y_pos + z * Constants::PROTON_MASS_U - h2o_mass) / (double)z);
          p.setIntensity(0.1 / (double)(z * z));
          if (aa2 == 'Q')           // pyroglutamic acid formation
          {
            p.setIntensity(0.5f);
          }
          if (z == 1 /* || y_pos > MIN_DOUBLE_MZ*/)
          {
            spec.push_back(p);
          }

          // NH3 loss
          if (y_NH3_loss || aa2 == 'Q' || aa2 == 'N' || aa2 == 'R' || aa2 == 'K')
          {
            y_NH3_loss = true;
            p.setPosition((y_pos + z * Constants::PROTON_MASS_U - nh3_mass) / (double)z);
            p.setIntensity(0.1 / (double)(z * z));

            if (z == 1 /*|| y_pos > MIN_DOUBLE_MZ*/)
            {
              spec.push_back(p);
            }
          }
        }
      }
    }

    // if Q1 abundant loss of water -> pyroglutamic acid formation

    //if (sequence[0] == 'Q' && prefix == 0 && suffix == 0)
    //{
      /*
      for (PeakSpectrum::Iterator it = spec.begin(); it != spec.end(); ++it)
      {
          it->setIntensity(it->getIntensity() * 0.5);
      }*/

      /*
      for (Size j = 0; j != max_isotope; ++j)
      {
  p.setPosition((precursor_weight + charge - 1 + j)/(double)charge);
  p.setIntensity(isotope_distributions_[(Int)p.getPosition()[0]][j] * 0.1);
  spec.push_back(p);
      }
      */
    //}
    spec.sortByPosition();
  }

  void CompNovoIdentificationBase::filterPermuts_(set<String> & permut)
  {
    set<String> tmp;
    for (set<String>::const_iterator it = permut.begin(); it != permut.end(); ++it)
    {
      if (tryptic_only_)
      {
        if ((*it)[it->size() - 1] == 'K' || (*it)[it->size() - 1] == 'R')
        {
          tmp.insert(*it);
        }
      }
      else
      {
        tmp.insert(*it);
      }
    }
    permut = tmp;
  }

  Size CompNovoIdentificationBase::countMissedCleavagesTryptic_(const String & peptide) const
  {
    Size missed_cleavages(0);

    if (peptide.size() < 2)
    {
      return 0;
    }
    for (Size i = 0; i != peptide.size() - 1; ++i)
    {
      if ((peptide[i] == 'R' || peptide[i] == 'K') && peptide[i + 1] != 'P')
      {
        ++missed_cleavages;
      }
    }

    return missed_cleavages;
  }

  void CompNovoIdentificationBase::permute_(const String& prefix, String s, set<String> & permutations)
  {
    if (s.size() <= 1)
    {
      permutations.insert(prefix + s);
    }
    else
    {
      for (String::Iterator p = s.begin(); p < s.end(); p++)
      {
        char c = *p;
        p = s.erase(p);
        permute_(prefix + c, s, permutations);
        s.insert(p, c);
      }
    }
  }

  void CompNovoIdentificationBase::getDecompositions_(vector<MassDecomposition> & decomps, double mass, bool no_caching)
  {
    //static Map<double, vector<MassDecomposition> > decomp_cache;
    if (!no_caching)
    {
      if (decomp_cache_.has(mass))
      {
        decomps = decomp_cache_[mass];
        return;
      }
    }

    mass_decomp_algorithm_.getDecompositions(decomps, mass);
    filterDecomps_(decomps);

    if (!no_caching)
    {
      decomp_cache_[mass]  = decomps;
    }
  }

  void CompNovoIdentificationBase::selectPivotIons_(vector<Size> & pivots, Size left, Size right, Map<double, CompNovoIonScoringBase::IonScore> & ion_scores, const PeakSpectrum & CID_spec, double precursor_weight, bool full_range)
  {
#ifdef SELECT_PIVOT_DEBUG
    cerr << "void selectPivotIons(pivots[" << pivots.size() << "], " << left << "[" << CID_spec[left].getPosition()[0] << "]" << ", " << right << "[" << CID_spec[right].getPosition()[0]  << "])" << endl;
#endif

    Size max_number_pivot(param_.getValue("max_number_pivot"));

    // TODO better heuristic, MAX_PIVOT dynamic from range
    if (right - left > 1)
    {
      right -= 1;
      left += 1;
      if (right - left < 1 || CID_spec[right].getPosition()[0] - CID_spec[left].getPosition()[0] < 57.0 - fragment_mass_tolerance_)
      {
        return;
      }
      // use more narrow window
      // diff between border and new pivot should be at least 57 - fragment_mass_tolerance (smallest aa)

      Size new_right(right), new_left(left);
      for (Size i = left - 1; i < right && CID_spec[i].getPosition()[0] - CID_spec[left - 1].getPosition()[0] < 57.0 - fragment_mass_tolerance_; ++i)
      {
        new_left = i;
      }

      for (Size i = right + 1; i > new_left &&
           CID_spec[right + 1].getPosition()[0] - CID_spec[i].getPosition()[0] < 57.0 - fragment_mass_tolerance_;
           --i)
      {
        new_right = i;
      }
#ifdef SELECT_PIVOT_DEBUG
      cerr << "new_left=" << new_left << "(" << CID_spec[new_left].getPosition()[0] << "), new_right=" << new_right << "(" << CID_spec[new_right].getPosition()[0] << ")" << endl;
#endif
      left = new_left;
      right = new_right;


      if (right - left <= 1)
      {
        return;
      }


      Size old_num_used(0);
      set<Size> used_pos;
      for (Size p = 0; p != min(right - left - 1, max_number_pivot); ++p)
      {
        double max(0);
        Size max_pos(0);

        bool found_pivot(false);
        for (Size i = left + 1; i < right; ++i)
        {
          double score = ion_scores[CID_spec[i].getPosition()[0]].score;
          double position = CID_spec[i].getPosition()[0];
#ifdef SELECT_PIVOT_DEBUG
          cerr << position << " " << precursor_weight << " " << full_range << " " << score;
#endif
          if (score >= max && used_pos.find(i) == used_pos.end())
          {
            // now check if a very similar ion is already selected +/- 3Da
            //bool has_similar(false);
            /*
for (set<Size>::const_iterator it = used_pos.begin(); it != used_pos.end(); ++it)
{
    if (fabs(CID_spec[*it].getPosition()[0] - CID_spec[i].getPosition()[0]) < 1.5)
    {
    has_similar = true;
    }
}*/

            // TODO this rule should be toggable
            if (!(full_range && (position < precursor_weight / 4.0 || position > precursor_weight / 4.0 * 3.0)))
            {
#ifdef SELECT_PIVOT_DEBUG
              cerr << " max score greater";
#endif
              max = score;
              max_pos = i;
              found_pivot = true;
            }
          }
#ifdef SELECT_PIVOT_DEBUG
          cerr << endl;
#endif
        }

        used_pos.insert(max_pos);

        // no pivot ion was added
        if (!found_pivot || (old_num_used == used_pos.size() && old_num_used != 0))
        {
          return;
        }
        else
        {
          old_num_used = used_pos.size();
        }

        pivots.push_back(max_pos);
        max = 0;
      }
    }
  }

  // s1 should be the original spectrum
  double CompNovoIdentificationBase::compareSpectra_(const PeakSpectrum & s1, const PeakSpectrum & s2)
  {
    double score(0.0);

    PeakSpectrum::ConstIterator it1 = s1.begin();
    PeakSpectrum::ConstIterator it2 = s2.begin();

    Size num_matches(0);
    while (it1 != s1.end() && it2 != s2.end())
    {
      double pos1(it1->getPosition()[0]), pos2(it2->getPosition()[0]);
      if (fabs(pos1 - pos2) < fragment_mass_tolerance_)
      {
        score += it1->getIntensity();
        ++num_matches;
      }

      if (pos1 <= pos2)
      {
        ++it1;
      }
      else
      {
        ++it2;
      }
    }

    if (num_matches == 0)
    {
      return 0;
    }

    score /= sqrt((double)num_matches);

    return score;
  }

  bool Internal::PermutScoreComparator(const CompNovoIdentificationBase::Permut & p1, const CompNovoIdentificationBase::Permut & p2)
  {
    return p1.getScore() > p2.getScore();
  }

  void CompNovoIdentificationBase::windowMower_(PeakSpectrum & spec, double windowsize, Size no_peaks)
  {
    PeakSpectrum copy(spec);
    vector<Peak1D> to_be_deleted;
    for (Size i = 0; i < spec.size(); ++i)
    {
      PeakSpectrum sub_spec;
      bool end(false);
      for (Size j = i; spec[j].getPosition()[0] - spec[i].getPosition()[0] < windowsize; )
      {
        sub_spec.push_back(spec[j]);
        if (++j == spec.size())
        {
          end = true;
          break;
        }
      }

      sub_spec.sortByIntensity(true);

      for (Size k = no_peaks; k < sub_spec.size(); ++k)
      {
        Peak1D p(sub_spec[k]);
        to_be_deleted.push_back(p);
      }

      if (end)
      {
        break;
      }
    }

    spec.clear(false);
    for (PeakSpectrum::ConstIterator it = copy.begin(); it != copy.end(); ++it)
    {
      if (find(to_be_deleted.begin(), to_be_deleted.end(), *it) == to_be_deleted.end())
      {
        spec.push_back(*it);
      }
    }

    spec.sortByPosition();

  }

  void CompNovoIdentificationBase::filterDecomps_(vector<MassDecomposition> & decomps)
  {
    Size max_number_aa_per_decomp(param_.getValue("max_number_aa_per_decomp"));
    vector<MassDecomposition> tmp;
    for (vector<MassDecomposition>::const_iterator it = decomps.begin(); it != decomps.end(); ++it)
    {
      if (it->getNumberOfMaxAA() <= max_number_aa_per_decomp)
      {
        tmp.push_back(*it);
      }
    }
    decomps = tmp;
  }

  void CompNovoIdentificationBase::initIsotopeDistributions_()
  {
    CoarseIsotopePatternGenerator solver(max_isotope_);
    for (Size i = 1; i <= max_mz_ * 2; ++i)
    {
      auto iso_dist = solver.estimateFromPeptideWeight((double)i);
      iso_dist.renormalize();
      vector<double> iso(max_isotope_, 0.0);

      for (Size j = 0; j != iso_dist.size(); ++j)
      {
        iso[j] = iso_dist.getContainer()[j].getIntensity();
      }
      isotope_distributions_[i] = iso;
    }
  }

  AASequence CompNovoIdentificationBase::getModifiedAASequence_(const String & sequence)
  {
    AASequence seq;
    for (String::ConstIterator it = sequence.begin(); it != sequence.end(); ++it)
    {
      if (name_to_residue_.has(*it))
      {
        seq += name_to_residue_[*it];
      }
      else
      {
        seq += AASequence::fromString(*it);
      }
    }

    return seq;
  }

  String CompNovoIdentificationBase::getModifiedStringFromAASequence_(const AASequence & sequence)
  {
    String seq;
    for (AASequence::ConstIterator it = sequence.begin(); it != sequence.end(); ++it)
    {
      if (residue_to_name_.has(&*it))
      {
        seq += residue_to_name_[&*it];
      }
      else
      {
        seq += it->getOneLetterCode();
      }
    }
    return seq;
  }

  void CompNovoIdentificationBase::updateMembers_()
  {
    // init residue mass table
    String residue_set(param_.getValue("residue_set"));

    set<const Residue *> residues = ResidueDB::getInstance()->getResidues(residue_set);
    for (set<const Residue *>::const_iterator it = residues.begin(); it != residues.end(); ++it)
    {
      aa_to_weight_[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
    }

    max_number_aa_per_decomp_ = param_.getValue("max_number_aa_per_decomp");
    tryptic_only_ = param_.getValue("tryptic_only").toBool();
    fragment_mass_tolerance_ = (double)param_.getValue("fragment_mass_tolerance");
    max_number_pivot_ = param_.getValue("max_number_pivot");
    decomp_weights_precision_ = (double)param_.getValue("decomp_weights_precision");
    min_mz_ = (double)param_.getValue("min_mz");
    max_mz_ = (double)param_.getValue("max_mz");
    max_decomp_weight_ = (double)param_.getValue("max_decomp_weight");
    max_subscore_number_ = param_.getValue("max_subscore_number");
    max_isotope_ = param_.getValue("max_isotope");

    name_to_residue_.clear();
    residue_to_name_.clear();

    // now handle the modifications
    ModificationDefinitionsSet mod_set(param_.getValue("fixed_modifications"), param_.getValue("variable_modifications"));
    set<ModificationDefinition> fixed_mods = mod_set.getFixedModifications();

    for (set<ModificationDefinition>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
    {
      ResidueModification mod = it->getModification();
      char aa = ' ';
      if (mod.getOrigin() == 'X')
      {
        cerr << "Warning: cannot handle modification " << mod.getName() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        aa = mod.getOrigin();
      }

      if (mod.getMonoMass() != 0)
      {
        aa_to_weight_[aa] = mod.getMonoMass();
      }
      else
      {
        if (mod.getDiffMonoMass() != 0)
        {
          aa_to_weight_[aa] += mod.getDiffMonoMass();
        }
        else
        {
          cerr << "Warning: cannot handle modification " << mod.getName() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
          continue;
        }
      }

      //cerr << "Setting fixed modification " << it->getModification() << " of amino acid '" << aa << "'; weight = " << aa_to_weight_[aa] << endl;

      const Residue* res = ResidueDB::getInstance()->getModifiedResidue(it->getModificationName());
      name_to_residue_[aa] = res;
      residue_to_name_[res] = aa;
    }

    const StringList mod_names(ListUtils::create<String>("a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z"));
    vector<String>::const_iterator actual_mod_name = mod_names.begin();
    set<ModificationDefinition> var_mods = mod_set.getVariableModifications();
    for (set<ModificationDefinition>::const_iterator it = var_mods.begin(); it != var_mods.end(); ++it)
    {
      ResidueModification mod = it->getModification();
      char aa = (*actual_mod_name)[0];
      char origin_aa = ' ';
      ++actual_mod_name;

      if (mod.getOrigin() == 'X')
      {
        cerr << "CompNovoIdentificationBase: Warning: cannot handle modification " << mod.getName() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        origin_aa = mod.getOrigin();
      }

      if (mod.getMonoMass() != 0)
      {
        aa_to_weight_[aa] = mod.getMonoMass();
      }
      else
      {
        if (mod.getDiffMonoMass() != 0)
        {
          aa_to_weight_[aa] = aa_to_weight_[origin_aa] + mod.getDiffMonoMass();
        }
        else
        {
          cerr << "CompNovoIdentificationBase: Warning: cannot handle modification " << mod.getName() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
          continue;
        }
      }

      //cerr << "Mapping variable modification " << it->getModification() << " to letter '" << aa << "' (@" << origin_aa << "); weight = " << aa_to_weight_[aa] << endl;
      const Residue* res = ResidueDB::getInstance()->getModifiedResidue(it->getModificationName());
      name_to_residue_[aa] = res;
      residue_to_name_[res] = aa;
    }

    /*
    cerr << "Following masses are used for identification: " << endl;

    for (Map<char, double>::const_iterator it = aa_to_weight_.begin(); it != aa_to_weight_.end(); ++it)
    {
        cerr << it->first << " " << precisionWrapper(it->second) << endl;
    }*/

    Param decomp_param(mass_decomp_algorithm_.getParameters());
    decomp_param.setValue("tolerance", fragment_mass_tolerance_);
    decomp_param.setValue("fixed_modifications", param_.getValue("fixed_modifications"));
    decomp_param.setValue("variable_modifications", param_.getValue("variable_modifications"));
    mass_decomp_algorithm_.setParameters(decomp_param);

    min_aa_weight_ = numeric_limits<double>::max();
    for (Map<char, double>::const_iterator it = aa_to_weight_.begin(); it != aa_to_weight_.end(); ++it)
    {
      if (min_aa_weight_ > it->second)
      {
        min_aa_weight_ = it->second;
      }
    }
  }

  void CompNovoIdentificationBase::getETDSpectrum_(PeakSpectrum & spec, const String & sequence, Size /* charge */, double prefix, double suffix)
  {
    if (isotope_distributions_.empty())
    {
      initIsotopeDistributions_();
    }

    Peak1D p;
    p.setIntensity(1.0f);

    double c_pos(17.0 + prefix);     // TODO high mass accuracy!!
    double z_pos(3.0 + suffix);
    //double b_pos(0.0 + prefix);
    //double y_pos(18.0 + suffix);
    // sometimes also b and y ions are in this spectrum

    #ifdef ETD_SPECTRUM_DEBUG
    cerr << "ETDSpectrum for " << sequence << " " << prefix << " " << suffix << endl;
    #endif

    for (Size i = 0; i != sequence.size() - 1; ++i)
    {
      char aa(sequence[i]);
      char aa_cterm(sequence[i + 1]);
      #ifdef ETD_SPECTRUM_DEBUG
      cerr << aa << " " << aa_cterm << endl;
      #endif

      c_pos += aa_to_weight_[aa];
      //b_pos += aa_to_weight_[aa];

      char aa2(sequence[sequence.size() - i - 1]);
      z_pos += aa_to_weight_[aa2];
      //y_pos += aa_to_weight_[aa2];

      #ifdef ETD_SPECTRUM_DEBUG
      cerr << b_pos << " " << c_pos << " " << y_pos << " " << z_pos << endl;
      #endif

      if (aa_cterm != 'P')
      {
        // c-ions
        if (c_pos + 1 >= min_mz_ && c_pos + 1 <= max_mz_)
        {
          //p.setIntensity(0.3);
          //p.setPosition(c_pos);
          //spec.push_back(p);
          for (Size j = 0; j != max_isotope_; ++j)
          {
            p.setIntensity(isotope_distributions_[(int)c_pos][j]);
            p.setPosition(c_pos + 1 + j);
            spec.push_back(p);
          }
        }
      }

      if (aa2 != 'P')
      {
        // z-ions
        if (z_pos >= min_mz_ && z_pos <= max_mz_)
        {
          p.setIntensity(0.3f);
          p.setPosition(z_pos);
          spec.push_back(p);

          for (Size j = 0; j != max_isotope_; ++j)
          {
            p.setIntensity(isotope_distributions_[(int)z_pos][j]);
            p.setPosition(z_pos + 1 + j);
            spec.push_back(p);
          }
        }
      }
    }

    spec.sortByPosition();

    #ifdef ETD_SPECTRUM_DEBUG
    for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
    {
      cerr << it->getPosition()[0] << " " << it->getIntensity() << endl;
    }
    #endif

    return;
  }

}
