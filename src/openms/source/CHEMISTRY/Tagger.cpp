// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Timo Sachsenberg, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  char Tagger::getAAByMass_(double m) const
  {
    // fast check for border cases
    if (m < min_gap_ || m > max_gap_) return ' ';

    const double delta = Math::ppmToMass(ppm_, m);
    auto left = mass2aa_.lower_bound(m - delta);
    //if (left == mass2aa_.end()) return ' '; // cannot happen, since we checked boundaries above

    if (fabs(left->first - m) < delta) return left->second;
    return ' ';
  }

  void Tagger::getTag_(std::string & tag, const std::vector<double>& mzs, const size_t i, std::vector<std::string>& tags, const size_t charge) const
  {
    const size_t N = mzs.size();
    size_t j = i + 1;
    while (j < N)
    {
      if (tag.size() == max_tag_length_) 
      { 
        return; // maximum tag size reached? - continue with next parent
      }
      const double gap = mzs[j] - mzs[i];
      if ((gap * charge) > max_gap_) 
      { 
        return; // already too far away - continue with next parent 
      } 

      const char aa = getAAByMass_(gap * charge);
      if (aa == ' ') { ++j; continue; } // can't extend tag

      tag += aa;

      if (tag.size() >= min_tag_length_)
      {
          tags.push_back(tag);
      }

      getTag_(tag, mzs, j, tags, charge);

      // if aa is "L", then also add "I" as an alternative residue and extend the tag again
      // this will add redundancy, (and redundant runtime) but we avoid dealing with J and ambiguous matching to I and L later on
      if (aa == 'L')
      {
        tag.pop_back();
        tag.push_back('I');

        if (tag.size() >= min_tag_length_)
        {
          tags.push_back(tag);
        }
        getTag_(tag, mzs, j, tags, charge);
      }
      tag.pop_back();  // remove last string
      ++j;
    }
  }

  Tagger::Tagger(size_t min_tag_length, double ppm, size_t max_tag_length, size_t min_charge, size_t max_charge, const StringList& fixed_mods, const StringList& var_mods)
    : ppm_{fabs(ppm)}, min_tag_length_{min_tag_length}, max_tag_length_{max_tag_length}, min_charge_{min_charge}, max_charge_{max_charge}
  {
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

    for (const auto& r : aas)
    {
      const char letter = r->getOneLetterCode()[0];
      const double mass = r->getMonoWeight(Residue::Internal);
      mass2aa_[mass] = letter;
    }

    // for fixed modifications, replace the unmodified residue with the modified one
    for (const auto& mod : fixed_mods)
    {
      const ResidueModification* rm = ModificationsDB::getInstance()->getModification(mod);
      Residue r = *(ResidueDB::getInstance()->getResidue(rm->getOrigin()));
      r.setModification(rm->getId());

      // remove the unmodified residue
      // this requires searching the map by value, but this is only done once when the Tagger is initialized
      for (std::map<double, char>::iterator it = mass2aa_.begin(); it != mass2aa_.end(); ++it)
      {
        if (it->second == rm->getOrigin())
        {
          mass2aa_.erase(it);
          break;
        }
      }
      const char name = rm->getOrigin();
      const double mass = r.getMonoWeight(Residue::Internal);
      mass2aa_[mass] = name;
    }

    // for variable modifications, add an additional instance of the residue with the modified mass to the list
    for (const auto& mod : var_mods)
    {
      const ResidueModification* rm = ModificationsDB::getInstance()->getModification(mod);
      Residue r = *(ResidueDB::getInstance()->getResidue(rm->getOrigin()));
      r.setModification(rm->getId());
      const char name = rm->getOrigin();
      const double mass = r.getMonoWeight(Residue::Internal);
      mass2aa_[mass] = name;
    }

    min_gap_ = mass2aa_.begin()->first - Math::ppmToMass(ppm, mass2aa_.begin()->first);
    max_gap_ = mass2aa_.rbegin()->first + Math::ppmToMass(ppm, mass2aa_.rbegin()->first);
  }

  void Tagger::getTag(const std::vector<double>& mzs, std::vector<std::string>& tags) const
  {
    // start peak
    if (min_tag_length_ > mzs.size()) return; // avoid segfault

#pragma omp parallel
    {
      std::vector<std::string> tags_local;
#pragma omp for schedule(guided)
      for (int i = 0; i < static_cast<int>(mzs.size() - min_tag_length_); ++i)
      {
        for (size_t charge = min_charge_; charge <= max_charge_; ++charge)
        {
          std::string tag;
          getTag_(tag, mzs, i, tags_local, charge);
        }
      } // end of loop over starting peaks
#pragma omp critical (join_tags)
      tags.insert(tags.end(), tags_local.begin(), tags_local.end());
    } // end of parallel section

    // make tags unique
    sort(tags.begin(), tags.end());
    auto last_unique_tag = unique(tags.begin(), tags.end());
    if (last_unique_tag != tags.end())
    {
      tags.erase(last_unique_tag, tags.end());
    }
  }

  void Tagger::getTag(const MSSpectrum& spec, std::vector<std::string>& tags) const
  {
    const size_t N = spec.size();
    if (N < min_tag_length_) { return; }
    // copy to double vector (speed)
    std::vector<double> mzs;
    mzs.reserve(N);
    for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
    getTag(mzs, tags);
  }
  void Tagger::setMaxCharge(size_t max_charge)
  {
    max_charge_ = max_charge;
  }
}
