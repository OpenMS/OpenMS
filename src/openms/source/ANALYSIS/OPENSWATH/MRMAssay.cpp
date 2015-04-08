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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>

namespace OpenMS
{
  MRMAssay::MRMAssay()
  {
  }

  MRMAssay::~MRMAssay()
  {
  }

  bool MRMAssay::isUIS_(const double fragment_ion, std::vector<double> ions, double mz_threshold)
  {
    size_t number_uis = 0;

    for (std::vector<double>::iterator i_it = ions.begin(); i_it != ions.end(); ++i_it)
    {
      if (*i_it - mz_threshold <= fragment_ion && *i_it + mz_threshold >= fragment_ion)
      {
        number_uis++;
      }
    }

    if (number_uis <= 1)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void MRMAssay::isSiteUIS_(const AASequence sequence, ReactionMonitoringTransition& tr)
  {
    std::vector<size_t> modified_residues;

    if (sequence.isModified())
    {
      for (size_t i = 0; i < sequence.size(); ++i)
      {
        if (sequence.isModified(i))
        {
          modified_residues.push_back(i + 1);
        }
      }
    }

    CVTermList annotation = tr.getProduct().getInterpretationList()[0];

    IntList modified_sites;
    StringList modified_sites_classes;

    for (std::vector<size_t>::iterator mod_it = modified_residues.begin(); mod_it != modified_residues.end(); ++mod_it)
    {
      // x, y or z ion
      if ((Size)annotation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt() == sequence.size() + 1 - *mod_it && (annotation.hasCVTerm("MS:1001228") || annotation.hasCVTerm("MS:1001220") || annotation.hasCVTerm("MS:1001230")))
      {
        modified_sites.push_back(*mod_it);
        modified_sites_classes.push_back("revdiag");
      }
      // neighboring (pos +1) x, y or z ion
      else if ((Size)annotation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt() == sequence.size() - *mod_it && (annotation.hasCVTerm("MS:1001228") || annotation.hasCVTerm("MS:1001220") || annotation.hasCVTerm("MS:1001230")))
      {
        modified_sites.push_back(*mod_it);
        modified_sites_classes.push_back("revnext");
      }
      // a, b or c ion
      else if ((Size)annotation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt() == *mod_it && (annotation.hasCVTerm("MS:1001229") || annotation.hasCVTerm("MS:1001224") || annotation.hasCVTerm("MS:1001231")))
      {
        modified_sites.push_back(*mod_it);
        modified_sites_classes.push_back("fwddiag");
      }
      // neighbouring (pos -1) ) a, b or c ion
      else if ((Size)annotation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt() == *mod_it - 1 && (annotation.hasCVTerm("MS:1001229") || annotation.hasCVTerm("MS:1001224") || annotation.hasCVTerm("MS:1001231")))
      {
        modified_sites.push_back(*mod_it);
        modified_sites_classes.push_back("fwdnext");
      }
    }

    tr.setMetaValue("site_identifying_transition", ListUtils::concatenate(modified_sites, ","));
    tr.setMetaValue("site_identifying_class", ListUtils::concatenate(modified_sites_classes, ","));
  }

  int MRMAssay::getSwath_(const std::vector<std::pair<double, double> > swathes, const double precursor_mz)
  {
    int swath = -1;

    // we go through all swaths in ascending order and if the transitions falls in overlap, only the upper swath will be used and checked.
    for (std::vector<std::pair<double, double> >::const_iterator it = swathes.begin(); it != swathes.end(); ++it)
    {
      if (precursor_mz >= it->first && precursor_mz <= it->second)
      {
        swath = it - swathes.begin();
      }
    }

    if (swath != -1)
    {
      return swath;
    }
    else
    {
      return -1;
    }
  }

  bool MRMAssay::isInSwath_(const std::vector<std::pair<double, double> > swathes, const double precursor_mz, const double product_mz)
  {
    int swath_idx = getSwath_(swathes, precursor_mz);

    if (swath_idx == -1) { return true; } // remove all transitions that are not in swath range
    else
    {
      std::pair<double, double> swath = swathes[getSwath_(swathes, precursor_mz)];

      if (product_mz >= swath.first && product_mz <= swath.second)
      {
        return true;
      }
      else { return false; }
    }
  }

  void MRMAssay::addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                                  int location, ResidueModification& rmod, const String& name)
  {
    TargetedExperiment::Peptide::Modification mod;
    String unimod_str = rmod.getUniModAccession();
    mod.location = location;
    mod.mono_mass_delta = rmod.getDiffMonoMass();
    mod.avg_mass_delta = rmod.getDiffAverageMass();
    // CV term with the full unimod accession number and name
    CVTerm unimod_name;
    unimod_name.setCVIdentifierRef("UNIMOD");
    unimod_name.setAccession(unimod_str.toUpper());
    unimod_name.setName(name);
    mod.addCVTerm(unimod_name);
    mods.push_back(mod);
  }

  void MRMAssay::reannotateTransitions(OpenMS::TargetedExperiment& exp, double mz_threshold, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_alternative_localizations, bool enable_reannotation, bool enable_losses, int round_decPow)
  {
    PeptideVectorType peptides;
    ProteinVectorType proteins;
    TransitionVectorType transitions;

    ModificationsDB* mod_db = ModificationsDB::getInstance();
    OpenMS::MRMIonSeries mrmis;

    size_t j = 0;
    for (size_t i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(tr.getPeptideRef());
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);

      std::vector<OpenMS::AASequence> alternative_peptide_sequences;

      if (enable_alternative_localizations)
      {
        alternative_peptide_sequences = combineModifications_(target_peptide_sequence);
      }
      else
      {
        alternative_peptide_sequences.push_back(target_peptide_sequence);
      }

      for (std::vector<OpenMS::AASequence>::iterator aa_it = alternative_peptide_sequences.begin(); aa_it != alternative_peptide_sequences.end(); ++aa_it)
      {
        std::vector<TargetedExperiment::Peptide::Modification> mods;

        if (!aa_it->getNTerminalModification().empty())
        {
          ResidueModification rmod = mod_db->getTerminalModification(aa_it->getNTerminalModification(), ResidueModification::N_TERM);
          addModification_(mods, -1, rmod, aa_it->getNTerminalModification());
        }
        if (!aa_it->getCTerminalModification().empty())
        {
          ResidueModification rmod = mod_db->getTerminalModification(aa_it->getCTerminalModification(), ResidueModification::C_TERM);
          addModification_(mods, aa_it->size(), rmod, aa_it->getCTerminalModification());
        }
        for (Size i = 0; i != aa_it->size(); ++i)
        {
          if (aa_it->isModified(i))
          {
            // search the residue in the modification database (if the sequence is valid, we should find it)
            TargetedExperiment::Peptide::Modification mod;
            ResidueModification rmod = mod_db->getModification(aa_it->getResidue(i).getOneLetterCode(),
                                                               aa_it->getResidue(i).getModification(), ResidueModification::ANYWHERE);
            addModification_(mods, i, rmod, aa_it->getResidue(i).getModification());
          }
        }

        target_peptide.mods = mods;
        target_peptide.setMetaValue("full_peptide_name", aa_it->toString());
        target_peptide.id = String(target_peptide.protein_refs[0]) + String("_") + TargetedExperimentHelper::getAASequence(target_peptide).toString() + String("_") + String(target_peptide.getChargeState()) + "_" + target_peptide.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString();

        mrmis.annotateTransition(tr, target_peptide, mz_threshold, enable_reannotation, fragment_types, fragment_charges, enable_losses);

        if (tr.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001240"))
        {
          LOG_DEBUG << "[unannotated] Skipping " << aa_it->toString() << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << " " << tr.getMetaValue("annotation") << std::endl;
          continue;
        }
        else
        {
          LOG_DEBUG << "[selected] " << aa_it->toString() << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << " " << tr.getMetaValue("annotation") << std::endl;
        }

        if (std::find(peptides.begin(), peptides.end(), target_peptide) == peptides.end())
        {
          peptides.push_back(target_peptide);
          LOG_DEBUG << "[selected] " <<  aa_it->toString() << std::endl;
        }

        tr.setPeptideRef(target_peptide.id);
        tr.setPrecursorMZ(Math::roundDecimal(target_peptide_sequence.getMonoWeight(Residue::Full, target_peptide.getChargeState()) / target_peptide.getChargeState(), round_decPow));
        tr.setNativeID(String(j) + String("_") +  String(target_peptide.protein_refs[0]) + String("_") + target_peptide.sequence + String("_") + String(tr.getPrecursorMZ()) + "_" + String(tr.getProductMZ()));

        j += 1;
        transitions.push_back(tr);
      }
    }

    exp.setTransitions(transitions);
    exp.setPeptides(peptides);
  }

  void MRMAssay::restrictTransitions(OpenMS::TargetedExperiment& exp, double lower_mz_limit, double upper_mz_limit, std::vector<std::pair<double, double> > swathes)
  {
    OpenMS::MRMIonSeries mrmis;
    MRMDecoy::PeptideVectorType peptides;
    MRMDecoy::ProteinVectorType proteins;
    MRMDecoy::TransitionVectorType transitions;

    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(tr.getPeptideRef());
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);

      if (tr.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001240"))
      {
        LOG_DEBUG << "[unannotated] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << " " << tr.getMetaValue("annotation") << std::endl;
        continue;
      }

      if (swathes.size() > 0)
      {
        if (MRMAssay::isInSwath_(swathes, tr.getPrecursorMZ(), tr.getProductMZ()))
        {
          LOG_DEBUG << "[swath] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << std::endl;
          continue;
        }
      }

      if (tr.getProductMZ() < lower_mz_limit || tr.getProductMZ() > upper_mz_limit)
      {
        LOG_DEBUG << "[mz_limit] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << std::endl;
        continue;
      }

      transitions.push_back(tr);
    }

    exp.setTransitions(transitions);
  }

  void MRMAssay::detectingTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions)
  {
    MRMDecoy::PeptideVectorType peptides;
    std::vector<String> peptide_ids;
    MRMDecoy::ProteinVectorType proteins;
    MRMDecoy::TransitionVectorType transitions;

    Map<String, MRMDecoy::TransitionVectorType> TransitionsMap;

    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      if (TransitionsMap.find(tr.getPeptideRef()) == TransitionsMap.end())
      {
        TransitionsMap[tr.getPeptideRef()];
      }

      TransitionsMap[tr.getPeptideRef()].push_back(tr);
    }

    for (Map<String, MRMDecoy::TransitionVectorType>::iterator m = TransitionsMap.begin();
         m != TransitionsMap.end(); ++m)
    {
      if (m->second.size() >= (Size)min_transitions)
      {
        std::vector<double> LibraryIntensity;
        for (MRMDecoy::TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          LibraryIntensity.push_back(boost::lexical_cast<double>(tr_it->getLibraryIntensity()));
        }

        // sort by intensity, reverse and delete all elements after max_transitions
        std::sort(LibraryIntensity.begin(), LibraryIntensity.end());
        std::reverse(LibraryIntensity.begin(), LibraryIntensity.end());
        if ((Size)max_transitions < LibraryIntensity.size())
        {
          std::vector<double>::iterator start_delete = LibraryIntensity.begin();
          std::advance(start_delete, max_transitions);
          LibraryIntensity.erase(start_delete, LibraryIntensity.end());
        }

        for (MRMDecoy::TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          ReactionMonitoringTransition tr = *tr_it;

          if ((std::find(LibraryIntensity.begin(), LibraryIntensity.end(), boost::lexical_cast<double>(tr.getLibraryIntensity())) != LibraryIntensity.end()) && tr.getDecoyTransitionType() != ReactionMonitoringTransition::DECOY)
          {
            tr.setMetaValue("detecting_transition", "true");
          }
          else
          {
            continue;
          }

          transitions.push_back(tr);
          if (std::find(peptide_ids.begin(), peptide_ids.end(), tr.getPeptideRef()) == peptide_ids.end())
          {
            peptide_ids.push_back(tr.getPeptideRef());
          }
        }
      }
    }

    std::vector<String> ProteinList;
    for (Size i = 0; i < exp.getPeptides().size(); ++i)
    {
      TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      if (std::find(peptide_ids.begin(), peptide_ids.end(), peptide.id) != peptide_ids.end())
      {
        peptides.push_back(peptide);
        for (Size j = 0; j < peptide.protein_refs.size(); ++j)
        {
          ProteinList.push_back(peptide.protein_refs[j]);
        }
      }
      else
      {
        LOG_DEBUG << "[peptide] Skipping " << peptide.id << std::endl;
      }
    }

    for (Size i = 0; i < exp.getProteins().size(); ++i)
    {
      OpenMS::TargetedExperiment::Protein protein = exp.getProteins()[i];
      if (find(ProteinList.begin(), ProteinList.end(), protein.id) != ProteinList.end())
      {
        proteins.push_back(protein);
      }
      else
      {
        LOG_DEBUG << "[protein] Skipping " << protein.id << std::endl;
      }
    }

    exp.setTransitions(transitions);
    exp.setPeptides(peptides);
    exp.setProteins(proteins);
  }

  void MRMAssay::uisTransitions(OpenMS::TargetedExperiment& exp, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses, bool enable_uis_scoring, bool enable_site_scoring, double mz_threshold, std::vector<std::pair<double, double> > swathes)
  {
    OpenMS::MRMDecoy mrmdg;
    OpenMS::MRMIonSeries mrmis;

    MRMDecoy::TransitionVectorType transitions;
    Map<size_t, std::vector<double> > IonMap;
    Map<String, MRMIonSeries::IonSeries> DecoyIonSeriesMap;

    for (size_t i = 0; i < exp.getPeptides().size(); ++i)
    {
      TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      OpenMS::AASequence peptide_sequence = TargetedExperimentHelper::getAASequence(peptide);
      int precursor_swath = getSwath_(swathes, peptide_sequence.getMonoWeight(Residue::Full, peptide.getChargeState()) / peptide.getChargeState());

      if (IonMap.find(precursor_swath) == IonMap.end())
      {
        IonMap[precursor_swath];
      }

      MRMIonSeries::IonSeries ionseries = mrmis.getIonSeries(peptide_sequence, peptide.getChargeState(), fragment_types, fragment_charges, enable_losses);
      std::vector<double> ions;
      std::transform(ionseries.begin(), ionseries.end(), std::back_inserter(ions), boost::bind(&MRMIonSeries::IonSeries::value_type::second, _1));

      IonMap[precursor_swath].insert(IonMap[precursor_swath].end(), ions.begin(), ions.end());

      DecoyIonSeriesMap[peptide_sequence.toString() + peptide.getChargeState()] = mrmis.getIonSeries(TargetedExperimentHelper::getAASequence(mrmdg.reversePeptide(peptide)), peptide.getChargeState(), fragment_types, fragment_charges, enable_losses);
    }

    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(tr.getPeptideRef());
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
      int target_precursor_swath = getSwath_(swathes, target_peptide_sequence.getMonoWeight(Residue::Full, target_peptide.getChargeState()) / target_peptide.getChargeState());

      if (tr.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001240"))
      {
        LOG_DEBUG << "[unannotated] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << " " << tr.getMetaValue("annotation") << std::endl;
        continue;
      }
      OpenMS::AASequence target_peptide_sequence_unmod = AASequence::fromString(TargetedExperimentHelper::getAASequence(target_peptide).toUnmodifiedString());

      if (!enable_uis_scoring)
      {
        tr.setMetaValue("identifying_transition", "false");
      }
      else if (isUIS_(tr.getProductMZ(), IonMap[target_precursor_swath], mz_threshold))
      {
        tr.setMetaValue("identifying_transition", "true");
      }
      else if (!isUIS_(tr.getProductMZ(), IonMap[target_precursor_swath], mz_threshold))
      {
        tr.setMetaValue("identifying_transition", "false");
      }

      if (enable_site_scoring)
      {
        isSiteUIS_(target_peptide_sequence, tr);
      }
      else
      {
        tr.setMetaValue("site_identifying_transition", "");
        tr.setMetaValue("site_identifying_class", "");
      }

      if (tr.getMetaValue("detecting_transition").toBool() || tr.getMetaValue("identifying_transition").toBool() || tr.getMetaValue("site_identifying_transition").toString() != "")
      {
        transitions.push_back(tr);
      }

      if (tr.getMetaValue("identifying_transition").toBool() || tr.getMetaValue("site_identifying_transition").toString() != "")
      {
        ReactionMonitoringTransition uisdecoy_tr = tr;
        uisdecoy_tr.setMetaValue("detecting_transition", "false");
        uisdecoy_tr.setName("UISDECOY_" + uisdecoy_tr.getName());

        std::pair<String, double> uisdecoy_ion = mrmis.getIon(DecoyIonSeriesMap[target_peptide_sequence.toString() + target_peptide.getChargeState()], tr.getMetaValue("annotation"));

        uisdecoy_tr.setProductMZ(uisdecoy_ion.second);
        uisdecoy_tr.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);

        transitions.push_back(uisdecoy_tr);
      }
    }

    exp.setTransitions(transitions);
  }

  void MRMAssay::insilicoTransitions(OpenMS::TargetedExperiment& exp, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses, int round_decPow)
  {
    OpenMS::MRMIonSeries mrmis;
    std::map<std::string, std::map<size_t, std::vector<String> > > TransitionsMzMap;
    std::map<std::string, std::map<size_t, TargetedExperimentHelper::RetentionTime> > TransitionsRtMap;
    std::map<std::string, std::map<size_t, String> > TransitionsPepMap;
    MRMAssay::TransitionVectorType insilico_transitions;

    for (size_t i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(tr.getPeptideRef());
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);

      if (TransitionsMzMap.find(target_peptide_sequence.toString()) == TransitionsMzMap.end())
      {
        TransitionsMzMap[target_peptide_sequence.toString()];
        TransitionsRtMap[target_peptide_sequence.toString()];
        TransitionsPepMap[target_peptide_sequence.toString()];
      }

      if (TransitionsMzMap[target_peptide_sequence.toString()].find(target_peptide.getChargeState()) == TransitionsMzMap[target_peptide_sequence.toString()].end())
      {
        TransitionsMzMap[target_peptide_sequence.toString()][target_peptide.getChargeState()];
      }

      TransitionsMzMap[target_peptide_sequence.toString()][target_peptide.getChargeState()].push_back(tr.getMetaValue("annotation"));
      TransitionsRtMap[target_peptide_sequence.toString()][target_peptide.getChargeState()] = tr.getRetentionTime();
      TransitionsPepMap[target_peptide_sequence.toString()][target_peptide.getChargeState()] = tr.getPeptideRef();
      insilico_transitions.push_back(tr);
    }

    for (std::map<std::string, std::map<size_t, std::vector<String> > >::iterator mz_it = TransitionsMzMap.begin(); mz_it != TransitionsMzMap.end(); ++mz_it)
    {
      for (std::map<size_t, std::vector<String> >::iterator pr_it = mz_it->second.begin(); pr_it != mz_it->second.end(); ++pr_it)
      {
        MRMIonSeries::IonSeries ionseries = mrmis.getIonSeries(AASequence::fromString(mz_it->first), pr_it->first, fragment_types, fragment_charges, enable_losses);

        for (MRMIonSeries::IonSeries::iterator ion_it = ionseries.begin(); ion_it != ionseries.end(); ++ion_it)
        {
          if (std::find(pr_it->second.begin(), pr_it->second.end(), ion_it->first) == pr_it->second.end())
          {
            ReactionMonitoringTransition trn;
            mrmis.annotateTransitionCV(trn, ion_it->first);
            trn.setName(String("insilico") + "_" + mz_it->first + "_" + String(pr_it->first) + "_" + ion_it->first + "_" + String(ion_it->second));
            trn.setNativeID(String("insilico") + "_" + mz_it->first + "_" + String(pr_it->first) + "_" + ion_it->first + "_" + String(ion_it->second));
            trn.setMetaValue("detecting_transition", "false");
            trn.setMetaValue("insilico_transition", "true");
            trn.setPrecursorMZ(Math::roundDecimal(AASequence::fromString(mz_it->first).getMonoWeight(Residue::Full, pr_it->first) / pr_it->first, round_decPow));
            trn.setProductMZ(ion_it->second);
            trn.setRetentionTime(TransitionsRtMap[mz_it->first][pr_it->first]);
            trn.setPeptideRef(TransitionsPepMap[mz_it->first][pr_it->first]);

            insilico_transitions.push_back(trn);
          }
        }
      }
    }
    exp.setTransitions(insilico_transitions);
  }

  std::vector<std::vector<size_t> > MRMAssay::nchoosekcombinations_(std::vector<size_t> n, size_t k)
  {
    std::vector<std::vector<size_t> > combinations;

    std::string bitmask(k, 1);
    bitmask.resize(n.size(), 0);

    do
    {
      std::vector<size_t> combination;
      for (size_t i = 0; i < n.size(); ++i)
      {
        if (bitmask[i])
        {
          combination.push_back(n[i]);
        }
      }
      combinations.push_back(combination);
    }
    while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
  }

  std::vector<OpenMS::AASequence> MRMAssay::combineModifications_(OpenMS::AASequence sequence)
  {
    std::map<OpenMS::String, size_t> mods;
    std::vector<OpenMS::AASequence> sequences;
    sequences.push_back(AASequence::fromString(sequence.toUnmodifiedString()));

    OpenMS::ModificationsDB* ptr = ModificationsDB::getInstance();

    for (size_t i = 0; i < sequence.size(); ++i)
    {
      if (sequence.isModified(i))
      {
        mods[sequence.getResidue(i).getModification()] += 1;
      }
    }

    for (std::map<OpenMS::String, size_t>::iterator mod_it = mods.begin(); mod_it != mods.end(); ++mod_it)
    {
      std::vector<size_t> mods_res;
      for (size_t i = 0; i < sequence.size(); ++i)
      {
        std::set<const ResidueModification*> modifiable_residues;
        ptr->searchModifications(modifiable_residues, sequence.getResidue(i).getOneLetterCode(), mod_it->first, ResidueModification::ANYWHERE);
        if (!modifiable_residues.empty())
        {
          mods_res.push_back(i);
        }
      }

      std::vector<std::vector<size_t> > mods_combs = nchoosekcombinations_(mods_res, mod_it->second);
      sequences = addModificationsSequences_(sequences, mods_combs, mod_it->first);
    }

    return sequences;
  }

  std::vector<OpenMS::AASequence> MRMAssay::addModificationsSequences_(std::vector<OpenMS::AASequence> sequences, std::vector<std::vector<size_t> > mods_combs, OpenMS::String modification)
  {
    std::vector<OpenMS::AASequence> modified_sequences;
    bool multi_mod_switch = false;

    for (std::vector<OpenMS::AASequence>::iterator sq_it = sequences.begin(); sq_it != sequences.end(); ++sq_it)
    {
      for (std::vector<std::vector<size_t> >::iterator mc_it = mods_combs.begin(); mc_it != mods_combs.end(); ++mc_it)
      {
        multi_mod_switch = false;
        OpenMS::AASequence temp_sequence = *sq_it;
        for (std::vector<size_t>::iterator pos_it = mc_it->begin(); pos_it != mc_it->end(); ++pos_it)
        {
          if (!temp_sequence.isModified(*pos_it))
          {
            temp_sequence.setModification(*pos_it, modification);
          }
          else
          {
            multi_mod_switch = true;
          }
        }
        if (!multi_mod_switch) { modified_sequences.push_back(temp_sequence); }
      }
    }

    return modified_sequences;
  }

}
