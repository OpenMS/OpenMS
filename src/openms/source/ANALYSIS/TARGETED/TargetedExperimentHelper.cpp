// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS::TargetedExperimentHelper
  {

    void setModification(int location, int max_size, const String& modification, OpenMS::AASequence& aas)
    {
      OPENMS_PRECONDITION(location >= -1 && location <= max_size, 
          (String("Location has invalid value") + (String)location).c_str() )

      if (location == -1)
      {
        aas.setNTerminalModification(modification);
      }
      else if (location == max_size)
      {
        aas.setCTerminalModification(modification);
      }
      else
      {
        aas.setModification(location, modification);
      }
    }

    OpenMS::AASequence getAASequence(const OpenMS::TargetedExperiment::Peptide& peptide)
    {

      // Note that the peptide.sequence is the "naked sequence" without any
      // modifications on it, therefore we have to populate the AASequence with
      // the correct modifications afterwards.
      OpenMS::ModificationsDB* mod_db = OpenMS::ModificationsDB::getInstance();
      OpenMS::AASequence aas = AASequence::fromString(peptide.sequence);

      // Populate the AASequence with the correct modifications derived from
      // the Peptide::Modification objects.
      for (std::vector<Peptide::Modification>::const_iterator it = peptide.mods.begin(); 
          it != peptide.mods.end(); ++it)
      {
        // Step 1: First look whether the UniMod ID is set (we don't use a CVTerm any more but a member)
        if (it->unimod_id != -1)
        {
          setModification(it->location, boost::numeric_cast<int>(peptide.sequence.size()), 
              "UniMod:" + String(it->unimod_id), aas);
          continue;
        }

        OPENMS_LOG_WARN << "Warning: No UniMod id set for modification on peptide " << peptide.sequence << 
          ". Will try to infer modification id by mass next." << std::endl;

        // compare with code in source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.cpp

        // Step 2: If the above step fails, try to find the correct
        // modification by using the mass difference
        const ResidueModification* mod = mod_db->getBestModificationByDiffMonoMass(
          it->mono_mass_delta, 1.0, peptide.sequence[it->location]);
        if (mod != nullptr)
        {
          setModification(it->location, boost::numeric_cast<int>(peptide.sequence.size()), mod->getId(), aas);
        }
        else
        {
          // could not find any modification ...
          std::cerr << "Warning: Could not determine modification with delta mass " <<
            it->mono_mass_delta << " for peptide " << peptide.sequence <<
            " at position " << it->location << std::endl;
          std::cerr << "Skipping this modification" << std::endl;
        }
      }
      return aas;
    }

  }
