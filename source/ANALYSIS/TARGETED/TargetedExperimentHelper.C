// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
  namespace TargetedExperimentHelper
  {

    void setModification(int location, int max_size, String modification, OpenMS::AASequence& aas)
    {
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
      OpenMS::ModificationsDB* mod_db = OpenMS::ModificationsDB::getInstance();
      OpenMS::AASequence aas = AASequence(peptide.sequence);

      for (std::vector<OpenMS::TargetedExperiment::Peptide::Modification>::const_iterator it = peptide.mods.begin(); it != peptide.mods.end(); ++it)
      {
        // Step 1: First look for a cv term that says which unimod nr it is...
        // compare with code in source/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.C
        int nr_modifications_added = 0;
        Map<String, std::vector<CVTerm> > cv_terms = it->getCVTerms();
        for (Map<String, std::vector<CVTerm> >::iterator li = cv_terms.begin(); li != cv_terms.end(); ++li)
        {
          std::vector<CVTerm> mods = (*li).second;
          for (std::vector<CVTerm>::iterator mo = mods.begin(); mo != mods.end(); ++mo)
          {
            // if we find a CV term that starts with UniMod, chances are we can use the UniMod accession number
            if (mo->getAccession().size() > 7 && mo->getAccession().prefix(7).toLower() == String("unimod:"))
            {
              nr_modifications_added++;
              setModification(it->location, boost::numeric_cast<int>(peptide.sequence.size()), "UniMod:" + mo->getAccession().substr(7), aas);
            }
          }
        }

        // Step 2: If the above step fails, try to find the correct
        // modification some other way, i.e. by using the mass difference
        if (nr_modifications_added == 0)
        {
          std::vector<String> mods;
          mod_db->getModificationsByDiffMonoMass(mods, peptide.sequence[it->location], it->mono_mass_delta, 0.0);
          for (std::vector<String>::iterator mo = mods.begin(); mo != mods.end(); ++mo)
          {
            nr_modifications_added++;
            setModification(it->location, boost::numeric_cast<int>(peptide.sequence.size()), *mo, aas);
          }
        }

        // In theory, each modification in the TraML should map to one modification added to the AASequence
        if (nr_modifications_added > 1)
        {
          std::cout << "Warning: More than one modification was found for peptide " << peptide.sequence << " at position " << it->location << std::endl;
        }
        else if (nr_modifications_added == 0)
        {
          std::cout << "Warning: Could not determine modification with delta mass " <<  it->mono_mass_delta << " for peptide " << peptide.sequence << " at position " << it->location << std::endl;
        }
      }
      return aas;
    }

  }
}
