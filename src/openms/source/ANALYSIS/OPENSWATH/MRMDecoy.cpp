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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <map>
#include <utility> //for pair
#include <string>
#include <vector>
#include <algorithm> // for sort

namespace OpenMS
{
  std::pair<String, DoubleReal> MRMDecoy::getDecoyIon(String ionid, boost::unordered_map<String, boost::unordered_map<String, DoubleReal> >& decoy_ionseries)
  {
    using namespace boost::assign;
    // Select SpectraST Style
    std::vector<String> SpectraST_order;
    SpectraST_order += "b", "b_loss", "y", "y_loss", "a", "b_isotopes", "b_isotopes_loss", "y_isotopes", "y_isotopes_loss", "a_isotopes";

    // Iterate over ion type and then ordinal
    std::pair<String, DoubleReal> ion;
    String unannotated = "unannotated";
    ion = make_pair(unannotated, -1);
    for (std::vector<String>::iterator iontype = SpectraST_order.begin(); iontype != SpectraST_order.end(); ++iontype)
    {
      for (boost::unordered_map<String, DoubleReal>::iterator ordinal = decoy_ionseries[*iontype].begin(); ordinal != decoy_ionseries[*iontype].end(); ++ordinal)
      {
        if (ordinal->first == ionid)
        {
          ion = make_pair(ordinal->first, ordinal->second);
        }
      }
    }
    return ion;
  }

  std::pair<String, double> MRMDecoy::getTargetIon(double ProductMZ, double mz_threshold, boost::unordered_map<String, boost::unordered_map<String, double> > target_ionseries, bool enable_losses, bool enable_isotopes)
  {
    // make sure to only use annotated transitions and to use the theoretical MZ
    using namespace boost::assign;
    // Select SpectraST Style
    std::vector<String> SpectraST_order;
    SpectraST_order += "b", "y", "a";

    if (enable_losses) {
      SpectraST_order += "b_loss", "y_loss";
    }
    if (enable_isotopes) {
      SpectraST_order += "b_isotopes", "y_isotopes", "a_isotopes";
    }
    if (enable_isotopes && enable_losses) {
      SpectraST_order += "b_isotopes_loss", "y_isotopes_loss";
    }

    // Iterate over ion type and then ordinal
    std::pair<String, double> ion;
    String unannotated = "unannotated";
    ion = make_pair(unannotated, -1);
    double closest_delta = std::numeric_limits<double>::max();
    for (std::vector<String>::iterator iontype = SpectraST_order.begin(); iontype != SpectraST_order.end(); ++iontype)
    {
      for (boost::unordered_map<String, double>::iterator ordinal = target_ionseries[*iontype].begin(); ordinal != target_ionseries[*iontype].end(); ++ordinal)
      {
        if (std::fabs(ordinal->second - ProductMZ) <= mz_threshold && std::fabs(ordinal->second - ProductMZ) <= closest_delta)
        {
          closest_delta = std::fabs(ordinal->second - ProductMZ);
          ion = make_pair(ordinal->first, ordinal->second);
          break;
        }
      }
    }
    return ion;
  }

  boost::unordered_map<String, boost::unordered_map<String, double> > MRMDecoy::getIonSeries(AASequence sequence, int precursor_charge, int max_isotope)
  {
    boost::unordered_map<String, boost::unordered_map<String, double> > ionseries;
    boost::unordered_map<String, double> bionseries, bionseries_isotopes, bionseries_loss,
                             bionseries_isotopes_loss, yionseries, yionseries_isotopes,
                             yionseries_loss, yionseries_isotopes_loss, aionseries,
                             aionseries_isotopes;

    for (int charge = 1; charge <= precursor_charge; ++charge)
    {
      for (Size i = 1; i < sequence.size(); ++i)
      {
        AASequence ion = sequence.getPrefix(i);
        double pos = ion.getMonoWeight(Residue::BIon, charge) / (double) charge;

        bionseries["b" + String(i) + "/" + String(charge) + "+"] = pos;

        IsotopeDistribution dist = ion.getFormula(Residue::BIon, charge).getIsotopeDistribution(max_isotope);
        UInt j(0);
        for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
        {
          if (j > 0)
          {
            bionseries_isotopes["b" + String(i) + "/" + String(charge) + "+" + "/" + "iso" + String(j)] = ((double) (
                                                                                                             pos + (double) j * Constants::NEUTRON_MASS_U) / (double) charge);
          }
        }
      }
    }
    ionseries["b"] = bionseries;
    ionseries["b_isotopes"] = bionseries_isotopes;

    for (boost::unordered_map<String, double>::iterator bit = bionseries.begin(); bit != bionseries.end(); ++bit)
    {
      bionseries_loss[bit->first + "-17"] = bit->second - 17;
      bionseries_loss[bit->first + "-18"] = bit->second - 18;
      bionseries_loss[bit->first + "-34"] = bit->second - 34;
      bionseries_loss[bit->first + "-35"] = bit->second - 35;
      bionseries_loss[bit->first + "-36"] = bit->second - 36;
      bionseries_loss[bit->first + "-44"] = bit->second - 44;
      bionseries_loss[bit->first + "-45"] = bit->second - 45;
      bionseries_loss[bit->first + "-46"] = bit->second - 46;
      bionseries_loss[bit->first + "-64"] = bit->second - 64;
      bionseries_loss[bit->first + "-98"] = bit->second - 98;
    }
    ionseries["b_loss"] = bionseries_loss;
    for (boost::unordered_map<String, double>::iterator bit = bionseries_isotopes.begin(); bit != bionseries_isotopes.end(); ++bit)
    {
      bionseries_isotopes_loss[bit->first + "-17"] = bit->second - 17;
      bionseries_isotopes_loss[bit->first + "-18"] = bit->second - 18;
      bionseries_isotopes_loss[bit->first + "-34"] = bit->second - 34;
      bionseries_isotopes_loss[bit->first + "-35"] = bit->second - 35;
      bionseries_isotopes_loss[bit->first + "-36"] = bit->second - 36;
      bionseries_isotopes_loss[bit->first + "-44"] = bit->second - 44;
      bionseries_isotopes_loss[bit->first + "-45"] = bit->second - 45;
      bionseries_isotopes_loss[bit->first + "-46"] = bit->second - 46;
      bionseries_isotopes_loss[bit->first + "-64"] = bit->second - 64;
      bionseries_isotopes_loss[bit->first + "-98"] = bit->second - 98;
    }
    ionseries["b_isotopes_loss"] = bionseries_isotopes_loss;

    for (int charge = 1; charge <= precursor_charge; ++charge)
    {
      for (Size i = 1; i < sequence.size(); ++i)
      {
        AASequence ion = sequence.getSuffix(i);
        double pos = ion.getMonoWeight(Residue::YIon, charge) / (double) charge;
        yionseries["y" + String(i) + "/" + String(charge) + "+"] = pos;

        IsotopeDistribution dist = ion.getFormula(Residue::YIon, charge).getIsotopeDistribution(max_isotope);
        UInt j(0);
        for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
        {
          if (j > 0)
          {
            yionseries_isotopes["y" + String(i) + "/" + String(charge) + "+" + "/" + "iso" + String(j)] = ((double) (
                                                                                                             pos + (double) j * Constants::NEUTRON_MASS_U) / (double) charge);
          }
        }
      }
    }
    ionseries["y"] = yionseries;
    ionseries["y_isotopes"] = yionseries_isotopes;

    for (boost::unordered_map<String, double>::iterator yit = yionseries.begin(); yit != yionseries.end(); ++yit)
    {
      yionseries_loss[yit->first + "-17"] = yit->second - 17;
      yionseries_loss[yit->first + "-18"] = yit->second - 18;
      yionseries_loss[yit->first + "-34"] = yit->second - 34;
      yionseries_loss[yit->first + "-35"] = yit->second - 35;
      yionseries_loss[yit->first + "-36"] = yit->second - 36;
      yionseries_loss[yit->first + "-44"] = yit->second - 44;
      yionseries_loss[yit->first + "-45"] = yit->second - 45;
      yionseries_loss[yit->first + "-46"] = yit->second - 46;
      yionseries_loss[yit->first + "-64"] = yit->second - 64;
      yionseries_loss[yit->first + "-98"] = yit->second - 98;
    }
    ionseries["y_loss"] = yionseries_loss;
    for (boost::unordered_map<String, double>::iterator yit = yionseries_isotopes.begin(); yit != yionseries_isotopes.end(); ++yit)
    {
      yionseries_isotopes_loss[yit->first + "-17"] = yit->second - 17;
      yionseries_isotopes_loss[yit->first + "-18"] = yit->second - 18;
      yionseries_isotopes_loss[yit->first + "-34"] = yit->second - 34;
      yionseries_isotopes_loss[yit->first + "-35"] = yit->second - 35;
      yionseries_isotopes_loss[yit->first + "-36"] = yit->second - 36;
      yionseries_isotopes_loss[yit->first + "-44"] = yit->second - 44;
      yionseries_isotopes_loss[yit->first + "-45"] = yit->second - 45;
      yionseries_isotopes_loss[yit->first + "-46"] = yit->second - 46;
      yionseries_isotopes_loss[yit->first + "-64"] = yit->second - 64;
      yionseries_isotopes_loss[yit->first + "-98"] = yit->second - 98;
    }
    ionseries["y_isotopes_loss"] = yionseries_isotopes_loss;

    for (int charge = 1; charge <= precursor_charge; ++charge)
    {
      for (Size i = 1; i < sequence.size(); ++i)
      {
        AASequence ion = sequence.getPrefix(i);
        double pos = ion.getMonoWeight(Residue::AIon, charge) / (double) charge;

        aionseries["a" + String(i) + "/" + String(charge) + "+"] = pos;

        IsotopeDistribution dist = ion.getFormula(Residue::AIon, charge).getIsotopeDistribution(max_isotope);
        UInt j(0);
        for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
        {
          if (j > 0)
          {
            aionseries_isotopes["a" + String(i) + "/" + String(charge) + "+" + "/" + "iso" + String(j)] = ((double) (
                                                                                                             pos + (double) j * Constants::NEUTRON_MASS_U) / (double) charge);
          }
        }
      }
    }
    ionseries["a"] = aionseries;
    ionseries["a_isotopes"] = aionseries_isotopes;

    return ionseries;
  }

  std::vector<std::pair<std::string::size_type, std::string> > MRMDecoy::find_all_tryptic(std::string sequence)
  {
    std::vector<std::pair<std::string::size_type, std::string> > idx;
    std::vector<std::string> pattern;
    pattern.push_back("K");
    pattern.push_back("R");
    pattern.push_back("P");

    for (Size i = 0; i < sequence.size(); i++)
    {
      for (Size j = 0; j < pattern.size(); j++)
      {
        if (sequence.substr(i, 1) == pattern[j])
        {
          std::pair<std::string::size_type, std::string> idx_pair(i, pattern[j]);
          idx.push_back(idx_pair);
        }
      }
    }
    return idx;
  }

  float MRMDecoy::AASequenceIdentity(const String& sequence, const String& decoy)
  {
    OPENMS_PRECONDITION(sequence.size() == decoy.size(), "Cannot compare two sequences of unequal length");

    std::vector<char> sequence_v(sequence.begin(), sequence.end());
    std::vector<char> decoy_v(decoy.begin(), decoy.end());
    int running = 0;
    for (Size i = 0; i < sequence_v.size(); i++)
    {
      if (sequence_v[i] == decoy_v[i])
      {
        running += 1;
      }
    }
    double identity = (double) running / sequence_v.size();
    return identity;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::shufflePeptide(
    OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed,
    int max_attempts)
  {
#ifdef DEBUG_MRMDECOY
    std::cout << " shuffle peptide " << peptide.sequence << std::endl;
    seed = 41;
#endif
    if (seed == -1)
    {
      seed = time(0);
    }
    OpenMS::TargetedExperiment::Peptide shuffled = peptide;

    boost::mt19937 generator(seed);
    boost::uniform_int<> uni_dist;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

    typedef std::vector<std::pair<std::string::size_type, std::string> > IndexType; 
    IndexType idx = MRMDecoy::find_all_tryptic(peptide.sequence);
    std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    int aa_size = 17;

    int attempts = 0;
    // loop: copy the original peptide, attempt to shuffle it and check whether difference is large enough
    while (MRMDecoy::AASequenceIdentity(peptide.sequence, shuffled.sequence) > identity_threshold && 
        attempts < max_attempts)
    {
      shuffled = peptide;
      std::vector<Size> peptide_index;
      for (Size i = 0; i < peptide.sequence.size(); i++)
      {
        peptide_index.push_back(i);
      }

      // we erase the indices where K/P/R are (from the back / in reverse order
      // to not delete indices we access later)
      for (IndexType::reverse_iterator it = idx.rbegin(); it != idx.rend(); ++it)
      {
        peptide_index.erase(peptide_index.begin() + it->first);
      }

      // shuffle the peptide index (without the K/P/R which we leave in place)
      std::random_shuffle(peptide_index.begin(), peptide_index.end(), pseudoRNG);

      // re-insert the missing K/P/R at the appropriate places
      for (IndexType::iterator it = idx.begin(); it != idx.end(); ++it)
      {
        peptide_index.insert(peptide_index.begin() + it->first, it->first);
      }

      // use the shuffled index to create the get the new peptide sequence and
      // then to place the modifications at their appropriate places (at the
      // same, shuffled AA where they were before).
      for (Size i = 0; i < peptide_index.size(); i++)
      {
        shuffled.sequence[i] = peptide.sequence[peptide_index[i]];
      }
      for (Size j = 0; j < shuffled.mods.size(); j++)
      {
        for (Size k = 0; k < peptide_index.size(); k++)
        {
          // C and N terminal mods are implicitly not shuffled because they live at positions -1 and sequence.size()
          if (boost::numeric_cast<int>(peptide_index[k]) == shuffled.mods[j].location)
          {
            shuffled.mods[j].location = boost::numeric_cast<int>(k);
            break;
          }
        }
      }

#ifdef DEBUG_MRMDECOY
      for (Size j = 0; j < shuffled.mods.size(); j++)
      {
        std::cout << " position after shuffling " << shuffled.mods[j].location << " mass difference " << shuffled.mods[j].mono_mass_delta << std::endl;
      }
#endif

      ++attempts;

      // If our attempts have failed so far, we will append two random AA to
      // the sequence and see whether we can achieve sufficient shuffling with
      // these additional AA added to the sequence.
      if (attempts % 10 == 9)
      {
        int pos = (pseudoRNG() % aa_size);
        peptide.sequence.append(aa[pos]);
        pos = (pseudoRNG() % aa_size);
        peptide.sequence.append(aa[pos]);

        // now make the shuffled peptide the same length as the new peptide
        shuffled = peptide;
      }
    }

    return shuffled;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::pseudoreversePeptide(
    OpenMS::TargetedExperiment::Peptide peptide)
  {
    OpenMS::TargetedExperiment::Peptide peptideorig = peptide;
    std::vector<Size> peptide_index;
    for (Size i = 0; i < peptide.sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    peptide.sequence = peptide.sequence.substr(0, peptide.sequence.size() - 1).reverse()
                       + peptide.sequence.substr(peptide.sequence.size() - 1, 1); // pseudo-reverse
    std::reverse(peptide_index.begin(), peptide_index.end() - 1);

    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      for (Size k = 0; k < peptide_index.size(); k++)
      {
        if (boost::numeric_cast<int>(peptide_index[k])  == peptide.mods[j].location)
        {
          peptide.mods[j].location = boost::numeric_cast<int>(k);
          break;
        }
      }
    }

    return peptide;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide(
    OpenMS::TargetedExperiment::Peptide peptide)
  {
    OpenMS::TargetedExperiment::Peptide peptideorig = peptide;
    std::vector<Size> peptide_index;
    for (Size i = 0; i < peptide.sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    peptide.sequence = peptide.sequence.reverse();
    std::reverse(peptide_index.begin(), peptide_index.end());

    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      for (Size k = 0; k < peptide_index.size(); k++)
      {
        if (boost::numeric_cast<int>(peptide_index[k]) == peptide.mods[j].location)
        {
          peptide.mods[j].location = boost::numeric_cast<int>(k);
          break;
        }
      }
    }

    return peptide;
  }

  bool MRMDecoy::has_CNterminal_mods(const OpenMS::TargetedExperiment::Peptide & peptide)
  {
    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      if (peptide.mods[j].location == -1 || peptide.mods[j].location == boost::numeric_cast<int>(peptide.sequence.size()) )
      {
        return true;
      }
    }
    return false;
  }

  void MRMDecoy::restrictTransitions(OpenMS::TargetedExperiment& exp, int min_transitions,
                                     int max_transitions)
  {
    OpenMS::TargetedExperiment restricted_exp;
    MRMDecoy::PeptideVectorType peptides;
    MRMDecoy::ProteinVectorType proteins;
    MRMDecoy::TransitionVectorType transitions;

    Map<String, MRMDecoy::TransitionVectorType> TransitionsMap;

    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];
      if (TransitionsMap.find(tr.getPeptideRef()) == TransitionsMap.end())
      {
        TransitionsMap[tr.getPeptideRef()];
      }
      TransitionsMap[tr.getPeptideRef()].push_back(tr);
    }

    for (Map<String, MRMDecoy::TransitionVectorType>::iterator m = TransitionsMap.begin();
         m != TransitionsMap.end(); m++)
    {
      if (m->second.size() >= (Size)min_transitions)
      {
        std::vector<double> LibraryIntensity;
        for (MRMDecoy::TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); tr_it++)
        {
          ReactionMonitoringTransition tr = *tr_it;
          LibraryIntensity.push_back(boost::lexical_cast<double>(tr.getLibraryIntensity()));
        }

        // sort by intensity, reverse and delete all elements after max_transitions 
        sort(LibraryIntensity.begin(), LibraryIntensity.end());
        reverse(LibraryIntensity.begin(), LibraryIntensity.end());
        if ((Size)max_transitions < LibraryIntensity.size() )
        {
          std::vector<double>::iterator start_delete = LibraryIntensity.begin();
          std::advance(start_delete, max_transitions);
          LibraryIntensity.erase(start_delete, LibraryIntensity.end() );
        }

        for (MRMDecoy::TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); tr_it++)
        {
          ReactionMonitoringTransition tr = *tr_it;
          if (std::find( LibraryIntensity.begin(), LibraryIntensity.end(),
                boost::lexical_cast<double>(tr.getLibraryIntensity()) ) != LibraryIntensity.end())
          {
            transitions.push_back(tr);
          }
        }
      }
    }

    std::vector<String> ProteinList;
    for (Size i = 0; i < exp.getPeptides().size(); i++)
    {
      TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      if (TransitionsMap.find(peptide.id) != TransitionsMap.end())
      {
        peptides.push_back(peptide);
        for (Size j = 0; j < peptide.protein_refs.size(); j++)
        {
          ProteinList.push_back(peptide.protein_refs[j]);
        }
      }
    }

    for (Size i = 0; i < exp.getProteins().size(); i++)
    {
      OpenMS::TargetedExperiment::Protein protein = exp.getProteins()[i];
      if (find(ProteinList.begin(), ProteinList.end(), protein.id) != ProteinList.end())
      {
        proteins.push_back(protein);
      }
    }

    restricted_exp.setTransitions(transitions);
    restricted_exp.setPeptides(peptides);
    restricted_exp.setProteins(proteins);

    exp = restricted_exp;
  }

  void MRMDecoy::generateDecoys(OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec,
                                String method, String decoy_tag, double identity_threshold, int max_attempts, 
                                double mz_threshold, bool theoretical, double mz_shift, bool exclude_similar, 
                                double similarity_threshold, bool remove_CNterminal_mods, double precursor_mass_shift,
                                bool enable_losses, bool enable_isotopes)  {
    MRMDecoy::PeptideVectorType peptides;
    MRMDecoy::ProteinVectorType proteins;
    MRMDecoy::TransitionVectorType decoy_transitions;
    for (Size i = 0; i < exp.getProteins().size(); i++)
    {
      OpenMS::TargetedExperiment::Protein protein = exp.getProteins()[i];
      protein.id = decoy_tag + protein.id;
      proteins.push_back(protein);
    }

    std::vector<String> exclusion_peptides;
    // Go through all peptides and apply the decoy method to the sequence
    // (pseudo-reverse, reverse or shuffle). Then set the peptides and proteins of the decoy
    // experiment.
    for (Size i = 0; i < exp.getPeptides().size(); i++)
    {
      OpenMS::TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      // continue if the peptide has C/N terminal modifications and we should exclude them
      if (remove_CNterminal_mods && MRMDecoy::has_CNterminal_mods(peptide)) {continue;}
      peptide.id = decoy_tag + peptide.id;
      OpenMS::String original_sequence = peptide.sequence;

      if (method == "pseudo-reverse")
      {
        peptide = MRMDecoy::pseudoreversePeptide(peptide);
      }
      else if (method == "reverse")
      {
        peptide = MRMDecoy::reversePeptide(peptide);
      }
      else if (method == "shuffle")
      {
        peptide = MRMDecoy::shufflePeptide(peptide, identity_threshold, -1, max_attempts);
      }
      for (Size i = 0; i < peptide.protein_refs.size(); i++)
      {
        peptide.protein_refs[i] = decoy_tag + peptide.protein_refs[i];
      }

      if (MRMDecoy::AASequenceIdentity(original_sequence, peptide.sequence) > identity_threshold)
      {
        if (!exclude_similar)
        { 
          std::cout << "Target sequence: " << original_sequence << " Decoy sequence: " << peptide.sequence  << " Sequence identity: " << MRMDecoy::AASequenceIdentity(original_sequence, peptide.sequence) << " Identity threshold: " << identity_threshold << std::endl;
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "AA Sequences are too similar. Either decrease identity_threshold and increase max_attempts for the shuffle method or set flag exclude_similar.");
        }
        else
        {
          exclusion_peptides.push_back(peptide.id);
        }
      }

      peptides.push_back(peptide);
    }
    dec.setPeptides(peptides);
    dec.setProteins(proteins);

    // hash of the peptide reference containing all transitions
    MRMDecoy::PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      peptide_trans_map[exp.getTransitions()[i].getPeptideRef()].push_back(&exp.getTransitions()[i]);
    }

    Size progress = 0;
    startProgress(0, exp.getTransitions().size(), "Creating decoys");
    for (MRMDecoy::PeptideTransitionMapType::iterator pep_it = peptide_trans_map.begin();
         pep_it != peptide_trans_map.end(); pep_it++)
    {
      String peptide_ref = pep_it->first;
      String decoy_peptide_ref = decoy_tag + pep_it->first; // see above, the decoy peptide id is computed deterministically from the target id
      const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(peptide_ref);
      // continue if the peptide has C/N terminal modifications and we should exclude them
      if (remove_CNterminal_mods && MRMDecoy::has_CNterminal_mods(target_peptide)) {continue;}

      const TargetedExperiment::Peptide decoy_peptide = dec.getPeptideByRef(decoy_peptide_ref);
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
      OpenMS::AASequence decoy_peptide_sequence = TargetedExperimentHelper::getAASequence(decoy_peptide);
      MRMDecoy::IonSeries decoy_ionseries = getIonSeries(decoy_peptide_sequence, decoy_peptide.getChargeState());
      MRMDecoy::IonSeries target_ionseries = getIonSeries(target_peptide_sequence, target_peptide.getChargeState());

      for (Size i = 0; i < pep_it->second.size(); i++)
      {
        setProgress(++progress);
        const ReactionMonitoringTransition tr = *(pep_it->second[i]);
        ReactionMonitoringTransition decoy_tr = tr; // copy the target transition

        decoy_tr.setNativeID(decoy_tag + tr.getNativeID());
        decoy_tr.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
        decoy_tr.setPrecursorMZ(tr.getPrecursorMZ() + precursor_mass_shift); // fix for TOPPView: Duplicate precursor MZ is not displayed.

        // determine the current annotation for the target ion and then select
        // the appropriate decoy ion for this target transition
        std::pair<String, double> targetion = getTargetIon(tr.getProductMZ(), mz_threshold, target_ionseries, enable_losses, enable_isotopes);
        std::pair<String, double> decoyion = getDecoyIon(targetion.first, decoy_ionseries);
        if (method == "shift")
        {
          decoy_tr.setProductMZ(decoyion.second + mz_shift);
        }
        else
        {
          decoy_tr.setProductMZ(decoyion.second);
        }
        decoy_tr.setPeptideRef(decoy_tag + tr.getPeptideRef());

        if (decoyion.second > 0)
        {
          if (similarity_threshold >=0)
          {
            if (std::fabs(tr.getProductMZ() - decoy_tr.getProductMZ()) < similarity_threshold)
            {
              exclusion_peptides.push_back(decoy_tr.getPeptideRef());
            }
          }
         decoy_transitions.push_back(decoy_tr);
        }
      } // end loop over transitions
    } // end loop over peptides
    endProgress();

    if (exclude_similar)
    {
      MRMDecoy::TransitionVectorType filtered_decoy_transitions;
      for ( MRMDecoy::TransitionVectorType::iterator tr_it = decoy_transitions.begin(); tr_it != decoy_transitions.end(); tr_it++)
      {
        if (std::find(exclusion_peptides.begin(), exclusion_peptides.end(), tr_it->getPeptideRef())==exclusion_peptides.end())
        {
          filtered_decoy_transitions.push_back(*tr_it);
        }
        else
        {
          std::cout << "Excluded: " << tr_it->getPeptideRef() << std::endl;
        }
      }
      dec.setTransitions(filtered_decoy_transitions);
    }
    else 
    {
      dec.setTransitions(decoy_transitions);
    }

    if (theoretical)
    {
      correctMasses(exp, mz_threshold, enable_losses, enable_isotopes);
    }
  }

  void MRMDecoy::correctMasses(OpenMS::TargetedExperiment& exp, double mz_threshold, bool enable_losses, bool enable_isotopes)
  {
    MRMDecoy::TransitionVectorType target_transitions;
    // hash of the peptide reference containing all transitions
    MRMDecoy::PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      peptide_trans_map[exp.getTransitions()[i].getPeptideRef()].push_back(&exp.getTransitions()[i]);
    }

    {
      Size progress = 0;
      startProgress(0, exp.getTransitions().size(), "Correcting masses (theoretical)");
      for (MRMDecoy::PeptideTransitionMapType::iterator pep_it = peptide_trans_map.begin();
           pep_it != peptide_trans_map.end(); pep_it++)
      {
        const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(pep_it->first);
        OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
        MRMDecoy::IonSeries target_ionseries = getIonSeries(target_peptide_sequence, target_peptide.getChargeState());

        for (Size i = 0; i < pep_it->second.size(); i++)
        {
          setProgress(++progress);
          const ReactionMonitoringTransition tr = *(pep_it->second[i]);

          // determine the current annotation for the target ion 
          std::pair<String, double> targetion = getTargetIon(tr.getProductMZ(), mz_threshold, target_ionseries, enable_losses, enable_isotopes);
          // correct the masses of the input experiment 
          {
            ReactionMonitoringTransition transition = *(pep_it->second[i]); // copy the transition
            if (targetion.second > 0)
            {
              transition.setProductMZ(targetion.second);
              target_transitions.push_back(transition);
            }
          }
        } // end loop over transitions
      } // end loop over peptides
      endProgress();
      exp.setTransitions(target_transitions);
    }
  }

}
