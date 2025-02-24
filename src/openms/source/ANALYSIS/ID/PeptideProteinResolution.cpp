// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <queue>
#include <unordered_set>
#include <algorithm>

using namespace OpenMS;
using namespace std;


namespace OpenMS
{
  std::ostream& operator<<(std::ostream& os, const ConnectedComponent& conn_comp)
  {
    os << "Proteins: ";
    for (std::set<Size>::const_iterator prot_it = conn_comp.prot_grp_indices.begin();
          prot_it != conn_comp.prot_grp_indices.end();
          ++prot_it)
    {
      os << *prot_it << ",";
    }
    os << std::endl;
    os << "Peptides: ";
    for (std::set<Size>::const_iterator pep_it = conn_comp.pep_indices.begin();
          pep_it != conn_comp.pep_indices.end();
          ++pep_it)
    {
      os << *pep_it << ",";
    }

    return os;    
  }


  // C'tor
  PeptideProteinResolution::PeptideProteinResolution(bool statistics) :
      /*indist_prot_grp_td_(),*/ 
      indist_prot_grp_to_pep_(), 
      pep_to_indist_prot_grp_(),
      prot_acc_to_indist_prot_grp_(), 
      statistics_(statistics)
  {
  }

  void PeptideProteinResolution::resolve(ProteinIdentification& protein,
                                            vector<PeptideIdentification>& peptides,
                                            bool resolve_ties,
                                            bool targets_first)
  {

    vector<ProteinIdentification::ProteinGroup>& groups = protein.getIndistinguishableProteins();
    vector<bool> indist_prot_grp_decoy(groups.size());
    unordered_map<string, Size> prot_acc_to_indist_prot_grp;

    if (groups.empty())
    {
      throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "No indistinguishable Groups annotated. Currently this class only resolves across groups.");
    }

    OPENMS_LOG_INFO << "Resolving peptides between " << protein.getHits().size() << " proteins in " << groups.size() << " indistinguishable groups." << std::endl;

    // I don't think we need to assume sortedness here
    //if (!skip_sort) sort(groups.begin(), groups.end());

    std::unordered_set<std::string> decoy_accs;
    for (const ProteinHit& p : protein.getHits())
    {
      if (p.metaValueExists("target_decoy") && p.getMetaValue("target_decoy") == "decoy")
      {
        decoy_accs.insert(p.getAccession());
      }
    }

    // Construct intermediate mapping of single protein accessions
    // to indist. protein groups
    for (vector<ProteinIdentification::ProteinGroup>::const_iterator group_it =
        groups.begin();
         group_it != groups.end(); ++group_it)
    {
      for (vector<String>::const_iterator acc_it = group_it->accessions.begin();
           acc_it != group_it->accessions.end(); ++acc_it)
      {
        Size idx = group_it - groups.begin();
        prot_acc_to_indist_prot_grp[*acc_it] = idx;
        if (decoy_accs.find(*acc_it) != decoy_accs.end())
        {
          indist_prot_grp_decoy[idx] = true;
        }
      }
    }

    // Go through PeptideIDs
    for (PeptideIdentification& pep : peptides)
    {

      vector<PeptideHit>& hits = pep.getHits();
      if (!hits.empty())
      {
        PeptideHit& best_hit = hits[0];
        const vector<PeptideEvidence>& pepev = best_hit.getPeptideEvidences();
        set<Size> bestNonDecoyGrpTie;
        set<Size> bestDecoyGrpTie;
        unordered_map<Size,set<Size>> grpIdxToEvIdx;

        Size ev_idx = 0;
        for (vector<PeptideEvidence>::const_iterator pepev_it = pepev.begin();
             pepev_it != pepev.end(); ++pepev_it, ++ev_idx)
        {
          String acc = pepev_it->getProteinAccession();

          auto found = prot_acc_to_indist_prot_grp.find(acc);
          if (found == prot_acc_to_indist_prot_grp.end())
          {
            throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Not all proteins present in an indistinguishable group. Make sure to add them as singletons.");
          }
          else
          {
            Size prot_group_index = found->second;
            auto it = grpIdxToEvIdx.emplace(prot_group_index, set<Size>());
            it.first->second.insert(ev_idx);
            //TODO work with a tolerance for doubles instead?
            if (indist_prot_grp_decoy[prot_group_index])
            {
              if (bestDecoyGrpTie.empty() ||
                  groups[prot_group_index].probability < groups[*bestDecoyGrpTie.begin()].probability)
              {
                bestDecoyGrpTie.clear();
                bestDecoyGrpTie.insert(prot_group_index);
              }
              else if (groups[prot_group_index].probability == groups[*bestDecoyGrpTie.begin()].probability)
              {
                bestDecoyGrpTie.insert(prot_group_index);
              }
            }
            else
            {
              if (bestNonDecoyGrpTie.empty() ||
                  groups[prot_group_index].probability < groups[*bestNonDecoyGrpTie.begin()].probability)
              {
                bestNonDecoyGrpTie.clear();
                bestNonDecoyGrpTie.insert(prot_group_index);
              }
              else if (groups[prot_group_index].probability == groups[*bestNonDecoyGrpTie.begin()].probability)
              {
                bestNonDecoyGrpTie.insert(prot_group_index);
              }
            }
          }
        }

        bool targets_first_resolve_ties = false;

        set<Size>* toResolve;
        set<Size> allGrpsSet;
        if (bestNonDecoyGrpTie.empty())
        {
          toResolve = &bestDecoyGrpTie;
        }
        else if (bestDecoyGrpTie.empty() || targets_first)
        {
          toResolve = &bestNonDecoyGrpTie;
        }
        else if (groups[*bestNonDecoyGrpTie.begin()].probability > groups[*bestDecoyGrpTie.begin()].probability)
        {
          toResolve = &bestNonDecoyGrpTie;
        }
        else if (groups[*bestDecoyGrpTie.begin()].probability > groups[*bestNonDecoyGrpTie.begin()].probability)
        {
          toResolve = &bestDecoyGrpTie;
        }
        else // both equal
        {
          if (resolve_ties && targets_first_resolve_ties)
          {
            toResolve = &bestNonDecoyGrpTie;
          }
          else // take all best groups
          {
            vector<Size> allGrps;
            merge(std::begin(bestNonDecoyGrpTie), std::end(bestNonDecoyGrpTie),
                  std::begin(bestDecoyGrpTie), std::end(bestDecoyGrpTie),
                  std::back_inserter(allGrps));
            allGrpsSet = set<Size>(allGrps.begin(),allGrps.end());
            toResolve = &allGrpsSet;
          }
        }

        set<Size> evToKeep;
        if (resolve_ties)
        {
          //TODO this tie resolution basically just takes the first group that occurred
          evToKeep = grpIdxToEvIdx[*toResolve->begin()];
          if (toResolve->size() > 1)
          {
           OPENMS_LOG_INFO << "Resolution: Peptide " << pep.getHits()[0].getSequence().toString() << " had groups:" << std::endl;

           OPENMS_LOG_INFO << "tgt: ";
            for (const auto& g : bestNonDecoyGrpTie)
            {
              OPENMS_LOG_INFO << g << "=" << groups[g].probability << ", ";
            }
           OPENMS_LOG_INFO << std::endl;
           OPENMS_LOG_INFO << "dec: ";
            for (const auto& g : bestDecoyGrpTie)
            {
              OPENMS_LOG_INFO << g << "=" << groups[g].probability << ", ";
            }
           OPENMS_LOG_INFO << std::endl;
           OPENMS_LOG_INFO << "Kept: " << *toResolve->begin() << std::endl;
          }
        }
        else
        {
          for (const auto& grp : *toResolve)
          {
            evToKeep.insert(grpIdxToEvIdx[grp].begin(),grpIdxToEvIdx[grp].end());
          }
        }

        vector<PeptideEvidence> newEv;
        newEv.reserve(evToKeep.size());
        for (const auto& idx : evToKeep)
        {
          newEv.push_back(pepev[idx]);
        }
        best_hit.setPeptideEvidences(newEv);
      }
      else
      {
       OPENMS_LOG_WARN << "Warning PeptideProteinResolution: Skipping spectrum without hits." << std::endl;
      }
    }
  }


  // Initialization of global variables (= graph)
  void PeptideProteinResolution::buildGraph(ProteinIdentification& protein,
                      const vector<PeptideIdentification>& peptides, bool skip_sort)
  {
    vector<ProteinIdentification::ProteinGroup>& groups = protein.getIndistinguishableProteins();

    if (groups.empty())
    {
      throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "No indistinguishable Groups annotated. Currently this class only resolves across groups.");
    }

   OPENMS_LOG_INFO << "Resolving peptides between " << protein.getHits().size() << " proteins in " << groups.size() << " indistinguishable groups." << std::endl;


    if (!skip_sort) sort(groups.begin(), groups.end());

    // TODO this is only needed for target_first option
    std::unordered_set<std::string> decoy_accs;
    for (const ProteinHit& p : protein.getHits())
    {
      if (p.metaValueExists("target_decoy") && p.getMetaValue("target_decoy") == "decoy")
      {
        decoy_accs.insert(p.getAccession());
      }
    }

    // Construct intermediate mapping of single protein accessions
    // to indist. protein groups
    for (vector<ProteinIdentification::ProteinGroup>::const_iterator group_it =
         groups.begin();
         group_it != groups.end(); ++group_it)
    {
      for (vector<String>::const_iterator acc_it = group_it->accessions.begin();
           acc_it != group_it->accessions.end(); ++acc_it)
      {
        Size idx = group_it - groups.begin();
        prot_acc_to_indist_prot_grp_[*acc_it] = idx;
      }
    }
    
    // Go through PeptideIDs and construct a bidirectional mapping
    for (vector<PeptideIdentification>::const_iterator pep_it = peptides.begin();
         pep_it != peptides.end();
         ++pep_it)
    {
      Size pep_index = pep_it - peptides.begin();

      const vector<PeptideHit>& hits = pep_it->getHits();
      if (!hits.empty())
      {
        PeptideHit best_hit = hits[0];
        const vector<PeptideEvidence>& pepev = best_hit.getPeptideEvidences();

        for (vector<PeptideEvidence>::const_iterator pepev_it = pepev.begin();
             pepev_it != pepev.end(); ++pepev_it)
        {
          String acc = pepev_it->getProteinAccession();
          auto found = prot_acc_to_indist_prot_grp_.find(acc);
          if (found == prot_acc_to_indist_prot_grp_.end())
          {
            throw Exception::MissingInformation(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Not all proteins present in an indistinguishable group. (" + acc + " not found). Make sure to add them as singletons.");
          }
          else
          {
            Size prot_group_index = found->second;
            pep_to_indist_prot_grp_[pep_index].insert(prot_group_index);
            indist_prot_grp_to_pep_[prot_group_index];
            indist_prot_grp_to_pep_[prot_group_index].insert(pep_index);
          }
        }
      }
      else
      {
       OPENMS_LOG_WARN << "Warning PeptideProteinResolution: Skipping spectrum without hits." << std::endl;
      }
    }
  }

  // "Main" function
  void PeptideProteinResolution::resolveGraph(ProteinIdentification& protein,
                  vector<PeptideIdentification>& peptides)
  {
    //Debugging
    Size old_size = indist_prot_grp_to_pep_.size();
    
    //Statistics
    ConnectedComponent most_peps;
    ConnectedComponent most_grps;
    ConnectedComponent most_both;
    
    // Traverse every connected component, remove visited "nodes" in each step
    while (!indist_prot_grp_to_pep_.empty())
    {
      if (statistics_ && (old_size - indist_prot_grp_to_pep_.size() > 1))
      {
        OPENMS_LOG_INFO << "resolved group of size "
        << old_size - indist_prot_grp_to_pep_.size() << " in last step "
        << endl;
        old_size = indist_prot_grp_to_pep_.size();
      }
      
      // We take any (= first) protein from map that is still left,
      // to start the next BFS from it
      Size root_prot_grp = indist_prot_grp_to_pep_.begin()->first;
      
      // do BFS, return connected proteins and peptides
      ConnectedComponent curr_component =
      PeptideProteinResolution::findConnectedComponent(root_prot_grp);
      // For debugging and statistics
      if (statistics_)
      {
        if (curr_component.prot_grp_indices.size() >
            most_grps.prot_grp_indices.size())
        {
          most_grps = curr_component;
        }
        
        if (curr_component.pep_indices.size() >
            most_peps.pep_indices.size())
        {
          most_peps = curr_component;
        }
        
        if ((curr_component.prot_grp_indices.size() +
             curr_component.pep_indices.size()) >
            (most_both.prot_grp_indices.size() +
             most_both.pep_indices.size()))
        {
          most_both = curr_component;
        }
        
        if (curr_component.prot_grp_indices.size() > 1)
        {
          OPENMS_LOG_INFO << "found group: " << endl;
          OPENMS_LOG_INFO << curr_component;
          OPENMS_LOG_INFO << endl << "Processing ..." << endl;
        }
      }
      
      // resolve shared peptides based on posterior probabilities
      // -> modifies PeptideIDs in peptides
      PeptideProteinResolution::resolveConnectedComponent(curr_component,
                                                           protein,
                                                           peptides);
      
      // mark proteins of this component as visited by removing them
      for (set<Size>::iterator grp_it =
           curr_component.prot_grp_indices.begin();
           grp_it != curr_component.prot_grp_indices.end();
           ++grp_it)
      {
        indist_prot_grp_to_pep_.erase(*grp_it);
      }
    }
    
    //TODO maybe extend statistics of connected components!
    if (statistics_)
    {
      OPENMS_LOG_INFO << endl << "Most protein groups in component:" << endl;
      OPENMS_LOG_INFO << most_grps;
      OPENMS_LOG_INFO << endl << "Most peptides in component:"<< endl;
      OPENMS_LOG_INFO << most_peps;
      OPENMS_LOG_INFO << endl << "Biggest component:" << endl;
      OPENMS_LOG_INFO << most_both;
    }
  }

  /*
   * Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
   * and peptides), switching from one to the other in each step.
   * Returns a Connected Component as set of group and peptide indices.
   */
  ConnectedComponent PeptideProteinResolution::findConnectedComponent(Size& root_prot_grp)
  {
    // init result
    ConnectedComponent conn_comp;
  
    // init queue, bool keeps track of if we need to use
    // ProteinGroup -> Peptide (true) or
    // Peptide -> ProteinGroup (false) as mapping
    queue<pair<bool, Size> > my_queue;
  
    // start with given root
    my_queue.push(make_pair(true, root_prot_grp));
  
    // check successes of insertions
    pair<set<Size>::iterator,bool> success;
  
    while (!my_queue.empty())
    {
      // save first element and pop
      pair<bool, Size> curr_node = my_queue.front();
      my_queue.pop();
    
      // initialize neighbors
      set<Size> neighbors;
    
      // Choose correct map, depending on if we deal with protGrp or peptide
      if (curr_node.first)
      {
        neighbors = indist_prot_grp_to_pep_[curr_node.second];
      }
      else
      {
        neighbors = pep_to_indist_prot_grp_[curr_node.second];
      }
    
      for (set<Size>::iterator nb_it = neighbors.begin();
         nb_it != neighbors.end();
         ++nb_it)
      {
        // If current node is protein, its neighbors are peptides and
        // vice versa -> look in corresponding "result" set and insert
        // if not present
        if (!curr_node.first)
        {
          success = conn_comp.prot_grp_indices.insert(*nb_it);
        }
        else
        {
          success = conn_comp.pep_indices.insert(*nb_it);
        }
      
        // If it was not seen yet, add it to the queue to process
        // its neighbors later. All neighbors are from the opposite type now
        if (success.second)
        {
          my_queue.push(make_pair(!curr_node.first, *nb_it));
        }
      }
    }
  
    return conn_comp;
  }

  /* TODO this does not produce correct results yet. Check again.
   * Resolves connected components based on Fido probabilities and adds them
   * as additional protein_groups to the output idXML.
   * Thereby greedily assigns shared peptides in this component uniquely to
   * the proteins of the current BEST INDISTINGUISHABLE protein group,
   * ready to be used in ProteinQuantifier then.
   * This is achieved by removing all other evidence from the input
   * PeptideIDs and iterating until
   * In accordance with Fido only the best hit (PSM) for an ID is considered.
   * Probability ties are _currently_ resolved by taking the first occurrence.
   */
/*  void PeptideProteinResolution::resolveConnectedComponentTargetsFirst(
      ConnectedComponent& conn_comp,
      ProteinIdentification& protein,
      vector<PeptideIdentification>& peptides,
      bool targets_first)
  {
    // Nothing to resolve in a singleton group (will not be added to output though)
    if (conn_comp.prot_grp_indices.size() == 1) return;

    // Add proteins from a connected component to ambiguity groups
    ProteinIdentification::ProteinGroup ambiguity_grp;

    vector<ProteinIdentification::ProteinGroup>& origin_groups = protein.getIndistinguishableProteins();

    // Save the max probability in this component to add it (should be first one, since groups were sorted and
    // lower index means higher score and set is sorted by index)
    ambiguity_grp.probability = origin_groups[*conn_comp.prot_grp_indices.begin()].probability;

    for (set<Size>::iterator grp_it = conn_comp.prot_grp_indices.begin();
         grp_it != conn_comp.prot_grp_indices.end();
         ++grp_it)
    {
      if (*grp_it >= origin_groups.size())
      {
       OPENMS_LOG_FATAL_ERROR << "Something went terribly wrong. "
                           "Group with index " << *grp_it << "doesn't exist. "
                                                             " ProteinPeptideResolution: Groups changed"
                                                             " after building data structures." << std::endl;
      }

      vector<String> accessions = origin_groups[*grp_it].accessions;

      // Put the accessions of the indist. groups into the subsuming
      // ambiguity group
      ambiguity_grp.accessions.insert(ambiguity_grp.accessions.end(),
                                      accessions.begin(),
                                      accessions.end());

      if (targets_first && indist_prot_grp_td_[*grp_it].first)
      {
        if (statistics_)
        {
         OPENMS_LOG_DEBUG << "Group: ";
          for (const String& s : origin_groups[*grp_it].accessions)
          {
           OPENMS_LOG_DEBUG << s << ", ";
          }
         OPENMS_LOG_DEBUG << " steals " << indist_prot_grp_to_pep_[*grp_it].size() << " peptides for itself." << std::endl;
        }
        // Update all the peptides the current best point to
        for (set<Size>::iterator pepid_it =
            indist_prot_grp_to_pep_[*grp_it].begin();
             pepid_it != indist_prot_grp_to_pep_[*grp_it].end(); ++pepid_it)
        {
          vector<PeptideHit> pep_id_hits = peptides[*pepid_it].getHits();
          vector<PeptideEvidence> best_hit_ev =
              pep_id_hits[0].getPeptideEvidences();

          // go through all the evidence of this peptide and remove all
          // proteins but the ones from the current indist. group
          for (vector<PeptideEvidence>::iterator pepev_it = best_hit_ev.begin();
               pepev_it != best_hit_ev.end();
            //don't increase index, will be done by case
              )
          {
            // if its accession is not in the current best group, remove evidence
            if (find(accessions.begin(),
                     accessions.end(),
                     pepev_it->getProteinAccession()) == accessions.end())
            {
              // we get valid iterator from erase with shifted objects
              pepev_it = best_hit_ev.erase(pepev_it);
              // also erase from the mapping of this class
              indist_prot_grp_to_pep_[prot_acc_to_indist_prot_grp_[pepev_it->getProteinAccession()]].erase(*pepid_it);
            }
            else
            { // iterate further
              ++pepev_it;
            }
          }
          // Set the remaining evidences as new evidence
          pep_id_hits[0].setPeptideEvidences(best_hit_ev);
          peptides[*pepid_it].setHits(pep_id_hits);
        }
      }

      if (targets_first) // we need a second run with only decoys to resolve potential remaining peptides
      {
        for (set<Size>::iterator grp_it = conn_comp.prot_grp_indices.begin();
             grp_it != conn_comp.prot_grp_indices.end();
             ++grp_it)
        {
          if (*grp_it >= origin_groups.size())
          {
           OPENMS_LOG_FATAL_ERROR << "Something went terribly wrong. "
                               "Group with index " << *grp_it << "doesn't exist. "
                                                                 " ProteinPeptideResolution: Groups changed"
                                                                 " after building data structures." << std::endl;
          }

          vector<String> accessions = origin_groups[*grp_it].accessions;

          // Put the accessions of the indist. groups into the subsuming
          // ambiguity group
          ambiguity_grp.accessions.insert(ambiguity_grp.accessions.end(),
                                          accessions.begin(),
                                          accessions.end());

          if (!indist_prot_grp_td_[*grp_it].first)
          {
            if (statistics_)
            {
             OPENMS_LOG_DEBUG << "Group: ";
              for (const String& s : origin_groups[*grp_it].accessions)
              {
               OPENMS_LOG_DEBUG << s << ", ";
              }
             OPENMS_LOG_DEBUG << " steals " << indist_prot_grp_to_pep_[*grp_it].size() << " peptides for itself." << std::endl;
            }

            // Update all the peptides the current best point to
            for (set<Size>::iterator pepid_it =
                indist_prot_grp_to_pep_[*grp_it].begin();
                 pepid_it != indist_prot_grp_to_pep_[*grp_it].end(); ++pepid_it)
            {
              vector<PeptideHit> pep_id_hits = peptides[*pepid_it].getHits();
              vector<PeptideEvidence> best_hit_ev =
                  pep_id_hits[0].getPeptideEvidences();

              // go through all the evidence of this peptide and remove all
              // proteins but the ones from the current indist. group
              for (vector<PeptideEvidence>::iterator pepev_it = best_hit_ev.begin();
                   pepev_it != best_hit_ev.end();
                //don't increase index, will be done by case
                  )
              {
                // if its accession is not in the current best group, remove evidence
                if (find(accessions.begin(),
                         accessions.end(),
                         pepev_it->getProteinAccession()) == accessions.end())
                {
                  // we get valid iterator from erase with shifted objects
                  pepev_it = best_hit_ev.erase(pepev_it);
                  // also erase from the mapping of this class
                  indist_prot_grp_to_pep_[prot_acc_to_indist_prot_grp_[pepev_it->getProteinAccession()]].erase(*pepid_it);
                }
                else
                { // iterate further
                  ++pepev_it;
                }
              }
              // Set the remaining evidences as new evidence
              pep_id_hits[0].setPeptideEvidences(best_hit_ev);
              peptides[*pepid_it].setHits(pep_id_hits);
            }
          }
        }
      }
    }

    //Finally insert ambiguity group
    protein.insertProteinGroup(ambiguity_grp);
  }*/


  void PeptideProteinResolution::resolveConnectedComponent(
      ConnectedComponent& conn_comp,
      ProteinIdentification& protein,
      vector<PeptideIdentification>& peptides)
  {
    // TODO think about ignoring decoy proteins (at least when resolving ties!)

    // Nothing to resolve in a singleton group (will not be added to output though)
    if (conn_comp.prot_grp_indices.size() <= 1) return;

    // Add proteins from a connected component to ambiguity groups
    ProteinIdentification::ProteinGroup ambiguity_grp;

    vector<ProteinIdentification::ProteinGroup>& origin_groups = protein.getIndistinguishableProteins();

    // Save the max probability in this component to add it (should be first one, since groups were sorted and
    // lower index means higher score and set is sorted by index)
    size_t best_grp_index = *conn_comp.prot_grp_indices.begin();
    ambiguity_grp.probability = origin_groups[best_grp_index].probability;
    
    // copy group indices so we can reorder them for tie resolution
    vector<Size> prot_grp_indices(conn_comp.prot_grp_indices.begin(), conn_comp.prot_grp_indices.end());

    // groups are currently only sorted by probability.
    // in the presence of ties we need to resolve them by the number of peptides.
    std::sort(prot_grp_indices.begin(), prot_grp_indices.end(), 
      [&](const Size & a, const Size & b) -> bool
      { 
        size_t as = indist_prot_grp_to_pep_[a].size();
        size_t bs = indist_prot_grp_to_pep_[b].size();
        return std::tie(origin_groups[a].probability, as) > std::tie(origin_groups[b].probability, bs);
      });   

    for (vector<Size>::iterator grp_it = prot_grp_indices.begin();
         grp_it != prot_grp_indices.end();
         ++grp_it)
    {
      if (*grp_it >= origin_groups.size())
      {
       OPENMS_LOG_FATAL_ERROR << "Something went terribly wrong. "
                              << "Group with index " << *grp_it << "doesn't exist. "
                              << " ProteinPeptideResolution: Groups changed"
                              << " after building data structures." << std::endl;
      }

      const vector<String>& accessions = origin_groups[*grp_it].accessions;

      // Put the accessions of the indist. groups into the subsuming
      // ambiguity group
      ambiguity_grp.accessions.insert(ambiguity_grp.accessions.end(),
                                      accessions.begin(),
                                      accessions.end());
      if (statistics_)
      {
        OPENMS_LOG_DEBUG << "Group: ";
        for (const String& s : accessions)
        {
          OPENMS_LOG_DEBUG << s << ", ";
        }
        OPENMS_LOG_DEBUG << " steals " << indist_prot_grp_to_pep_[*grp_it].size() << " peptides for itself." << std::endl;
      }

      // Update all the peptides the current best point to
      for (set<Size>::iterator pepid_it = indist_prot_grp_to_pep_[*grp_it].begin();
           pepid_it != indist_prot_grp_to_pep_[*grp_it].end(); ++pepid_it)
      {
        vector<PeptideHit> pep_id_hits = peptides[*pepid_it].getHits();
        vector<PeptideEvidence> best_hit_ev =
            pep_id_hits[0].getPeptideEvidences();

        // Go through all _remaining_ proteins of the component and remove this
        // peptide from their mapping
        vector<Size>::iterator grp_it_cont = grp_it;
        ++grp_it_cont;
        for (; grp_it_cont != prot_grp_indices.end();
                              ++grp_it_cont)
        {
          indist_prot_grp_to_pep_[*grp_it_cont].erase(*pepid_it);
        }

        // go through all the evidence of this peptide and remove all
        // proteins but the ones from the current indist. group
        for (vector<PeptideEvidence>::iterator pepev_it = best_hit_ev.begin();
             pepev_it != best_hit_ev.end();
          //don't increase index, will be done by case
            )
        {
          // if its accession is not in the current best group, remove evidence
          if (find(accessions.begin(),
                   accessions.end(),
                   pepev_it->getProteinAccession()) == accessions.end())
          {
            // we get valid iterator from erase with shifted objects
            pepev_it = best_hit_ev.erase(pepev_it);
          }
          else
          { // iterate further
            ++pepev_it;
          }
        }
        // Set the remaining evidences as new evidence
        pep_id_hits[0].setPeptideEvidences(best_hit_ev);
        peptides[*pepid_it].setHits(pep_id_hits);
      }
    }

    //Finally insert ambiguity group
    protein.insertProteinGroup(ambiguity_grp);
  }


  void PeptideProteinResolution::run(vector<ProteinIdentification>& inferred_protein_ids, 
    vector<PeptideIdentification>& inferred_peptide_ids)
  {
    PeptideProteinResolution ppr;
    ppr.buildGraph(inferred_protein_ids[0], inferred_peptide_ids);
    ppr.resolveGraph(inferred_protein_ids[0], inferred_peptide_ids);    
    IDFilter::removeUnreferencedProteins(inferred_protein_ids, inferred_peptide_ids);
    IDFilter::updateProteinGroups(inferred_protein_ids[0].getIndistinguishableProteins(), inferred_protein_ids[0].getHits());
    IDFilter::updateProteinGroups(inferred_protein_ids[0].getProteinGroups(), inferred_protein_ids[0].getHits());
  }

}

