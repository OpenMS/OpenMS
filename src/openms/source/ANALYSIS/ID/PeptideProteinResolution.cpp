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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>

#include <queue>
#include <ostream>

using namespace OpenMS;
using namespace std;


ConnectedComponent::ConnectedComponent() : prot_grp_indices(std::set<Size>()),
    pep_indices(std::set<Size>())
    {}

std::ostream& OpenMS::operator << (std::ostream& os, const ConnectedComponent& conn_comp)
    {
      os << "Proteins: ";
      for (std::set<Size>::iterator prot_it = conn_comp.prot_grp_indices.begin();
           prot_it != conn_comp.prot_grp_indices.end();
           ++prot_it)
      {
        os << *prot_it << ",";
      }
      os << std::endl;
      os << "Peptides: ";
      for (std::set<Size>::iterator pep_it = conn_comp.pep_indices.begin();
           pep_it != conn_comp.pep_indices.end();
           ++pep_it)
      {
        os << *pep_it << ",";
      }
      
      return os;
      
    }


// C'tor
PeptideProteinResolution::PeptideProteinResolution(bool statistics) :
indist_prot_grp_to_pep_(), pep_to_indist_prot_grp_(),
prot_acc_to_indist_prot_grp_(), statistics_(statistics)
{}

// Initialization of global variables (= graph)
void PeptideProteinResolution::buildGraph(const ProteinIdentification& protein,
                      const vector<PeptideIdentification>& peptides)
{
    // Construct intermediate mapping of single protein accessions
    // to indist. protein groups
    for (vector<ProteinIdentification::ProteinGroup>::const_iterator group_it =
         protein.getIndistinguishableProteins().begin();
         group_it != protein.getIndistinguishableProteins().end(); ++group_it)
    {
      for (vector<String>::const_iterator acc_it = group_it->accessions.begin();
           acc_it != group_it->accessions.end(); ++acc_it)
      {
        prot_acc_to_indist_prot_grp_[*acc_it] =
        group_it - protein.getIndistinguishableProteins().begin();
      }
    }
    
    // Go through PeptideIDs and construct a bidirectional mapping
    for (vector<PeptideIdentification>::const_iterator pep_it = peptides.begin();
         pep_it != peptides.end();
         ++pep_it)
    {
      Size pep_index = pep_it - peptides.begin();
      
      PeptideHit best_hit = pep_it->getHits()[0];
      const vector<PeptideEvidence> pepev = best_hit.getPeptideEvidences();
      
      for (vector<PeptideEvidence>::const_iterator pepev_it = pepev.begin();
           pepev_it != pepev.end(); ++pepev_it)
      {
        String acc = pepev_it->getProteinAccession();
        Size prot_group_index = prot_acc_to_indist_prot_grp_[acc];
        pep_to_indist_prot_grp_[pep_index].insert(prot_group_index);
        indist_prot_grp_to_pep_[prot_group_index];
        indist_prot_grp_to_pep_[prot_group_index].insert(pep_index);
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
        LOG_INFO << "resolved group of size "
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
          LOG_INFO << "found group: " << endl;
          LOG_INFO << curr_component;
          LOG_INFO << endl << "Processing ..." << endl;
        }
      }
      
      // resolve shared peptides based on Fido probabilities
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
      LOG_INFO << endl << "Most protein groups in component:" << endl;
      LOG_INFO << most_grps;
      LOG_INFO << endl << "Most peptides in component:"<< endl;
      LOG_INFO << most_peps;
      LOG_INFO << endl << "Biggest component:" << endl;
      LOG_INFO << most_both;
    }
}


/*
 * Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
 * and peptides), switching from one to the other in each step.
 * Returns a Connected Component as set of group and peptide indices.
 */
ConnectedComponent
PeptideProteinResolution::findConnectedComponent(Size& root_prot_grp)
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

/*
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
void PeptideProteinResolution::resolveConnectedComponent(
                                ConnectedComponent& conn_comp,
                                ProteinIdentification& protein,
                                vector<PeptideIdentification>& peptides)
{
  /* TODO for resolving ties:
   while grpit not at end
   while grpit.probability does not change:
   save index with max nr of peptides in mapping and compare
   grpit++
   resolve(with max index grp)
   */
  
  // Add proteins from a connected component to ambiguity groups
  ProteinIdentification::ProteinGroup ambiguity_grp;
  
  // Save the max probability in this component to add it (should be first one)
  double max_prob = 0.0;
  
  // Go through protein groups (sorted by probability -> higher index
  // means worse probability)
  bool first_change = true;
  
  for (set<Size>::iterator grp_it = conn_comp.prot_grp_indices.begin();
       grp_it != conn_comp.prot_grp_indices.end();
       ++grp_it)
  {
    
    // Take first probability -> best
    if (first_change)
    {
      max_prob = protein.getIndistinguishableProteins()[*grp_it].probability;
      first_change = false;
    }
    
    ambiguity_grp.probability = max_prob;
    
    vector<String> accessions =
    protein.getIndistinguishableProteins()[*grp_it].accessions;
    
    // Put the accessions of the indist. groups into the subsuming
    // ambiguity group
    ambiguity_grp.accessions.insert(ambiguity_grp.accessions.end(),
                                    accessions.begin(),
                                    accessions.end());
    
    // Update all the peptides the current best point to
    for (set<Size>::iterator pepid_it =
         indist_prot_grp_to_pep_[*grp_it].begin();
         pepid_it != indist_prot_grp_to_pep_[*grp_it].end(); ++pepid_it)
    {
      
      vector<PeptideHit> pep_id_hits = peptides[*pepid_it].getHits();
      vector<PeptideEvidence> best_hit_ev =
      pep_id_hits[0].getPeptideEvidences();
      
      
      // Go through all _remaining_ proteins of the component and remove this
      // peptide from their mapping
      set<Size>::iterator grp_it_cont = grp_it;
      ++grp_it_cont;
      for (/*grp_it_cont*/; grp_it_cont != conn_comp.prot_grp_indices.end();
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
