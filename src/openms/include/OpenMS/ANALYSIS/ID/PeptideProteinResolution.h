// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <vector>
#include <set>

namespace OpenMS
{

  /// Represents a connected component of the bipartite graph
  /// Holds indices of peptides and (indist.) protein groups in them
  struct OPENMS_DLLAPI ConnectedComponent
  {
    std::set<Size> prot_grp_indices;
    std::set<Size> pep_indices;
    
    // Default constructor
    ConnectedComponent();

    /// Overloaded operator '<<' for ConnectedComponents
    friend std::ostream& operator << (std::ostream& os, const ConnectedComponent& conn_comp);
  };

  /**
   @brief Resolves shared peptides based on protein scores
   
   Resolves connected components of the bipartite protein-peptide graph based
   on protein probabilities/scores and adds them as additional protein_groups
   to the protein identification run processed.
   Thereby greedily assigns shared peptides in this component uniquely to the
   proteins of the current @em best @em indistinguishable protein group, until
   every peptide is uniquely assigned. This effectively allows more peptides to
   be used in ProteinQuantifier at the cost of potentially additional noise in
   the peptides quantities.
   In accordance with most state-of-the-art protein inference tools, only the
   best hit (PSM) for a peptide ID is considered.  Probability ties are
   currently resolved by taking the first occurring protein of the component.

   @TODO Implement probability tie resolution.
   @improvement The class could provide iterator for ConnectedComponents in the
   future. One could extend the graph to include all PeptideHits (not only the
   best). It becomes a tripartite graph with larger connected components then.
   Maybe extend it to work with MS1 features. Separate resolution and adding
   groups to output.
   
   @ingroup Analysis_ID
   */
  class OPENMS_DLLAPI PeptideProteinResolution
  {
  
  private:
    /// to build bipartite graph as two maps (adjacency "lists"):
    /// ProtGroups-Indices <-> PepID-Indices
    /// so we get bidirectional connectivity
    /// We always take first PepHit from PepID, because those are usually used
    /// for inference
    typedef std::map<Size, std::set<Size> > IndexMap_;

    /// if the protein group at index i contains a target (first) and/or decoy (second)
    // TODO WIP for better tie resolution
    // std::vector<std::pair<bool,bool>> indist_prot_grp_td_;
    /// mapping indist. protein group indices -> peptide identification indices
    IndexMap_ indist_prot_grp_to_pep_;
    /// mapping indist. protein group indices <- peptide identification indices
    IndexMap_ pep_to_indist_prot_grp_;
    
    /** represents the middle layer of an implicit tripartite graph:
    consists of single protein accessions and their mapping to the (indist.)
    group's indices */
    std::map<String, Size> prot_acc_to_indist_prot_grp_;
    
    /// log debug information?
    bool statistics_;
    
  public:
    /// Constructor
    /// @param statistics Specifies if the class stores/outputs info about statistics    
    PeptideProteinResolution(bool statistics = false);


    /// A peptide-centric reimplementation of the resolution process. Can be used statically
    /// without building a bipartite graph first.
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param peptides vector of ProteinIdentifications with links to the proteins
    /// @param resolve_ties If ties should be resolved or multiple best groups reported
    /// @param targets_first If target groups should get picked first no matter the posterior
    /// @todo warning: all peptides are used (not filtered for matching protein ID run yet).
    static void resolve(ProteinIdentification& protein,
                        std::vector<PeptideIdentification>& peptides,
                        bool resolve_ties,
                        bool targets_first);

    /// Initialize and store the graph (= maps), needs sorted groups for
    /// correct functionality. Therefore sorts the indist. protein groups
    /// if not skipped.
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param peptides vector of ProteinIdentifications with links to the proteins
    /// @param skip_sort Skips sorting of groups, nothing is modified then.
    void buildGraph(ProteinIdentification& protein,
                    const std::vector<PeptideIdentification>& peptides,
                    bool skip_sort = false);
      
    /// Applies resolveConnectedComponent to every component of the graph and
    /// is able to write statistics when specified. Parameters will
    /// both be mutated in this method.
    /// @param protein ProteinIdentification object storing IDs and groups
    /// @param peptides vector of ProteinIdentifications with links to the proteins
    /// @todo warning: all peptides are used (not filtered for matching protein ID run yet).
    void resolveGraph(ProteinIdentification& protein,
                      std::vector<PeptideIdentification>& peptides);
    
    /// Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
    /// and peptides), switching from one to the other in each step.
    /// @param root_prot_grp Starts the BFS at this protein group index
    /// @return Returns a Connected Component as set of group and peptide indices.
    ConnectedComponent findConnectedComponent(Size& root_prot_grp);
    

    /*! Resolves connected components based on posterior probabilities and adds them
     * as additional protein_groups to the output idXML.
     * Thereby greedily assigns shared peptides in this component uniquely to
     * the proteins of the current BEST INDISTINGUISHABLE protein group,
     * ready to be used in ProteinQuantifier then.
     * This is achieved by removing all other evidence from the input
     * PeptideIDs and iterating until each peptide is uniquely assigned.
     * In accordance with Fido only the best hit (PSM) for an ID is considered.
     * Probability ties are _currently_ resolved by taking the first occurrence.
     * @param conn_comp The component to be resolved
     * @param protein ProteinIdentification object storing IDs and groups
     * @param peptides vector of ProteinIdentifications with links to the proteins
     */
    void resolveConnectedComponent(ConnectedComponent& conn_comp,
                                    ProteinIdentification& protein,
                                    std::vector<PeptideIdentification>& peptides);
};
  
} //namespace OpenMS

