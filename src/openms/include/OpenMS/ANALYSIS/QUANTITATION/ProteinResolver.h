// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_PROTEINRESOLVER_H
#define OPENMS_ANALYSIS_QUANTITATION_PROTEINRESOLVER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{
  /**
    @brief Helper class for peptide and protein quantification based on feature data annotated with IDs

    This class is used by @ref TOPP_ProteinResolver. See there for further documentation.

    @htmlinclude OpenMS_ProteinResolver.parameters

    @ingroup Analysis_QUANTITATION

  */
  class OPENMS_DLLAPI ProteinResolver :
    public DefaultParamHandler
  {

public:

    //default constructor
    ProteinResolver();

    //copy constructor
    ProteinResolver(const ProteinResolver & rhs);

    //assignment operator
    ProteinResolver & operator=(const ProteinResolver & rhs);

    //destructor
    virtual ~ProteinResolver();


    struct ProteinEntry;
    struct PeptideEntry;
    struct ISDGroup;
    struct MSDGroup;
    struct ResolverResult;

    /// represents a protein from FASTA file
    struct ProteinEntry
    {
      std::list<PeptideEntry *> peptides;
      bool traversed;
      FASTAFile::FASTAEntry * fasta_entry;
      enum type  {primary, secondary, primary_indistinguishable, secondary_indistinguishable} protein_type;
      double weight;    //monoisotopic
      float coverage;    //in percent
      //if Protein is indistinguishable all his fellows are in the list indis
      std::list<ProteinEntry *> indis;
      Size index;
      Size  msd_group;     //index
      Size  isd_group;     //index
      Size number_of_experimental_peptides;
    };

    /// represents a peptide. First in silico. If experimental is set to true it is MS/MS derived.
    struct PeptideEntry
    {
      std::list<ProteinEntry *> proteins;
      bool traversed;
      String sequence;
      Size peptide_identification;
      Size peptide_hit;
      Size index;
      Size  msd_group;     //index
      Size isd_group;     //index
      bool experimental;
      float intensity;
      String origin;
    };

    /// representation of an msd group. Contains peptides, proteins and a pointer to its ISD group
    struct MSDGroup
    {
      std::list<ProteinEntry *> proteins;
      std::list<PeptideEntry *> peptides;
      Size index;
      ISDGroup * isd_group;
      Size number_of_decoy;
      Size number_of_target;
      Size number_of_target_plus_decoy;
      float intensity;     // intensity of the MSD Group. Defined as the median of the peptide intensities.
    };

    struct ISDGroup
    {
      std::list<ProteinEntry *> proteins;
      std::list<PeptideEntry *> peptides;
      Size index;
      std::list<Size> msd_groups;
    };

    struct ResolverResult
    {
      String identifier;
      std::vector<ISDGroup> * isds;
      std::vector<MSDGroup> * msds;
      std::vector<ProteinEntry> * protein_entries;
      std::vector<PeptideEntry> * peptide_entries;
      std::vector<Size> * reindexed_peptides;
      std::vector<Size> * reindexed_proteins;
      enum type  {PeptideIdent, Consensus} input_type;
      std::vector<PeptideIdentification> * peptide_identification;
      ConsensusMap * consensus_map;
    };

    /**
      @brief Computing protein groups from peptide identifications OR consensus map.

      Computes ISD and MSD groups.

      @param consensus ConsensusMap in case consensusXML file is given as input
    */
    void resolveConsensus(ConsensusMap & consensus);

    /**
      @brief Computing protein groups from peptide identifications OR consensus map.

      Computes ISD and MSD groups.

      @param peptide_identifications Vector of PeptideIdentification in case idXML is given as input
    */
    void resolveID(std::vector<PeptideIdentification> & peptide_identifications);

    /**
      @brief NOT IMPLEMENTED YET

      @param protein_nodes
      @param peptide_nodes
      @param reindexed_proteins
      @param reindexed_peptides
      @param peptide_identifications
      @param output
    */
    // void writeProteinsAndPeptidesmzTab(std::vector<ProteinEntry>& protein_nodes, std::vector<PeptideEntry>& peptide_nodes, std::vector<Size>& reindexed_proteins, std::vector<Size>& reindexed_peptides, std::vector<PeptideIdentification>& peptide_identifications, String& output  );
    /**
      @brief Writing peptide table into text file

      @param peptides
      @param reindexed_peptides
      @param identifications
      @param output_file
    */
    // void writePeptideTable(std::vector<PeptideEntry> & peptides, std::vector<Size> & reindexed_peptides, std::vector<PeptideIdentification> & identifications, String & output_file); // not implemented
    /**
      @brief Writing peptide table into text file

      @param peptides
      @param reindexed_peptides
      @param consensus
      @param output
    */
    // void writePeptideTable(std::vector<PeptideEntry> & peptides, std::vector<Size> & reindexed_peptides, ConsensusMap & consensus, String & output_file); // not implemented
    /**
      @brief Writing protein table into text file

      @param proteins
      @param reindexed_proteins
      @param output_file
    */
    // void writeProteinTable(std::vector<ProteinEntry> & proteins, std::vector<Size> & reindexed_proteins, String & output_file); // not implemented
    /**
      @brief Writing protein groups into text file

      @param isd_groups ISD groups
      @param msd_groups MSD groups
      @param output_file Path of output file
    */
    // void writeProteinGroups(std::vector<ISDGroup> & isd_groups, std::vector<MSDGroup> & msd_groups, String & output_file); // not implemented

    /**
      @brief brief

      @param msd_groups
      @param consensus
    */
    void countTargetDecoy(std::vector<MSDGroup> & msd_groups, ConsensusMap & consensus);

    /**
      @brief brief

      @param msd_groups
      @param peptide_nodes
    */
    void countTargetDecoy(std::vector<MSDGroup> & msd_groups, std::vector<PeptideIdentification> & peptide_nodes);

    void clearResult();

    void setProteinData(std::vector<FASTAFile::FASTAEntry> & protein_data);

    const  std::vector<ResolverResult> & getResults();

    /// overloaded functions -- return a const reference to a PeptideIdentification object or a peptideHit either from a consensusMap or a vector<PeptideIdentification>
    static const PeptideIdentification & getPeptideIdentification(const ConsensusMap & consensus, const PeptideEntry * peptide);
    static const PeptideHit & getPeptideHit(const ConsensusMap & consensus, const PeptideEntry * peptide);
    static const PeptideIdentification & getPeptideIdentification(const std::vector<PeptideIdentification> & peptide_nodes, const PeptideEntry * peptide);
    static const PeptideHit & getPeptideHit(const std::vector<PeptideIdentification> & peptide_nodes, const PeptideEntry * peptide);

private:

    std::vector<ResolverResult> resolver_result_;
    std::vector<FASTAFile::FASTAEntry> protein_data_;

    void computeIntensityOfMSD_(std::vector<MSDGroup> & msd_groups);

    /// traverse protein and peptide nodes for building MSD groups
    void traverseProtein_(ProteinEntry * prot_node, MSDGroup & group);
    void traversePeptide_(PeptideEntry * pep_node, MSDGroup & group);
    /// searches given sequence in all nodes and returns its index or nodes.size() if not found.
    Size findPeptideEntry_(String seq, std::vector<PeptideEntry> & nodes);
    /// helper function for findPeptideEntry.
    Size binarySearchNodes_(String & seq, std::vector<PeptideEntry> & nodes, Size start, Size end);
    /// includes all MS/MS derived peptides into the graph --idXML
    Size includeMSMSPeptides_(std::vector<PeptideIdentification> & peptide_identifications, std::vector<PeptideEntry> & peptide_nodes);
    /// TODO include run information for each peptide
    /// includes all MS/MS derived peptides into the graph --consensusXML
    Size includeMSMSPeptides_(ConsensusMap & consensus, std::vector<PeptideEntry> & peptide_nodes);
    /// Proteins and Peptides get reindexed, based on whether they belong to msd groups or not. Indexes of Proteins which are in an ISD group but in none of the MSD groups will not be used anymore.
    void reindexingNodes_(std::vector<MSDGroup> & msd_groups, std::vector<Size> & reindexed_proteins, std::vector<Size> & reindexed_peptides);
    /// marks Proteins which have a unique peptide as primary. Uses reindexed vector, thus reindexingNodes has to be called before.
    void primaryProteins_(std::vector<PeptideEntry> & peptide_nodes, std::vector<Size> & reindexed_peptides);
    void buildingMSDGroups_(std::vector<MSDGroup> & msd_groups, std::vector<ISDGroup> & isd_groups);
    void buildingISDGroups_(std::vector<ProteinEntry> & protein_nodes, std::vector<PeptideEntry> & peptide_nodes,
                            std::vector<ISDGroup> & isd_groups);
    // disabled/buggy
    //ProteinResolver::indistinguishableProteins(vector<MSDGroup>& msd_groups);

  }; // class

} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_PROTEINRESOLVER_H
