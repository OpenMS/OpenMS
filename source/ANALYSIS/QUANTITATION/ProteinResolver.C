// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

//#include <algorithm>

using std::map;
using std::vector;
using std::list;
using std::set;

//using namespace std;

namespace OpenMS
{

  ProteinResolver::ProteinResolver()
    : DefaultParamHandler("ProteinResolver"), resolver_result_(), protein_data_()
  {
    defaults_.setValue("resolver:missed_cleavages", 2, "Number of allowed missed cleavages");
    defaults_.setMinInt("resolver:missed_cleavages", 0);
    defaults_.setValue("resolver:min_length", 6, "Minimum length of peptide");
    defaults_.setMinInt("resolver:min_length", 1);
    defaults_.setValue("resolver:enzyme", "Trypsin", "Digestion enzyme");
    defaults_.setValidStrings("resolver:enzyme", StringList::create("Trypsin"));

    defaults_.setSectionDescription("resolver", "Additional options for algorithm");

    defaultsToParam_();
  }

  ProteinResolver::ProteinResolver(const ProteinResolver& cp)
    : DefaultParamHandler(cp),
    resolver_result_(cp.resolver_result_),
    protein_data_(cp.protein_data_)
  {
  }


  ProteinResolver& ProteinResolver::operator= (const ProteinResolver& rhs)
  {
    if (this == &rhs) return *this;

    DefaultParamHandler::operator= (rhs);
    this->resolver_result_ = rhs.resolver_result_;
    this->protein_data_ = rhs.protein_data_;

    return *this;
  }


  ProteinResolver::~ProteinResolver()
  {
    clearResult();
  }

  void ProteinResolver::clearResult()
  {
    for (vector<ResolverResult>::iterator it = resolver_result_.begin(); it != resolver_result_.end(); ++it )
    {
      delete it->isds;
      delete it->msds;
      delete it->peptide_entries;
      delete it->protein_entries;
      delete it->reindexed_peptides;
      delete it->reindexed_proteins;
    }
    resolver_result_.clear();
  }

  void ProteinResolver::resolveID(vector<PeptideIdentification>& peptide_identifications)
  {
    vector<ProteinEntry>* protein_nodes = new vector<ProteinEntry>;
    protein_nodes->resize(protein_data_.size());
    vector<PeptideEntry>* peptide_nodes = new vector<PeptideEntry>;
    vector<ISDGroup>* isd_groups = new vector<ISDGroup>;
    vector<MSDGroup>* msd_groups = new vector<MSDGroup>;
    vector<Size>* reindexed_proteins = new vector<Size>;
    vector<Size>* reindexed_peptides = new vector<Size>;

    // building ISD Groups
    buildingISDGroups_(*protein_nodes, *peptide_nodes, *isd_groups);

    // Including all MSMS derived peptides into the graph
    Size found_peptides;
    found_peptides = includeMSMSPeptides_(peptide_identifications,*peptide_nodes);

    // building MSD Groups
    buildingMSDGroups_(*msd_groups, *isd_groups);

    // calculations + reindexing
    reindexingNodes_(*msd_groups, *reindexed_proteins, *reindexed_peptides);
    primaryProteins_(*peptide_nodes, *reindexed_peptides);
    //TODO indistinguishableProteins(msd_groups);

    //
    countTargetDecoy(*msd_groups, peptide_identifications);

    ResolverResult result;
    //result.identifier = file_identifier;
    result.isds = isd_groups;
    result.msds = msd_groups;
    result.peptide_entries = peptide_nodes;
    result.protein_entries = protein_nodes;
    result.reindexed_peptides = reindexed_peptides;
    result.reindexed_proteins = reindexed_proteins;
    result.input_type = ResolverResult::PeptideIdent;
    result.peptide_identification = &peptide_identifications;

    resolver_result_.push_back(result);
  }

  void ProteinResolver::resolveConsensus(ConsensusMap& consensus)
  {
    vector<ProteinEntry>* protein_nodes = new vector<ProteinEntry>;
    protein_nodes->resize(protein_data_.size());
    vector<PeptideEntry>* peptide_nodes = new vector<PeptideEntry>;
    vector<ISDGroup>* isd_groups = new vector<ISDGroup>;
    vector<MSDGroup>* msd_groups = new vector<MSDGroup>;
    vector<Size>* reindexed_proteins = new vector<Size>;
    vector<Size>* reindexed_peptides = new vector<Size>;

    // building ISD Groups
    buildingISDGroups_(*protein_nodes, *peptide_nodes, *isd_groups);

    // Including all MSMS derived peptides into the graph
    Size found_peptides;
    found_peptides = includeMSMSPeptides_(consensus,*peptide_nodes);

    // building MSD Groups
    buildingMSDGroups_(*msd_groups, *isd_groups);

    // calculations + reindexing
    reindexingNodes_(*msd_groups, *reindexed_proteins, *reindexed_peptides);

    // compute intensity of a msd group
    computeIntensityOfMSD_(*msd_groups);

    primaryProteins_(*peptide_nodes, *reindexed_peptides);
    //TODO indistinguishableProteins(msd_groups);

    //
    countTargetDecoy(*msd_groups, consensus);

    ResolverResult result;
    //result.identifier = file_identifier;
    result.isds = isd_groups;
    result.msds = msd_groups;
    result.peptide_entries = peptide_nodes;
    result.protein_entries = protein_nodes;
    result.reindexed_peptides = reindexed_peptides;
    result.reindexed_proteins = reindexed_proteins;
    result.input_type = ResolverResult::Consensus;
    result.consensus_map = &consensus;

    resolver_result_.push_back(result);
  }

  void ProteinResolver::computeIntensityOfMSD_(vector<MSDGroup>& msd_groups)
  {
    // iteriert ueber alles msd gruppe
    for(vector<MSDGroup>::iterator group = msd_groups.begin(); group != msd_groups.end(); ++group)
    {
      DoubleList intensities;
      // iterierere ueber peptide entry (peptide identification), intensitaet (summe der einzelintensitaeten)
      for(list<PeptideEntry*>::iterator pep = group->peptides.begin(); pep != group->peptides.end(); ++pep)
      {
        intensities.push_back((*pep)->intensity);
      }
      // median von der list ist itensity der msd group
      group->intensity = Math::median(intensities.begin(),intensities.end());
    }
  }

  void ProteinResolver::countTargetDecoy(vector<MSDGroup>& msd_groups, ConsensusMap& consensus)
  {
    for(vector<MSDGroup>::iterator group = msd_groups.begin(); group != msd_groups.end(); ++group)
    {
      for(list<PeptideEntry*>::iterator pep = group->peptides.begin(); pep != group->peptides.end(); ++pep)
      {
        String tmp = getPeptideHit(consensus, *pep).getMetaValue("target_decoy");

        if(tmp == "target") ++(group->number_of_target);
        else if(tmp == "decoy") ++(group ->number_of_decoy);
        else /*if(tmp == "target+decoy")*/ ++(group->number_of_target_plus_decoy);
      }
    }
  }

  void ProteinResolver::countTargetDecoy(vector<MSDGroup>& msd_groups, vector<PeptideIdentification>& peptide_nodes)
  {
    for(vector<MSDGroup>::iterator group = msd_groups.begin(); group != msd_groups.end(); ++group)
    {
      for(list<PeptideEntry*>::iterator pep = group->peptides.begin(); pep != group->peptides.end(); ++pep)
      {
        String tmp = getPeptideHit(peptide_nodes, *pep).getMetaValue("target_decoy");

        if(tmp == "target") ++(group->number_of_target);
        else if(tmp == "decoy") ++(group ->number_of_decoy);
        else /*if(tmp == "target+decoy")*/ ++(group->number_of_target_plus_decoy);
      }
    }
  }

  //travers Protein and peptide nodes. Once for building ISD groups and once for building MSD groups
  void ProteinResolver::traversProtein_(ProteinEntry* prot_node, ISDGroup& group)
  {
    group.proteins.push_back(prot_node);
    prot_node->isd_group = group.index;
    for(list<PeptideEntry*>::iterator i = prot_node->peptides.begin(); i != prot_node->peptides.end();++i)
    {
      if(!(*i)->traversed)
      {
        (*i)->traversed = true;
        traversPeptide_((*i), group);
      }
    }
  }

  void ProteinResolver::traversProtein_(ProteinEntry* prot_node, MSDGroup& group)
  {
    group.proteins.push_back(prot_node);
    prot_node->msd_group = group.index;
    for(list<PeptideEntry*>::iterator i = prot_node->peptides.begin(); i != prot_node->peptides.end();++i)
    {
      if((*i)->experimental) ++(prot_node->number_of_experimental_peptides);
      if((*i)->traversed)
      {
        (*i)->traversed = false;
        if((*i)->experimental)
        {
          traversPeptide_((*i), group);
        }
      }
    }
  }

  void ProteinResolver::traversPeptide_(PeptideEntry* pep_node, ISDGroup& group)
  {
    group.peptides.push_back(pep_node);
    pep_node->isd_group = group.index;
    for(list<ProteinEntry*>::iterator i = pep_node->proteins.begin(); i != pep_node->proteins.end();++i)
    {
      if(!(*i)->traversed)
      {
        (*i)->traversed = true;
        traversProtein_((*i), group);
      }
    }
  }

  void ProteinResolver::traversPeptide_(PeptideEntry* pep_node, MSDGroup& group)
  {
    group.peptides.push_back(pep_node);
    pep_node->msd_group = group.index;
    for(list<ProteinEntry*>::iterator i = pep_node->proteins.begin(); i != pep_node->proteins.end();++i)
    {
      if((*i)->traversed)
      {
        (*i)->traversed = false;
        traversProtein_((*i), group);
      }
    }
  }

  //searches given sequence in all  nodes and returns its index or nodes.size() if not found.
  Size ProteinResolver::findPeptideEntry_(String seq, vector<PeptideEntry>& nodes)
  {
    if(nodes.size() == 0) return 0;
    return binarySearchNodes_(seq, nodes, 0, nodes.size()-1);
  }

  //helper function for findPeptideEntry.
  Size ProteinResolver::binarySearchNodes_(String& seq, vector<PeptideEntry>& nodes, Size start, Size end )
  {

    // Termination condition: start index greater than end index
    if(start > end) return -1;

    Size compare_value = (start + end)/2;
    String& node_sequence = nodes[compare_value].sequence;
    int compar = seq.compare(node_sequence);
    if(start == end)
    {
      if(compar != 0) return nodes.size();
      else 					 return compare_value;
    }
    else if(compar <  0)
    {
      --compare_value;
      if (compare_value < start) compare_value = start;

      return binarySearchNodes_(seq, nodes,start,compare_value);
    }
    else if(compar > 0)
    {
      ++compare_value;
      if(compare_value > end) compare_value = end;

      return binarySearchNodes_(seq,nodes,compare_value, end);
    }
    else
    {
      return compare_value;
    }
  }

  //includes all MSMS derived peptides into the graph --IdXML
  Size ProteinResolver::includeMSMSPeptides_(vector<PeptideIdentification> &peptide_identifications, vector<PeptideEntry> &peptide_nodes )
  {
    Size found_peptide = 0;
    for(Size pep = 0; pep != peptide_identifications.size();++pep)
    {
      Size peptide_entry = findPeptideEntry_(peptide_identifications[pep].getHits().front().getSequence().toUnmodifiedString(),peptide_nodes);

      if(peptide_entry != peptide_nodes.size())
      {
        if(!peptide_nodes[peptide_entry].experimental)
        {
          ++found_peptide;
        }
        peptide_nodes[peptide_entry].peptide_identification = pep;
        peptide_nodes[peptide_entry].peptide_hit = 0; //only top hit is used at the moment
        peptide_nodes[peptide_entry].experimental = true;
      }
    }
    return found_peptide;
  }

  //TODO include run information for each peptide
  //includes all MSMS derived peptides into the graph --ConsensusXML
  Size ProteinResolver::includeMSMSPeptides_(ConsensusMap & consensus, vector<PeptideEntry> &peptide_nodes )
  {
    Size found_peptide = 0;
    for(Size pep = 0; pep != consensus.size(); ++pep)
    {
      ConsensusFeature& feature = consensus.at(pep);

      // get all peptide identifications
      const vector<PeptideIdentification> & pep_id  = feature.getPeptideIdentifications();


      for(Size cons_pep = 0; cons_pep < pep_id.size(); ++cons_pep)
      {
        String seq = pep_id.at(cons_pep).getHits().front().getSequence().toUnmodifiedString();
        Size peptide_entry = findPeptideEntry_(seq, peptide_nodes);

        if(peptide_entry != peptide_nodes.size())
        {
          if(!peptide_nodes.at(peptide_entry).experimental)
          {
            ++found_peptide;
          }
          //should be changed -- for consensus peptide_identification is the consensus and peptide_hit is the PeptideIdentification. PeptideHit is only top hit at the moment
          peptide_nodes.at(peptide_entry).peptide_identification = pep;
          peptide_nodes.at(peptide_entry).peptide_hit = cons_pep; //only top hit is used at the moment
          peptide_nodes.at(peptide_entry).experimental = true;
          // get intensity of the feature
          peptide_nodes.at(peptide_entry).intensity = feature.getIntensity();
          peptide_nodes.at(peptide_entry).origin = feature.getMetaValue("file_origin");
        }
      }
    }
    return found_peptide;
  }

  //overloaded functions -- return a const reference to a PeptideIdentification object or a peptideHit either from a consensusMap or a vector<PeptideIdentification>
  const PeptideIdentification& ProteinResolver::getPeptideIdentification(const ConsensusMap& consensus, const PeptideEntry* peptide)
  {
    return (consensus[peptide->peptide_identification].getPeptideIdentifications()[peptide->peptide_hit]);
  }
  const PeptideHit& ProteinResolver::getPeptideHit(const ConsensusMap& consensus,const PeptideEntry* peptide)
  {
    return getPeptideIdentification(consensus,peptide).getHits().front();
  }

  const PeptideIdentification& ProteinResolver::getPeptideIdentification(const vector<PeptideIdentification>& peptide_nodes, const PeptideEntry* peptide)
  {
    return peptide_nodes[peptide->peptide_identification];
  }
  const PeptideHit& ProteinResolver::getPeptideHit(const vector<PeptideIdentification>& peptide_nodes,const PeptideEntry* peptide)
  {
    return getPeptideIdentification(peptide_nodes,peptide).getHits().front();
  }


  //Proteins and Peptides get reindexed, based on whether they belong to msd groups or not. Indexes of Proteins which are in an ISD group but in none of the MSD groups will not be used anymore.
  void ProteinResolver::reindexingNodes_(vector<MSDGroup>& msd_groups, vector<Size> &reindexed_proteins,vector<Size> &reindexed_peptides )
  {
    Size new_prot_index(0);
    Size new_pep_index(0);
    for(vector<MSDGroup>::iterator msd = msd_groups.begin();msd != msd_groups.end(); ++msd)
    {
      for(list<ProteinEntry*>::iterator prot = msd->proteins.begin(); prot != msd->proteins.end(); ++prot)
      {
        reindexed_proteins.push_back((*prot)->index);
        (*prot)->index = new_prot_index;
        ++new_prot_index;
      }

      for(list<PeptideEntry*>::iterator pep = msd->peptides.begin(); pep != msd->peptides.end(); ++pep)
      {
        reindexed_peptides.push_back((*pep)->index);
        (*pep)->index = new_pep_index;
        ++new_pep_index;
      }
    }
  }

  //marks Proteins which have a unique peptide as primary. Uses reindexed vector, thus reindexingNodes has to be called before.
  void ProteinResolver::primaryProteins_(vector<PeptideEntry> &peptide_nodes,  vector<Size>& reindexed_peptides)
  {
    //primary proteins
    Size pc = 0;
    for(vector<Size>::iterator pep = reindexed_peptides.begin(); pep != reindexed_peptides.end(); ++pep)
    {
      ++pc;
      if(peptide_nodes[*pep].proteins.size() == 1)
      {
        peptide_nodes[*pep].proteins.front()->protein_type = ProteinEntry::primary;
      }
    }
  }

  void ProteinResolver::buildingISDGroups_(vector<ProteinEntry>& protein_nodes, vector<PeptideEntry>& peptide_nodes,
                                           vector<ISDGroup>& isd_groups)
  {
    EnzymaticDigestion digestor;
    String enzyme_name = param_.getValue("resolver:enzyme");
    digestor.setEnzyme(digestor.getEnzymeByName(enzyme_name));
    UInt min_size = param_.getValue("resolver:min_length");
    UInt missed_cleavages = param_.getValue("resolver:missed_cleavages");
    digestor.setMissedCleavages(missed_cleavages);


    //-------------------------------------------------------------
    // building ISD Groups
    //-------------------------------------------------------------

    vector<AASequence> temp_peptides;
    map<String,set<Size> > peptides;

    for (Size i = 0; i < protein_data_.size(); ++i)
    {
      protein_nodes[i].fasta_entry = &protein_data_[i];
      protein_nodes[i].traversed = false;
      protein_nodes[i].index = i;
      protein_nodes[i].protein_type = ProteinEntry::secondary;
      protein_nodes[i].weight = AASequence(protein_data_[i].sequence).getMonoWeight();
      protein_nodes[i].coverage = 0.;
      protein_nodes[i].number_of_experimental_peptides = 0;
      digestor.digest(AASequence(protein_data_[i].sequence), temp_peptides);
      for (Size j = 0; j < temp_peptides.size(); ++j)
      {
        if (temp_peptides[j].size() >= min_size)
        {
          peptides[temp_peptides[j].toUnmodifiedString()].insert(i);
        }
      }
    }
    // important to resize
    peptide_nodes.resize(peptides.size());
    vector<PeptideEntry>::iterator pep_node = peptide_nodes.begin();
    Size peptide_counter = 0;

    for(map<String,set<Size> >::iterator i  = peptides.begin(); i != peptides.end(); ++i, ++pep_node, ++peptide_counter)
    {
      pep_node->index = peptide_counter;
      pep_node->traversed = false;
      pep_node->sequence = (*i).first;
      pep_node->experimental = false;
      for (set<Size>::iterator j = (*i).second.begin(); j != (*i).second.end();++j)
      {
        pep_node->proteins.push_back(&protein_nodes[*j]);
        protein_nodes[*j].peptides.push_back(&*pep_node);
      }
    }
    //ISDGraph constructed
    Size isd_group_counter = 0;
    Size i = 0;
    for(vector<ProteinEntry>::iterator prot_node = protein_nodes.begin(); prot_node != protein_nodes.end();++prot_node)
    {
      ++i;
      if(! prot_node -> traversed)
      {
        prot_node -> traversed = true;
        ISDGroup group;
        group.index = isd_group_counter;
        ++isd_group_counter;
        traversProtein_(&*prot_node,group);
        isd_groups.push_back(group);
      }
    }
  }

  void ProteinResolver::buildingMSDGroups_(vector<MSDGroup>& msd_groups, vector<ISDGroup>& isd_groups)
  {
    //-------------------------------------------------------------
    // building MSDGroups
    //-------------------------------------------------------------
    Size msd_group_counter = 0;
    for(Size isd_group = 0; isd_group != isd_groups.size(); ++isd_group)
    {
      for(list<ProteinEntry*>::iterator prot = isd_groups[isd_group].proteins.begin(); prot != isd_groups[isd_group].proteins.end(); ++prot)
      {
        ProteinEntry* prot_node = (*prot);
        if(prot_node -> traversed)
        {
          prot_node -> traversed = false;
          MSDGroup msd_group;
          msd_group.index = msd_group_counter;
          msd_group.isd_group = &(isd_groups[isd_group]);
          msd_group.number_of_target = 0;
          msd_group.number_of_decoy= 0;
          msd_group.number_of_target_plus_decoy = 0;
          traversProtein_(&*prot_node,msd_group);
          if(msd_group.peptides.size()> 0)
          {
            msd_groups.push_back(msd_group);
            isd_groups[isd_group].msd_groups.push_back(msd_group_counter);
            ++msd_group_counter;
          }
        }
      }
    }
  }

  void ProteinResolver::setProteinData(vector<FASTAFile::FASTAEntry>& protein_data)
  {
    protein_data_ = protein_data;
  }

  const vector< ProteinResolver::ResolverResult >& ProteinResolver::getResults()
  {
    return resolver_result_;
  }

  //not tested
  //ProteinResolver::indistinguishableProteins(vector<MSDGroup>& msd_groups)
  //{
  //indistinguishable proteins
  //TODO check which proteins are indistinguishable from eachother. This code seems to be instable.
  /*
  Size mg = 0;
  cout<<"\nindistinguishable\n";
  for(vector<MSDGroup>::iterator group = msd_groups.begin(); group != msd_groups.end(); ++group)
  {
    cout<<mg<<" msd groups of " <<msd_groups.size() << " done\r";
    ++mg;
    multimap<Size,ProteinEntry*> same_number_of_peptides;
    for(list<ProteinEntry*>::iterator prot = group->proteins.begin(); prot != group->proteins.end(); ++prot)
    {
      same_number_of_peptides.insert(pair<Size,ProteinEntry*> ((*prot)->peptides.size(),(*prot)));
    }
    Size i = same_number_of_peptides.begin()->first;
    do{
      if(i == 0) continue;
      multimap<Size,ProteinEntry*>::iterator it1,it2, itlow,itup;
      itlow = same_number_of_peptides.lower_bound(i);
      itup = same_number_of_peptides.upper_bound(i);

      //indistinguishable(itlow,itup);
      for(it1 = itlow; it1 != itup; ++it1)
      {
        if(it1->second->protein_type != ProteinEntry::secondary) continue;
        bool same_peptides = true;
        for(it2 = (++itlow); it2 != itup;++it2)
        {
          bool indi = false;
          for(list<PeptideEntry*>::iterator pep1 = it1->second->peptides.begin();pep1 != it1->second->peptides.end(); ++pep1)
          {
            for(list<PeptideEntry*>::iterator pep2 = it2->second->peptides.begin();pep2 != it2->second->peptides.end(); ++pep2)
            {
              //hope that pointer arithmetic works as expected
            //	if(test_vec[(*pep1)->peptide_identification].getHits()[(*pep1)->peptide_hit].getSequence()== test_vec[(*pep2)->peptide_identification].getHits()[(*pep2)->peptide_hit].getSequence())
              if((*pep1) == (*pep2))
              {
                indi = true;
                break;
              }
            }
            if(!indi)
            {
              same_peptides = false;
              break;
            }
          }
          if(same_peptides)
          {
            it1->second->protein_type = ProteinEntry::secondary_indistinguishable;
            it2->second->protein_type = ProteinEntry::secondary_indistinguishable;
            it1->second->indis.push_back(it2->second);
            it2->second->indis.push_back(it1->second);
          }
        }
        //TODO gurantee that all indistinguishable have all others in their list
        //at the moment the key idea is that only the first Protein has all of them. The others have just a pointer to the first one.
      }
    ++i;
    }while( i <= same_number_of_peptides.rbegin()->first);
  }*/
  //TODO exclusively at least one peptides that they share -> primary indistinguishable!
  //Find out which are primary but indistinguishable
  /*cout<<"\nprimary indis\n";
  mg = 0;
  for(vector<MSDGroup>::iterator group  = msd_groups.begin(); group != msd_groups.end(); ++group)
  {
    cout<<mg<< " msd group of " << msd_groups.size() <<" done\r";
    ++mg;
    for(list<ProteinEntry*>::iterator prot_node = group->proteins.begin(); prot_node != group->proteins.end(); ++prot_node)
    {
      ProteinEntry* prot = (*prot_node);
      if(!prot->traversed)
      {
        prot->traversed = true;
        if(prot->indis.size() > 0)
        {
          Size number_of_proteins = prot->indis.size()+1;
          bool primary_indistinguishable;
          for(list<PeptideEntry*>::iterator pep =  prot->peptides.begin(); pep != prot->peptides.end(); ++pep)
          {
            PeptideEntry* p = (*pep);
            if(p->proteins.size() == number_of_proteins)
            {
              primary_indistinguishable = true;
              break;
            }
          }
          if(primary_indistinguishable)prot->protein_type = ProteinEntry::primary_indistinguishable;
          for(list<ProteinEntry*>::iterator indi = prot->indis.begin(); indi != prot->indis.end(); ++indi)
          {
            (*indi)->index = prot->index;
            if(primary_indistinguishable) (*indi)->protein_type = ProteinEntry::primary_indistinguishable;
            (*indi)->traversed = true;
          }
        }
      }
    }
  }
 }*/

}


