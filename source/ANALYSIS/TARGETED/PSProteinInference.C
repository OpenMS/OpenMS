// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//
#include <OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/CONCEPT/LogStream.h>


#include <algorithm>

using namespace std;
//#define PIS_DEBUG
#undef PIS_DEBUG
namespace OpenMS
{
  PSProteinInference::PSProteinInference():solver_(LPWrapper::SOLVER_GLPK)
  {
  }
  
  PSProteinInference::~PSProteinInference()
  {
  }

  
  Size PSProteinInference::findMinimalProteinList(const std::vector<PeptideIdentification>& peptide_ids)
  {
    // // first map peptides to proteins
    // std::map<String,vector<PeptideIdentification> > pep_prot_map;
    // for(Size i = 0; i < ids.size();++i)
    //   {
    //     // consider only first peptide hit -> should be filtered before
    //     if(ids[i].getHits().empty() || ids[i].getHits().size() > 1)
    //       {
    //         LOG_FATAL_ERROR << "peptide id contains more than 1 peptide hit -> filter for best hits before using PSProteinInference!";
    //         throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peptide Id contains more than 1 peptide hit", String(ids[i].getHits().size()));
    //       }
        
    //     const vector<String> & accs = ids[i].getHits()[0].getProteinAccessions();
    //     String seq = ids[i].getHits()[0].getSequence().toUnmodifiedString();
    //     DoubleReal score = ids[i].getHits()[0].getScore();
    //     bool higher_better = ids[i].isHigherScoreBetter();
    //     for(Size a = 0; a < accs.size(); ++a)
    //       {
    //         // first check if protein has already a peptide with the same sequence stored
    //         std::vector<PeptideIdentification>::iterator it = pep_prot_map[accs[a]].begin();
    //         bool insert = true;
    //         for(;it != pep_prot_map[accs[a]].end();++it)
    //           {
    //             if(it->getHits()[0].getSequence().toUnmodifiedString() == seq) // find peptide sequence
    //               {
    //                 if((higher_better && it->getHits()[0].getScore() > score)
    //                    || (!higher_better && it->getHits()[0].getScore() < score))// if there is the same peptide with better score
    //                   {
    //                     insert = false; // new pep id is not inserted
    //                     break;
    //                   }
    //                 else if((higher_better && it->getHits()[0].getScore() < score)
    //                         || (!higher_better && it->getHits()[0].getScore() > score)) // if there is the same peptide with a worse score
    //                   {
    //                     *it = ids[i];
    //                     insert = false; // it is replaced by the new one
    //                     break;
    //                   }
    //               }
    //           }
    //         if(insert) pep_prot_map[accs[a]].push_back(ids[i]);
    //       }
    //   }


    LPWrapper problem;
    problem.setSolver(solver_);
    set<String> all_accs;
    problem.setObjectiveSense(LPWrapper::MIN);
    minimal_protein_list_accessions_.clear();
      
    // first get all protein accessions:
    for(Size p = 0; p < peptide_ids.size();++p)
      {
        const vector<String> & accs = peptide_ids[p].getHits()[0].getProteinAccessions();
//         std::cout << peptide_ids[p].getHits()[0].getSequence()<<"Peptide Id with "
//                   << accs.size() << " protein accs.\n";
        for(Size a = 0; a < accs.size(); ++a)
          {
            all_accs.insert(accs[a]);
          }
      }
    
    // add variable for each protein:
    for(set<String>::const_iterator p = all_accs.begin(); p != all_accs.end();++p)
      {
        // create column with boundaries 0-1 and integer/binary variable
        Size index = problem.addColumn();
        problem.setColumnBounds(index,0.,1.,LPWrapper::DOUBLE_BOUNDED);
        problem.setColumnName(index,*p);
        problem.setColumnType(index,LPWrapper::BINARY);
        problem.setObjective(index,1.);
      }

    // now go through all peptide ids and add an constraint for each
    for(Size p = 0; p < peptide_ids.size();++p)
      {
        // consider only first peptide hit -> should be filtered before
        if(peptide_ids[p].getHits().size() > 1)
          {
            LOG_FATAL_ERROR << "peptide id contains more than 1 peptide hit -> filter for best hits before using PSProteinInference!";
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peptide Id contains more than 1 peptide hit", String(peptide_ids[p].getHits().size()));
          }
        
        // get column indices for all corresponding proteins
        vector<Int> indices;
        const vector<String> & accs = peptide_ids[p].getHits()[0].getProteinAccessions();
        for(Size a = 0; a < accs.size(); ++a)
          {
            indices.push_back((Int)problem.getColumnIndex(accs[a]));
          }
        vector<DoubleReal> values(indices.size(),1.);

        // enter constraint
        problem.addRow(indices,values,String(p) + peptide_ids[p].getHits()[0].getSequence().toString(),1.,1.,LPWrapper::LOWER_BOUND_ONLY);
      }

    // solve problem
    LPWrapper::SolverParam param;
    problem.solve(param);
    //    problem.writeProblem("prote_inference.mps","MPS");
    // get solution and store accession 
    for(Int c = 0; c < problem.getNumberOfColumns();++c)
      {
        if(problem.getColumnValue(c) == 1)
          {
            minimal_protein_list_accessions_.push_back(problem.getColumnName(c)); // enter protein accession
          }
      }
    
    return  minimal_protein_list_accessions_.size();
  }
  
  void PSProteinInference::calculateProteinProbabilities(const std::vector<PeptideIdentification>& ids)
  {
    accessions_.clear(); probabilities_.clear();
    
    // first map peptides to proteins
    std::map<String,vector<PeptideIdentification> > pep_prot_map;
    for(Size i = 0; i < ids.size();++i)
      {
        // consider only first peptide hit -> should be filtered before
        if(ids[i].getHits().empty() || ids[i].getHits().size() > 1)
          {
            LOG_FATAL_ERROR << "peptide id contains more than 1 peptide hit -> filter for best hits before using PSProteinInference!";
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peptide Id contains more than 1 peptide hit", String(ids[i].getHits().size()));
          }
        
        const vector<String> & accs = ids[i].getHits()[0].getProteinAccessions();
        String seq = ids[i].getHits()[0].getSequence().toUnmodifiedString();
        DoubleReal score = ids[i].getHits()[0].getScore();
        bool higher_better = ids[i].isHigherScoreBetter();
        for(Size a = 0; a < accs.size(); ++a)
          {
            // first check if protein has already a peptide with the same sequence stored
            std::vector<PeptideIdentification>::iterator it = pep_prot_map[accs[a]].begin();
            bool insert = true;
            for(;it != pep_prot_map[accs[a]].end();++it)
              {
                if(it->getHits()[0].getSequence().toUnmodifiedString() == seq) // find peptide sequence
                  {
                    if((higher_better && it->getHits()[0].getScore() >= score)
                       || (!higher_better && it->getHits()[0].getScore() <= score))// if there is the same peptide with better score
                      {
                        insert = false; // new pep id is not inserted
                        break;
                      }
                    else if((higher_better && it->getHits()[0].getScore() < score)
                            || (!higher_better && it->getHits()[0].getScore() > score)) // if there is the same peptide with a worse score
                      {
                        *it = ids[i];
                        insert = false; // it is replaced by the new one
                        break;
                      }
                    insert = false;
                    break;
                  }
              }
            if(insert) pep_prot_map[accs[a]].push_back(ids[i]);
          }
      }

    // enter accessions in member vector
    std::map<String,vector<PeptideIdentification> >::const_iterator map_it = pep_prot_map.begin();
    for(;map_it != pep_prot_map.end(); ++map_it) accessions_.push_back(map_it->first);
    
    // now calculate protein probabilities
    probabilities_.assign(accessions_.size(),1.);
    for(Size a = 0; a < accessions_.size();++a)
      {
        if(pep_prot_map.find(accessions_[a]) == pep_prot_map.end())
          {
            probabilities_[a] = 1.;
            continue;
          }
        //        std::cout << accessions_[a] << " "<< pep_prot_map[accessions_[a]].size()<<std::endl;
        for(Size p = 0; p < pep_prot_map[accessions_[a]].size();++p)
          {
            // std::cout <<pep_prot_map[accessions_[a]][p].getHits()[0].getSequence() << " "
            //           <<pep_prot_map[accessions_[a]][p].getHits()[0].getScore() <<std::endl;
            bool higher_better = pep_prot_map[accessions_[a]][p].isHigherScoreBetter();
            //std::cout << higher_better << " "<<(pep_prot_map[accessions_[a]][p].getHits()[0].getScore()) << " "<<probabilities_[a]<<"  -> ";
            if(higher_better) probabilities_[a] *= (1.-pep_prot_map[accessions_[a]][p].getHits()[0].getScore());
            else probabilities_[a] *= (pep_prot_map[accessions_[a]][p].getHits()[0].getScore());
            //std::cout << probabilities_[a]<<std::endl;
          }
        //        std::cout << "\n";
      }

    for(Size a = 0; a < probabilities_.size();++a)  probabilities_[a] = 1. - probabilities_[a];

  }


//   DoubleReal PSProteinInference::getProteinProbability(const String& acc,const std::vector<String>& accessions, const std::vector<DoubleReal>& probabilities)
//   {
//     std::vector<String>::const_iterator it = std::find(accessions.begin(),accessions.end(),acc);
//     if(it == accessions.end())  return 0.;
//     return probabilities[distance(accessions.begin(),it)];
//   }
  

  DoubleReal PSProteinInference::getProteinProbability(const String& acc)
  {
    std::vector<String>::iterator it = std::find(accessions_.begin(),accessions_.end(),acc);
    if(it == accessions_.end())  return 0.;
    return probabilities_[distance(accessions_.begin(),it)];
  }

  bool PSProteinInference::isProteinInMinimalList(const String& acc)
  {
    return find(minimal_protein_list_accessions_.begin(),minimal_protein_list_accessions_.end(),acc) != minimal_protein_list_accessions_.end();
  }

  Int PSProteinInference::getNumberOfProtIds(DoubleReal protein_id_threshold)
  {
    Int num = 0;
    for(Size i = 0; i < minimal_protein_list_accessions_.size();++i)
      {
        // std::cout << minimal_protein_list_accessions_[i] << " "<<getProteinProbability(minimal_protein_list_accessions_[i])
        //           << " "<<protein_id_threshold;
        if(getProteinProbability(minimal_protein_list_accessions_[i]) > protein_id_threshold) ++num; //std::cout << "--->";}
        // std::cout <<std::endl;
      }
    return num;
  }

  Int PSProteinInference::getNumberOfProtIdsPeptideRule(Int min_peptides,	std::map<String,std::set<String> >& prot_id_counter)
  {
    Int num = 0;
    for(Size i = 0; i < minimal_protein_list_accessions_.size();++i)
      {
        // std::cout << minimal_protein_list_accessions_[i] << " "<<getProteinProbability(minimal_protein_list_accessions_[i])
        //           << " "<<protein_id_threshold;
        if(prot_id_counter[minimal_protein_list_accessions_[i]].size() >= (Size)min_peptides) ++num; //std::cout << "--->";}
        // std::cout <<std::endl;
      }
    return num;
  }

  
}//namespace
