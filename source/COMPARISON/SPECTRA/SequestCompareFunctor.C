// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/SequestCompareFunctor.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/CLUSTERING/SpectrumGenerator.h>

#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  SequestCompareFunctor::SequestCompareFunctor()
    : CompareFunctor()
  {
		name_ = SequestCompareFunctor::getName();
    defaults_.setValue("isobaric", 1);
    defaults_.setValue("maxdeltaCN", 0.1);
    defaults_.setValue("tolerance", 0.1);
		param_ = defaults_;
    usebins_ = false;
  }
 
  SequestCompareFunctor::SequestCompareFunctor(const SequestCompareFunctor& source)
    : CompareFunctor(source)
  {
  }
  
  SequestCompareFunctor& SequestCompareFunctor::operator=(const SequestCompareFunctor& source)
  {
    CompareFunctor::operator=(source);
    return *this;
  }
  
  SequestCompareFunctor::~SequestCompareFunctor()
  {
  }

  //very ineffective since the Identificationes are taken from DB every time
  //DSpectrum<1>::identification or argument of operator() should be changed (mutable or not const)
  double SequestCompareFunctor::operator()(const ClusterSpectrum& csa, const ClusterSpectrum& csb) const
  {
    double deltaCN = (double)param_.getValue("maxdeltaCN");
    
    //e/peptides of different charge states show different (but somewhat similar) fragmentations
    if (csa.getParentionCharge() != csb.getParentionCharge() ) return 0;
    
    //when this was written there was just 1 identification per peaklist and no
    //way of identifying it
    const Identification* dbsa = 0;
    if (csa.getIdentification().size() ) dbsa = &*csa.getIdentification().begin();
    const Identification* dbsb = 0;
    if (csb.getIdentification().size() ) dbsb = &*csb.getIdentification().begin();
    
    if ( (! dbsa || !dbsb) || (!dbsa->getPeptideHits().size() || !dbsb->getPeptideHits().size()) )
    {
      return 0;
    }
    //compare the PeptideHits
    //either for == or alignment, or ...
    
    //simple comparison, if the two sequences match score = 1
    double scorea = 0;
    double scoreb = 0;
    //in case the peptidehits are not ordered
    for (vector<PeptideHit>::const_iterator lit = dbsa->getPeptideHits().begin(); lit != dbsa->getPeptideHits().end(); ++lit)
    {
      scorea = max(scorea, (double)lit->getScore());
    }
    for (vector<PeptideHit>::const_iterator lit = dbsb->getPeptideHits().begin(); lit != dbsb->getPeptideHits().end(); ++lit)
    {
      scoreb = max(scoreb, (double)lit->getScore());
    }
    vector<String> asequences;
    vector<String> bsequences;
    
    //look for top scoreing sequences
    for (vector<PeptideHit>::const_iterator lit = dbsa->getPeptideHits().begin(); lit != dbsa->getPeptideHits().end(); ++lit)
    {
      if ( (scorea - lit->getScore())/scorea < deltaCN)
      {
        asequences.push_back(lit->getSequence());
      }
    }
    for (vector<PeptideHit>::const_iterator lit = dbsb->getPeptideHits().begin(); lit != dbsb->getPeptideHits().end(); ++lit)
    {
      if ( (scoreb - lit->getScore())/scoreb < deltaCN)
      {
        bsequences.push_back(lit->getSequence());
      }
    }

    double score = 0;
    //if a pair of the topscoreing sequences matches
    for (uint i = 0; i < asequences.size();++i )
    {
      for (uint j = 0; j < bsequences.size(); ++j)
      {
        if ( fabs((double)param_.getValue("isobaric") - 1 ) < 1e-8 ) {
          if (asequences[i] == bsequences[j]) score = 1;
        }
        else {
          score = matchIsobaric(asequences[i],bsequences[j]);
        }
      }
    }
    return score;
  }

  double SequestCompareFunctor::matchIsobaric(const String& seq1, const String& seq2) const
  {
    double tolerance = (double)param_.getValue("tolerance");
    uint matches = 0;
    uint mismatches = 0;
    String::const_iterator sit1 = seq1.begin();
    String::const_iterator sit2 = seq2.begin();
    if (sit1 == seq1.end() || sit2 == seq2.end()) return 0;
    vector<double> aamasses1;
    vector<double> aamasses2;
    double mass1 = SpectrumGenerator::instance()->getPeptidemass(seq1);
    double mass2 = SpectrumGenerator::instance()->getPeptidemass(seq2);
    while ( sit1 != seq1.end() )
    {
      aamasses1.push_back(SpectrumGenerator::instance()->residuemass(sit1,seq1));
      ++sit1;
    }
    while ( sit2 != seq2.end() )
    {
      aamasses2.push_back(SpectrumGenerator::instance()->residuemass(sit2,seq2));
      ++sit2;
    }
    
    vector<double>::const_iterator cvit1 = aamasses1.begin();
    vector<double>::const_iterator cvit2 = aamasses2.begin();

    double bfragment1 = 1 + *cvit1++;
    double bfragment2 = 1 + *cvit2++;
    
    while ( cvit1 != aamasses1.end() && cvit2 != aamasses2.end() )
    {
      if ( abs(bfragment1 - bfragment2) < tolerance )
      {
        matches++;
        bfragment1 += *cvit1++;
        bfragment2 += *cvit2++;
      }
      else if ( bfragment1 > bfragment2 )
      {
        mismatches++;
        bfragment2 += *cvit2++; 
      }
      else if ( bfragment2 > bfragment1 )
      {
        mismatches++;
        bfragment1 += *cvit1++;
      }
    }

    while ( cvit1 != aamasses1.end() ) 
    {
      ++mismatches;
      ++cvit1;
    }
    while ( cvit2 != aamasses2.end() ){
      ++mismatches;
      ++cvit2;
    }

    vector<double>::reverse_iterator rcvit1 = aamasses1.rbegin();
    vector<double>::reverse_iterator rcvit2 = aamasses2.rbegin();
    
    double yfragment1 = mass1 - *rcvit1++ - 15;
    double yfragment2 = mass2 - *rcvit2++ - 15;
    
    while ( rcvit1 != aamasses1.rend() && rcvit2 != aamasses2.rend() )
    {
      if ( abs(yfragment1 - yfragment2) < tolerance )
      {
        matches++;
        yfragment1 += *rcvit1++;
        yfragment2 += *rcvit2++;
      }
      else if ( yfragment1 > yfragment2 )
      {
        mismatches++;
        yfragment2 += *rcvit2++; 
      }
      else if ( yfragment2 > yfragment1 )
      {
        mismatches++;
        yfragment1 += *rcvit1++;
      }
    }

    while ( rcvit1 != aamasses1.rend() ) 
    {
      ++mismatches;
      ++rcvit1;
    }
    while ( rcvit2 != aamasses2.rend() ) 
    {
      ++rcvit2;
      ++mismatches;
    }
    return (double) matches/(matches+mismatches);
  }
}
