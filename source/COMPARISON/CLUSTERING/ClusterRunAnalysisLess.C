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
// $Id: ClusterRunAnalysisLess.C,v 1.4 2006/03/28 10:07:05 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <iostream>
#include <cassert>

using namespace std;

namespace OpenMS
{
  ClusterExperiment::ClusterRunAnalysisLess::ClusterRunAnalysisLess(String cfig,String param)
    : configurablename_(cfig),paramname_(param),requirements_()
  {
  }
  
  ClusterExperiment::ClusterRunAnalysisLess::~ClusterRunAnalysisLess()
  {
  }
  
  void ClusterExperiment::ClusterRunAnalysisLess::setRequirement(String param, double value)
  {
    requirements_.insert(make_pair(param,value));
  }
  
  bool ClusterExperiment::ClusterRunAnalysisLess::operator()(const ClusterRun* ap, const ClusterRun* bp)
  {
    //check if both meet the requirements
    //if one doesnt, he is smaller
    //if both dont a is smaller
    
    //check if both contain the required AnalysisFunctor
    const ClusterRun& a = *ap;
    const ClusterRun& b = *bp;
    const FactoryProduct* ac = 0;
    const FactoryProduct* bc = 0;
    uint ai;
    uint bi;
    for (ai = 0; ai < a.size(); ++ai)
    {
      ac = a[ai].anafuncp();
      if (ac->getName() == configurablename_)
      {
        // we have found an AnalysisFunctor with the right name
        // no we look if the requirements are met
        for (map<String, double>::const_iterator cmit = requirements_.begin(); cmit != requirements_.end(); ++cmit)
        {
          for (Param::ConstIterator cmit2 = ac->getParam().begin(); cmit2 != ac->getParam().end(); ++cmit2)
          {
            if (cmit->first == cmit2->first)
            {
              if ( fabs(cmit->second - (double)cmit2->second ) > 1e-8 )
              {
                ac = 0;
              }
            }
          }
        }
        if (ac) break;
      }
    }
    if ( !ac )
    {
      cerr << "a doesnt fit\n";
      return false;
    }
    
    for (bi = 0; bi < b.size(); ++bi)
    {
      bc = b[bi].anafuncp();
      if (bc->getName() == configurablename_)
      {
        for (map<String, double>::const_iterator cmit = requirements_.begin(); cmit != requirements_.end(); ++cmit)
        {
          for (Param::ConstIterator cmit2 = bc->getParam().begin(); cmit2 != bc->getParam().end(); ++cmit2)
          {
            if (cmit->first == cmit2->first)
            {
              if ( fabs((double)cmit->second - (double)cmit2->second) > 1e-8 )
              {
                bc = 0;
              }
            }
          }
        }
        if (bc) break;
      }
    }
    if ( !bc )
    {
      cerr << "b doesnt fit\n";
      return false;
    }

    //comparison
    //both AnalysisFunctors take the same number of params
    map<String,double>::const_iterator cmita = a[ai].results().find(paramname_);
    map<String,double>::const_iterator cmitb = b[bi].results().find(paramname_);
    if (cmita == a[ai].results().end() || cmitb == b[bi].results().end())
    {
      return false;
    }
    return cmita->second < cmitb->second;
  }

}
