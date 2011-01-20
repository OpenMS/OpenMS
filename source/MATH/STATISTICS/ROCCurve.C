// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/MATH/STATISTICS/ROCCurve.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	namespace Math
	{
	  ROCCurve::ROCCurve()
	    : score_clas_pairs_(),pos_(0),neg_(0)
	  {
	  }
	
	  ROCCurve::~ROCCurve()
	  {
	  }
	
	  ROCCurve::ROCCurve( const ROCCurve& source )
	    : score_clas_pairs_(source.score_clas_pairs_),pos_(source.pos_), neg_(source.neg_)
	  {
	  }
	
	  ROCCurve& ROCCurve::operator = (const ROCCurve& source)
	  {
			if (this != &source)
			{
	    	score_clas_pairs_ = source.score_clas_pairs_;
	    	pos_ = source.pos_;
	    	neg_ = source.neg_;
			}
	    return *this;
	  }
	
	  void ROCCurve::insertPair(double score, bool clas)
	  {
	    score_clas_pairs_.push_back(std::make_pair(score,clas));
	    if ( clas )
	    {
	      ++pos_;
	    }
	    else
	    {
	      ++neg_;
	    }
	  }
	
	  double ROCCurve::AUC()
	  {
      if (score_clas_pairs_.size()==0)
      {
	      cerr << "ROCCurve::AUC() : unsuitable dataset (no positives or no negatives)\n";
        return 0.5;
      }

	    score_clas_pairs_.sort(simsortdec());
	    // value that is not in score_clas_pairs_
	    double prevsim = score_clas_pairs_.begin()->first + 1;
	    UInt truePos = 0;
	    UInt falsePos = 0;
      std::vector<DPosition<2> > polygon;
	    for ( list<pair<double,bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit )
	    {
	      if ( fabs(cit->first - prevsim) > 1e-8 )
	      {
	        polygon.push_back(DPosition<2> ((double)falsePos/neg_,(double)truePos/pos_));
	      }
	      if ( cit->second )
	      {
	        ++truePos;
	      }
	      else
	      {
	        ++falsePos;
	      }
	    }
      polygon.push_back(DPosition<2>(1,1));
      std::sort(polygon.begin(), polygon.end());
      DPosition<2> last(0,0);
      DoubleReal area(0);
      for (std::vector<DPosition<2> >::const_iterator it=polygon.begin(); it!=polygon.end();++it)
      {
        area += (it->getX() - last.getX()) * (it->getY());
        last = *it;
      }
	    return area;
	  }
	
	  std::vector<std::pair<double,double > > ROCCurve::curve(UInt resolution)
	  {
	    score_clas_pairs_.sort(simsortdec());
	    vector<pair<double,double> > result;
	    UInt position = 0;
	    UInt truePos = 0;
	    UInt falsePos = 0;
	    for ( list<pair<double,bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit )
	    {
	      if ( cit->second )
	      {
	        ++truePos;
	      }
	      else
	      {
	        ++falsePos;
	      }
	      if ( ((double)++position/score_clas_pairs_.size())*resolution > result.size() )
	      {
	        result.push_back(make_pair((double)falsePos/neg_,(double)truePos/pos_));
	      }
	    }
	    return result;
	  }
	
	  /** 
	  \param fraction
	  \return cutoff for classifying <i>fraction</i> of the positives right <br> 
	  */
	  double ROCCurve::cutoffPos(double fraction)
	  {
	    score_clas_pairs_.sort(simsortdec());
	    UInt truePos = 0;
	    for ( list<pair<double,bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit )
	    {
	      if ( cit->second )
	      {
	        if ( (double)truePos++/pos_ > fraction )
	        {
	          return cit->first;
	        }
	      }
	    }
	    return -1;
	  }
	
	  /** 
	  \param fraction
	  \return cutoff for classifying <i>fraction</i> of the negatives right <br> 
	  */
	  double ROCCurve::cutoffNeg(double fraction)
	  {
	    score_clas_pairs_.sort(simsortdec());
	    UInt trueNeg = 0;
	    for ( list<pair<double,bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit )
	    {
	      if ( cit->second )
	      {
	        if ( (double)trueNeg++/neg_ > 1-fraction )
	        {
	          return cit->first;
	        }
	      }
	    }
	    return -1;
	  }
	}
}
