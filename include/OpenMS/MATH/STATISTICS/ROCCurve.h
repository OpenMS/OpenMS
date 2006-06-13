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
// $Id: ROCCurve.h,v 1.4 2006/03/28 16:19:59 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_MATH_STATISTICS_ROCCURVE_H
#define OPENMS_MATH_STATISTICS_ROCCURVE_H

#include <OpenMS/config.h>

#ifdef CGAL_DEF
  #include <CGAL/Cartesian.h>
  #include <CGAL/Polygon_2.h>
#endif


#include <list>
#include <vector>

namespace OpenMS
{
  /**
  	@brief ROCCurves show the tradeoff in sensitivity and specitivity for binary classifiers using different cutoff values <br>
  	
  	@todo add to Math namespace (Andreas)
  	
  	@ingroup Math
  */
  class ROCCurve
  {
  private:
  #ifdef CGAL_DEF
    typedef CGAL::Point_2<CGAL::Cartesian<double> > Point;
    typedef CGAL::Polygon_2<CGAL::Cartesian<double> > Polygon;
  #endif
    /**
    predicate for sort()
    */
    class simsortdec
    {
    public:
      bool operator()(const std::pair<double,bool>& a, const std::pair<double,bool>& b)
      {
        return b.first < a.first;
      }
    };
  public:
    /** @name constructors, destructor, assignment operator <br> */
    //@{
    ROCCurve();
    ~ROCCurve();
    ROCCurve( const ROCCurve& source );

    ROCCurve& operator=(const ROCCurve& source);
    //@}

    /** @brief insert score,type pair <br> */
    void insertPair(double score, bool clas);

    /** @brief Area Under Curve <br> */
    double AUC();

    /** @brief some points in the ROC Curve <br> */
    std::vector<std::pair<double,double> > curve(uint resolution = 10);

    double cutoffPos(double fraction = 0.95);

    double cutoffNeg(double fraction = 0.95);

  private:
    std::list<std::pair<double,bool> > score_clas_pairs_;
    uint pos_;
    uint neg_;
  };
}
#endif // OPENMS_MATH_STATISTICS_ROCCURVE_H
