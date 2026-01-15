// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <list>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
      @brief ROCCurves show the trade-off in sensitivity and specificity for binary classifiers using different cutoff values

    [This class is buggy and usage is discouraged!]

      @ingroup Math
    */
    class OPENMS_DLLAPI ROCCurve
    {
public:

      // @name Constructors and Destructors
      // @{
      /// default constructor
      ROCCurve();

      ///constructor with value, class pairs
      explicit ROCCurve(const std::vector<std::pair<double,bool>> & pairs);

      /// destructor
      virtual ~ROCCurve() = default;

      /// copy constructor
      ROCCurve(const ROCCurve & source);
      // @}

      // @name Operators
      // @{
      /// assignment operator
      ROCCurve & operator=(const ROCCurve & source);
      // @}

      // @name Accessors
      // @{
      /// insert score, type pair
      void insertPair(double score, bool clas);

      /// returns Area Under Curve
      double AUC();

      /// returns ROC-N score (e.g. ROC-50). Returns -1 if not enough false positives were found
      double rocN(Size N);

      /// some points in the ROC Curve
      std::vector<std::pair<double, double> > curve(UInt resolution = 10);

      /// returns the score at which you would need to set a cutoff to have fraction positives
      /// returns -1 if there are not enough positives
      double cutoffPos(double fraction = 0.95);

      /// returns the score at which you would need to set a cutoff to have fraction positives
      /// returns -1 if there are not enough positives
      double cutoffNeg(double fraction = 0.95);

      // @}

private:
      /// sorts data and caches if sorted
      inline void sort();

      /// counts global pos and neg
      inline void count();

      /// calculates area with trapezoidal rule
      /// @param x1,x2,y1,y2
      inline double trapezoidal_area(double x1, double x2, double y1, double y2);

      /// predicate for sort()
      class OPENMS_DLLAPI simsortdec
      {
public:

        bool operator()(const std::pair<double, bool> & a, const std::pair<double, bool> & b)
        {
          return b.first < a.first;
        }

      };


      std::vector<std::pair<double, bool> > score_clas_pairs_;

      UInt pos_{};

      UInt neg_{};

      bool sorted_{};
    };
  }
}
