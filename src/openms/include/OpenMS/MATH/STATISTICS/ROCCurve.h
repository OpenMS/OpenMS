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

      /// counts global pos and negs
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

      UInt pos_;

      UInt neg_;

      bool sorted_;
    };
  }
}
