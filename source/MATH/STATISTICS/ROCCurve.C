// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
    ROCCurve::ROCCurve() :
      score_clas_pairs_(), pos_(0), neg_(0)
    {
    }

    ROCCurve::~ROCCurve()
    {
    }

    ROCCurve::ROCCurve(const ROCCurve & source) :
      score_clas_pairs_(source.score_clas_pairs_), pos_(source.pos_), neg_(source.neg_)
    {
    }

    ROCCurve & ROCCurve::operator=(const ROCCurve & source)
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
      score_clas_pairs_.push_back(std::make_pair(score, clas));
      if (clas)
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
      if (score_clas_pairs_.empty())
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
      for (list<pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
      {
        if (fabs(cit->first - prevsim) > 1e-8)
        {
          polygon.push_back(DPosition<2>((double)falsePos / neg_, (double)truePos / pos_));
        }
        if (cit->second)
        {
          ++truePos;
        }
        else
        {
          ++falsePos;
        }
      }
      polygon.push_back(DPosition<2>(1, 1));
      std::sort(polygon.begin(), polygon.end());
      DPosition<2> last(0, 0);
      DoubleReal area(0);
      for (std::vector<DPosition<2> >::const_iterator it = polygon.begin(); it != polygon.end(); ++it)
      {
        area += (it->getX() - last.getX()) * (it->getY());
        last = *it;
      }
      return area;
    }

    std::vector<std::pair<double, double> > ROCCurve::curve(UInt resolution)
    {
      score_clas_pairs_.sort(simsortdec());
      vector<pair<double, double> > result;
      UInt position = 0;
      UInt truePos = 0;
      UInt falsePos = 0;
      for (list<pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
      {
        if (cit->second)
        {
          ++truePos;
        }
        else
        {
          ++falsePos;
        }
        if (((double)++position / score_clas_pairs_.size()) * resolution > result.size())
        {
          result.push_back(make_pair((double)falsePos / neg_, (double)truePos / pos_));
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
      for (list<pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
      {
        if (cit->second)
        {
          if ((double)truePos++ / pos_ > fraction)
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
      for (list<pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
      {
        if (cit->second)
        {
          if ((double)trueNeg++ / neg_ > 1 - fraction)
          {
            return cit->first;
          }
        }
      }
      return -1;
    }

  }
}
