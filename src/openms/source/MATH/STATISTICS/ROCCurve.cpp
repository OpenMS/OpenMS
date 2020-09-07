// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
//
#include <OpenMS/MATH/STATISTICS/ROCCurve.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

#include <cmath>
#include <numeric>

namespace OpenMS
{
  namespace Math
  {
    ROCCurve::ROCCurve() :
      score_clas_pairs_(), pos_(0), neg_(0)
    {
    }

    ROCCurve::ROCCurve(const ROCCurve & source) :
      score_clas_pairs_(source.score_clas_pairs_), pos_(source.pos_), neg_(source.neg_)
    {
    }

    ROCCurve::ROCCurve(const std::vector<std::pair<double,bool>> & pairs) :
        score_clas_pairs_(pairs)
    {
      pos_ = std::accumulate(score_clas_pairs_.begin(), score_clas_pairs_.end(), 0u,
                             [&](const UInt& x, const std::pair<double,bool>& y )
                              {
                                return x + y.second;
                              }
      );
      neg_ = static_cast<UInt>(score_clas_pairs_.size()) - pos_;
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
      score_clas_pairs_.emplace_back(std::make_pair(score, clas));
      if (clas)
      {
        ++pos_;
      }
      else
      {
        ++neg_;
      }
      sorted_ = false;
    }

    double ROCCurve::AUC()
    {
      if (score_clas_pairs_.empty())
      {
        std::cerr << "ROCCurve::AUC() : unsuitable dataset (no positives or no negatives)\n";
        return 0.5;
      }

      sort();
      double prevscore = -std::numeric_limits<double>::infinity();
      UInt prevTP = 0;
      UInt prevFP = 0;
      UInt truePos = 0;
      UInt falsePos = 0;
      double area = 0.0;
      for (auto const& pair : score_clas_pairs_)
      {
        //since it is sorted we do not need abs here.
        //TODO not sure if this makes much difference to an equality
        if ((pair.first - prevscore) > 1e-8)
        {
          area += trapezoidal_area(falsePos,prevFP,truePos,prevTP);
          prevscore = pair.first;
          prevTP = truePos;
          prevFP = falsePos;
        }
        pair.second ? ++truePos : ++falsePos;
      }
      area += trapezoidal_area(falsePos,prevFP,truePos,prevTP);
      // scale to unit square
      area /= truePos * falsePos;

      // update internals
      pos_ = truePos;
      neg_ = falsePos;

      return area;
    }

    double ROCCurve::rocN(Size N)
    {
      if (score_clas_pairs_.size() < N)
      {
        std::cerr << "ROCCurve::rocN() : unsuitable dataset (not enough false positives)\n";
        return -1;
      }

      sort();
      count();
      // value that is not in score_clas_pairs_
      double prevsim = score_clas_pairs_.begin()->first + 1;
      UInt truePos = 0;
      UInt falsePos = 0;
      std::vector<DPosition<2> > polygon;
      for (std::vector<std::pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end() && falsePos <= N; ++cit)
      {
        if (fabs(cit->first - prevsim) > 1e-8)
        {
          polygon.emplace_back(DPosition<2>((double)falsePos / neg_, (double)truePos / pos_));
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
      polygon.emplace_back(DPosition<2>(1, 1));
      std::sort(polygon.begin(), polygon.end());
      DPosition<2> last(0, 0);
      double area(0);
      for (std::vector<DPosition<2> >::const_iterator it = polygon.begin(); it != polygon.end(); ++it)
      {
        area += (it->getX() - last.getX()) * (it->getY());
        last = *it;
      }

      if (falsePos < N)
      {
        std::cerr << "ROCCurve::rocN() : unsuitable dataset (not enough false positives)\n";
        return -1;
      }
      return area;
    }

    std::vector<std::pair<double, double>> ROCCurve::curve(UInt resolution)
    {
      sort();
      count();
      std::vector<std::pair<double, double> > result;
      UInt position = 0;
      UInt truePos = 0;
      UInt falsePos = 0;
      for (auto const& pair : score_clas_pairs_)
      {
        pair.second ? ++truePos : ++falsePos;
        if (((double)++position / score_clas_pairs_.size()) * resolution > result.size())
        {
          result.emplace_back(std::make_pair((double)falsePos / neg_, (double)truePos / pos_));
        }
      }
      return result;
    }

    double ROCCurve::cutoffPos(double fraction)
    {
      sort();
      count();
      UInt truePos = 0;
      for (std::vector<std::pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
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

    double ROCCurve::cutoffNeg(double fraction)
    {
      sort();
      count();
      UInt trueNeg = 0;
      for (std::vector<std::pair<double, bool> >::const_iterator cit = score_clas_pairs_.begin(); cit != score_clas_pairs_.end(); ++cit)
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


    void ROCCurve::sort()
    {
      if (!sorted_)
      {
        std::sort(score_clas_pairs_.begin(),score_clas_pairs_.end(), simsortdec());
        sorted_ = true;
      }
    }

    void ROCCurve::count()
    {
      if (pos_ == 0 && neg_ == 0)
      {
        pos_ = std::accumulate(score_clas_pairs_.begin(), score_clas_pairs_.end(), 0u,
                               [&](const UInt& x, const std::pair<double,bool>& y )
                               {
                                 return x + y.second;
                               }
        );
        neg_ = static_cast<UInt>(score_clas_pairs_.size()) - pos_;
      }
    }

    double ROCCurve::trapezoidal_area(double x1, double x2, double y1, double y2)
    {
      double base = fabs(x1 - x2);
      double avgHeight = (y1+y2)/2.0;
      return base * avgHeight;
    }


  }
}
