// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/EuclideanSimilarity.h>

namespace OpenMS
{
  EuclideanSimilarity::EuclideanSimilarity() :
    scale_(1)
  {
  }

  EuclideanSimilarity::EuclideanSimilarity(const EuclideanSimilarity & source) :
    scale_(source.scale_)
  {
  }

  EuclideanSimilarity::~EuclideanSimilarity()
  {
  }

  EuclideanSimilarity & EuclideanSimilarity::operator=(const EuclideanSimilarity & source)
  {
    if (this != &source)
    {
      scale_ = source.scale_;
    }
    return *this;
  }

  Real EuclideanSimilarity::operator()(const std::pair<Real, Real> & c) const
  {
    return operator()(c, c);
  }

  // calculates euclidean distance between two points
  Real EuclideanSimilarity::operator()(const std::pair<Real, Real> & a, const std::pair<Real, Real> & b) const
  {
    if (scale_ == 0)
    {
      //unapplicable scaling
      throw Exception::DivisionByZero(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    return 1 - (sqrt((a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second)) / scale_);
  }

  void EuclideanSimilarity::setScale(Real x)
  {
    scale_ = x;
  }

}
