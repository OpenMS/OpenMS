// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>

using namespace std;

namespace OpenMS
{

  MultiplexFilterResult::MultiplexFilterResult()
  {
  }

  void MultiplexFilterResult::addFilterResultPeak(double mz, double rt, vector<double> mz_shifts, std::vector<double> intensities, vector<MultiplexFilterResultRaw> raw_data_points)
  {
    MultiplexFilterResultPeak peak(mz, rt, mz_shifts, intensities, raw_data_points);
    result_.push_back(peak);
  }

  MultiplexFilterResultPeak MultiplexFilterResult::getFilterResultPeak(int i) const
  {
    return result_[i];
  }

  MultiplexFilterResultRaw MultiplexFilterResult::getFilterResultRaw(int i, int j) const
  {
    return result_[i].getFilterResultRaw(j);
  }

  double MultiplexFilterResult::getMZ(int i) const
  {
    return result_[i].getMZ();
  }

  vector<double> MultiplexFilterResult::getMZ() const
  {
    vector<double> mz;
    for (unsigned i = 0; i < result_.size(); ++i)
    {
      mz.push_back(result_[i].getMZ());
    }
    return mz;
  }

  double MultiplexFilterResult::getRT(int i) const
  {
    return result_[i].getRT();
  }

  vector<double> MultiplexFilterResult::getRT() const
  {
    vector<double> rt;
    for (unsigned i = 0; i < result_.size(); ++i)
    {
      rt.push_back(result_[i].getRT());
    }
    return rt;
  }

  int MultiplexFilterResult::size() const
  {
    return result_.size();
  }

}
