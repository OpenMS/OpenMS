// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>

#include <OpenMS/CONCEPT/Factory.h>

#include <set>

namespace OpenMS
{

  BaseGroupFinder::BaseGroupFinder() :
    DefaultParamHandler("BaseGroupFinder")
  {
  }

  BaseGroupFinder::~BaseGroupFinder()
  {
  }

  void BaseGroupFinder::registerChildren()
  {
    Factory<BaseGroupFinder>::registerProduct(
      SimplePairFinder::getProductName(), &SimplePairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      LabeledPairFinder::getProductName(), &LabeledPairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      StablePairFinder::getProductName(), &StablePairFinder::create);
    Factory<BaseGroupFinder>::registerProduct(
      QTClusterFinder::getProductName(), &QTClusterFinder::create);
  }

  void BaseGroupFinder::checkIds_(const std::vector<ConsensusMap> & maps) const
  {
    std::set<Size> used_ids;
    for (Size i = 0; i < maps.size(); ++i)
    {
      const ConsensusMap & map = maps[i];
      for (ConsensusMap::FileDescriptions::const_iterator it = map.getFileDescriptions().begin(); it != map.getFileDescriptions().end(); ++it)
      {
        if (used_ids.find(it->first) != used_ids.end())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "file ids have to be unique");
        }
        else
        {
          used_ids.insert(it->first);
        }
      }
    }
  }

}
