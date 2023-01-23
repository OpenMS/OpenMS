// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/FeatureOverlapFilter.h>

#include <Quadtree.h>
#include <Box.h>
#include <unordered_set>

namespace OpenMS
{
  void FeatureOverlapFilter::filter(FeatureMap& fmap)
  {

    fmap.sortByOverallQuality(true);
    const auto getBox = [](const Feature* f)
    {
        const auto& bb = f->getConvexHull().getBoundingBox();
        return quadtree::Box<float>(bb.minY(),bb.minX(),bb.maxY()-bb.minY(),bb.maxX()-bb.minX());
    };

    float minMZ = fmap.getMinMZ();
    float maxMZ = fmap.getMaxMZ();
    float minRT = fmap.getMinRT();
    float maxRT = fmap.getMaxRT();

    quadtree::Box<float> fullExp(minMZ-1,minRT-1,maxMZ-minMZ+2,maxRT-minRT+2);
    auto quadtree = quadtree::Quadtree<Feature*, decltype(getBox)>(fullExp, getBox);
    for (auto& f : fmap)
        quadtree.add(&f);

    std::unordered_set<Size> removed_uids;
    for (const auto& f : fmap)
    {
      if (removed_uids.count(f.getUniqueId()) == 0)
      {
        for (const auto& overlap : quadtree.query(getBox(&f)))
        {
          if (overlap != &f)
          {
            removed_uids.insert(overlap->getUniqueId());
            overlap->setOverallQuality(-1.); // used to filter
          }
        }
      }
    }

    const auto lowQuality = [](const Feature& f)
    {
        return f.getOverallQuality() < 0;
    };
    fmap.erase(std::remove_if(fmap.begin(), fmap.end(), lowQuality), fmap.end());

  }
}