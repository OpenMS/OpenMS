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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{

  GridFeature::GridFeature(const BaseFeature & feature, Size map_index,
                           Size feature_index) :
    feature_(feature),
    map_index_(map_index),
    feature_index_(feature_index),
    annotations_()
  {
    const vector<PeptideIdentification> & peptides =
      feature.getPeptideIdentifications();
    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty())
        continue;                                    // shouldn't be the case
      annotations_.insert(pep_it->getHits()[0].getSequence());
    }
  }

  GridFeature::~GridFeature()
  {
  }

  const BaseFeature & GridFeature::getFeature() const
  {
    return feature_;
  }

  Size GridFeature::getMapIndex() const
  {
    return map_index_;
  }

  Size GridFeature::getFeatureIndex() const
  {
    return feature_index_;
  }

  Int GridFeature::getID() const
  {
    return (Int)feature_index_;
  }

  const set<AASequence> & GridFeature::getAnnotations() const
  {
    return annotations_;
  }

  double GridFeature::getRT() const
  {
    return feature_.getRT();
  }

  double GridFeature::getMZ() const
  {
    return feature_.getMZ();
  }

}
