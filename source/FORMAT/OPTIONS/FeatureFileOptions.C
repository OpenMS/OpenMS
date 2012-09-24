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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  FeatureFileOptions::FeatureFileOptions() :
    loadConvexhull_(true),
    loadSubordinates_(true),
    metadata_only_(false),
    has_rt_range_(false),
    has_mz_range_(false),
    has_intensity_range_(false)
  {
  }

  FeatureFileOptions::~FeatureFileOptions()
  {
  }

  void FeatureFileOptions::setLoadConvexHull(bool convex)
  {
    loadConvexhull_ = convex;
  }

  bool FeatureFileOptions::getLoadConvexHull() const
  {
    return loadConvexhull_;
  }

  void FeatureFileOptions::setLoadSubordinates(bool sub)
  {
    loadSubordinates_ = sub;
  }

  bool FeatureFileOptions::getLoadSubordinates() const
  {
    return loadSubordinates_;
  }

  void FeatureFileOptions::setMetadataOnly(bool only)
  {
    metadata_only_ = only;
  }

  bool FeatureFileOptions::getMetadataOnly() const
  {
    return metadata_only_;
  }

  void FeatureFileOptions::setRTRange(const DRange<1> & range)
  {
    rt_range_ = range;
    has_rt_range_ = true;
  }

  bool FeatureFileOptions::hasRTRange() const
  {
    return has_rt_range_;
  }

  const DRange<1> & FeatureFileOptions::getRTRange() const
  {
    return rt_range_;
  }

  void FeatureFileOptions::setMZRange(const DRange<1> & range)
  {
    mz_range_ = range;
    has_mz_range_ = true;
  }

  bool FeatureFileOptions::hasMZRange() const
  {
    return has_mz_range_;
  }

  const DRange<1> & FeatureFileOptions::getMZRange() const
  {
    return mz_range_;
  }

  void FeatureFileOptions::setIntensityRange(const DRange<1> & range)
  {
    intensity_range_ = range;
    has_intensity_range_ = true;
  }

  bool FeatureFileOptions::hasIntensityRange() const
  {
    return has_intensity_range_;
  }

  const DRange<1> & FeatureFileOptions::getIntensityRange() const
  {
    return intensity_range_;
  }

} // namespace OpenMS
